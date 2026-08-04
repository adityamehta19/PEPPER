"""
utils.py
========
Shared utilities for MetalMotifExtractor.

Provides
--------
- Configuration loading from YAML with CLI override merging
- Centralised logging initialisation (file + console handlers)
- Amino-acid three-letter → one-letter conversion
- Lightweight retry wrapper for HTTP requests

Biological rationale
--------------------
EF-hand motifs are helix–loop–helix calcium-binding structural elements.
The loop portion, typically 12 residues long, coordinates Ca²⁺ through
oxygen-containing side chains at positions 1, 3, 5, 7, 9 and 12 of the loop.
A 15-residue extraction window (3 flanking residues on each side) captures
enough secondary-structure context for downstream sequence analysis.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import requests
import yaml

# ---------------------------------------------------------------------------
# Three-to-one letter amino-acid mapping
# ---------------------------------------------------------------------------
AA3TO1: dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # selenocysteine / selenomethionine variants common in X-ray structures
    "SEC": "U", "MSE": "M",
}

# ---------------------------------------------------------------------------
# All biologically relevant metal ion PDB HETATM residue names.
# Used when --metal any is specified to search for every metal in a structure.
# ---------------------------------------------------------------------------
ALL_METALS: list[str] = [
    "CA", "MG", "ZN", "MN", "FE", "FE2", "FE3",
    "CU", "CU1", "NI", "CO", "CD", "HG",
    "NA", "K", "LI", "SR", "BA",
    "LA", "GD", "TB", "EU", "YB", "SM", "LU", "ND", "CE",
]


def residue_one_letter(resname: str) -> str:
    """Return the one-letter code for a three-letter residue name.

    Parameters
    ----------
    resname : str
        Three-letter residue name (e.g. ``"ALA"``).

    Returns
    -------
    str
        Single character (``"X"`` for unknown residues).
    """
    return AA3TO1.get(resname.strip().upper(), "X")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_DEFAULT_CONFIG: dict[str, Any] = {
    "metal": "CA",
    "input_mode": "uniprot",
    "input_file": "input.xlsx",
    "max_resolution": 2.0,
    "include_nmr": False,
    "allowed_methods": ["X-RAY DIFFRACTION"],
    "distance_cutoff": 4.0,
    "motif_length": 15,
    "use_biological_assembly": False,
    "max_workers": 8,
    "ssl_verify": True,
    "proxy": "",
    "any_metal": False,
    "uniprot_query": '(ft_domain:"EF-hand") AND (reviewed:true)',
    "output_dir": "outputs",
    "pdb_dir": "pdb_files",
    "motif_dir": "motif_pdbs",
    "metadata_dir": "metadata",
    "log_dir": "logs",
}


def load_config(config_path: str | Path = "config.yaml",
                overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    """Load configuration from a YAML file and apply optional CLI overrides.

    Parameters
    ----------
    config_path : str or Path
        Path to the YAML configuration file.
    overrides : dict, optional
        Key-value pairs that override values from the file (typically
        parsed from CLI arguments).  ``None`` values inside *overrides*
        are ignored so that unset CLI flags do not silently wipe YAML
        values.

    Returns
    -------
    dict
        Merged configuration dictionary.
    """
    cfg: dict[str, Any] = dict(_DEFAULT_CONFIG)

    path = Path(config_path)
    if path.exists():
        with path.open("r", encoding="utf-8") as fh:
            file_cfg = yaml.safe_load(fh) or {}
        cfg.update(file_cfg)
    else:
        logging.getLogger(__name__).warning(
            "Config file %s not found; using built-in defaults.", config_path
        )

    if overrides:
        cfg.update({k: v for k, v in overrides.items() if v is not None})

    # Normalise metal to a list of upper-case residue names.
    # Handles all of:
    #   metal: "CA"                   → ['CA']   (single string in YAML)
    #   metal: ["CA", "MG"]           → ['CA','MG']  (proper YAML list)
    #   metal: "['CA', 'MG']"         → ['CA','MG']  (Python-list-as-string,
    #                                                  common copy-paste mistake)
    #   CLI --metal CA MG             → ['CA','MG']  (list from argparse nargs="+")
    #   CLI --metal CA (old parser)   → ['CA']       (single string from argparse)
    raw_metal = cfg["metal"]

    if isinstance(raw_metal, list):
        # Already a proper list (from argparse nargs="+" or correct YAML)
        metals = [str(m).strip().upper() for m in raw_metal if str(m).strip()]

    elif isinstance(raw_metal, str):
        s = raw_metal.strip()
        if s.startswith("[") and s.endswith("]"):
            # Looks like a Python list literal written as a YAML string.
            # Parse it safely: strip brackets, split on commas, remove quotes.
            import re
            tokens = re.findall(r"[A-Za-z0-9]+", s[1:-1])
            metals = [t.upper() for t in tokens if t]
            if metals:
                logging.getLogger(__name__).warning(
                    "config.yaml has metal: %r — looks like a Python list "
                    "written as a quoted string. Parsed as %s. "
                    "Fix your config.yaml to use proper YAML list syntax: "
                    "metal: %s",
                    raw_metal, metals,
                    str(metals).replace("'", '"'),
                )
            else:
                logging.getLogger(__name__).error(
                    "Could not parse metal value %r — using default 'CA'.", raw_metal
                )
                metals = ["CA"]
        else:
            # Plain single metal name string e.g. "CA" or "ZN"
            metals = [s.upper()]
    else:
        metals = [str(raw_metal).strip().upper()]

    if not metals:
        raise ValueError(
            f"'metal' resolved to an empty list from raw value {raw_metal!r}. "
            "Provide at least one metal residue name, e.g. --metal CA"
        )

    # Special value "ANY" — search for all known biologically relevant metals.
    # Skips the metal presence filter at discovery time so every EF-hand
    # structure is accepted regardless of which metal it contains in that
    # specific crystal form.
    if metals == ["ANY"]:
        cfg["any_metal"] = True
        cfg["metal"] = ALL_METALS
        logging.getLogger(__name__).info(
            "--metal any: will search for all %d known metal types. "
            "Metal presence filter at discovery stage is disabled — "
            "all structures from the protein family query are accepted.",
            len(ALL_METALS),
        )
    else:
        cfg["any_metal"] = False
        cfg["metal"] = metals

    # Apply SSL and proxy settings immediately so all subsequent
    # http_get calls respect them.
    configure_ssl(bool(cfg.get("ssl_verify", True)))
    configure_proxy(str(cfg.get("proxy", "")))

    return cfg


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(log_dir: str | Path, level: int = logging.INFO) -> logging.Logger:
    """Initialise root logger with console and rotating file handlers.

    Parameters
    ----------
    log_dir : str or Path
        Directory where ``log.txt`` is written.
    level : int
        Logging level (default ``logging.INFO``).

    Returns
    -------
    logging.Logger
        Configured root logger.
    """
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "log.txt"

    fmt = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
    date_fmt = "%Y-%m-%d %H:%M:%S"

    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=date_fmt,
        handlers=[
            logging.FileHandler(log_file, mode="a", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger()


# Module-level SSL verification flag.
_SSL_VERIFY: bool = True

# Module-level proxy setting.
# Set via configure_proxy() when config is loaded.
# Format: "http://host:port" or "" for no proxy.
_PROXIES: dict | None = None


def configure_ssl(verify: bool) -> None:
    """Set the global SSL verification flag for all HTTP requests.

    Parameters
    ----------
    verify : bool
        If False, SSL certificate verification is disabled for all
        requests made through :func:`http_get`.  This is necessary
        when running behind an institutional proxy or VPN that performs
        SSL inspection (man-in-the-middle) and presents its own
        certificate instead of the remote server's certificate.

    Security note
    -------------
    Disabling SSL verification removes protection against genuine
    man-in-the-middle attacks on public networks.  Only set this to
    False on trusted institutional networks where SSL inspection is
    performed by your IT infrastructure.  Never disable it on
    untrusted or public Wi-Fi.
    """
    global _SSL_VERIFY
    _SSL_VERIFY = verify
    if not verify:
        # Suppress the urllib3 InsecureRequestWarning that fires on
        # every request when verify=False, since the user has already
        # acknowledged the implication by setting ssl_verify: false.
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
        logging.getLogger(__name__).warning(
            "SSL certificate verification is DISABLED. "
            "All HTTPS connections will proceed without verifying the "
            "server certificate. This is appropriate behind an SSL-inspecting "
            "institutional proxy but reduces security on untrusted networks."
        )


def configure_proxy(proxy: str) -> None:
    """Set the global proxy for all HTTP requests.

    Parameters
    ----------
    proxy : str
        Proxy URL, e.g. ``"http://fwg.iitgn.ac.in:80"`` or ``""`` for none.

    Notes
    -----
    Python's ``requests`` library does not automatically read the system
    proxy on all platforms.  When a browser can reach external sites but
    command-line Python cannot, the proxy must be set explicitly here.
    Set ``proxy`` in ``config.yaml`` to your institution's proxy address.
    """
    global _PROXIES
    log = logging.getLogger(__name__)
    if proxy:
        _PROXIES = {"http": proxy, "https": proxy}
        log.info("HTTP proxy configured: %s", proxy)
    else:
        _PROXIES = None


def http_get(url: str,
             params: dict | None = None,
             timeout: int = 30,
             retries: int = 3,
             backoff: float = 2.0,
             **kwargs: Any) -> requests.Response:
    """Perform a GET request with automatic exponential-backoff retries.

    SSL verification is controlled by the global ``_SSL_VERIFY`` flag,
    which is set from ``ssl_verify`` in ``config.yaml``.  If your
    institution uses an SSL-inspecting proxy (common in universities,
    hospitals, and corporate networks), set ``ssl_verify: false`` in
    ``config.yaml`` to bypass certificate verification.

    Parameters
    ----------
    url : str
        Request URL.
    params : dict, optional
        Query parameters.
    timeout : int
        Per-attempt timeout in seconds.
    retries : int
        Maximum number of attempts.
    backoff : float
        Multiplier for successive wait intervals (starts at 1 s).
    **kwargs
        Forwarded to :func:`requests.get`.

    Returns
    -------
    requests.Response
        Successful response object.

    Raises
    ------
    requests.HTTPError
        If all retries are exhausted.
    """
    log = logging.getLogger(__name__)
    wait = 1.0
    last_exc: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(
                url, params=params, timeout=timeout,
                verify=_SSL_VERIFY, proxies=_PROXIES, **kwargs
            )
            resp.raise_for_status()
            return resp
        except requests.exceptions.SSLError as exc:
            # Give a clear, actionable diagnosis on SSL failures instead
            # of the generic urllib3 message which is hard to interpret.
            last_exc = exc
            if attempt == 1:
                log.error(
                    "SSL certificate verification failed for %s.\n"
                    "\n"
                    "  This usually means your network (institutional VPN,\n"
                    "  university/corporate proxy, or firewall) is intercepting\n"
                    "  HTTPS traffic and presenting its own certificate instead\n"
                    "  of the remote server's certificate.\n"
                    "\n"
                    "  To fix this, add the following line to config.yaml:\n"
                    "      ssl_verify: false\n"
                    "\n"
                    "  Or pass --no_ssl_verify on the command line.\n"
                    "\n"
                    "  If the problem persists, try running:\n"
                    "      python -c \"import ssl; print(ssl.get_default_verify_paths())\"\n"
                    "  to find your system's certificate store location, and\n"
                    "  ask your IT department to add the proxy CA certificate.\n",
                    url,
                )
            # Do not retry SSL errors — they will not resolve on retry
            raise requests.HTTPError(
                f"SSL certificate error for {url}. "
                "Set ssl_verify: false in config.yaml if behind an "
                "SSL-inspecting proxy."
            ) from exc
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < retries:
                log.warning(
                    "Request failed (attempt %d/%d): %s — retrying in %.1fs",
                    attempt, retries, exc, wait,
                )
                time.sleep(wait)
                wait *= backoff
            else:
                log.error("All %d attempts failed for %s: %s", retries, url, exc)

    raise requests.HTTPError(f"Failed after {retries} attempts: {url}") from last_exc


def safe_json(resp: requests.Response, context: str = "") -> dict | list | None:
    """Parse a response as JSON with clear diagnostics on failure.

    Parameters
    ----------
    resp : requests.Response
        HTTP response object.
    context : str
        Human-readable description of what this request was for,
        used in log messages.

    Returns
    -------
    dict or list or None
        Parsed JSON, or None if the response was empty or not valid JSON.
    """
    log = logging.getLogger(__name__)

    if not resp.text or not resp.text.strip():
        log.warning(
            "Empty response body%s (status %d, url=%s).\n"
            "  This can happen when an SSL-inspecting proxy intercepts the\n"
            "  request and returns an empty or HTML error page instead of\n"
            "  the real API response.  The tool will continue where possible.",
            f" for {context}" if context else "",
            resp.status_code,
            resp.url,
        )
        return None

    # Check for HTML response (proxy error page instead of JSON)
    stripped = resp.text.strip()
    if stripped.startswith("<"):
        log.warning(
            "Received HTML instead of JSON%s (status %d, url=%s).\n"
            "  First 200 chars: %r\n"
            "  This usually means an institutional proxy or firewall is\n"
            "  intercepting the request and returning an error page.\n"
            "  Ensure ssl_verify: false is set if behind a proxy, and\n"
            "  that the target domain is not blocked by your network.",
            f" for {context}" if context else "",
            resp.status_code,
            resp.url,
            stripped[:200],
        )
        return None

    try:
        return resp.json()
    except Exception as exc:
        log.warning(
            "JSON decode error%s (status %d, url=%s): %s.\n"
            "  First 200 chars of response: %r",
            f" for {context}" if context else "",
            resp.status_code,
            resp.url,
            exc,
            stripped[:200],
        )
        return None


def ensure_dirs(*paths: str | Path) -> None:
    """Create directories (including parents) if they do not exist.

    Parameters
    ----------
    *paths
        One or more directory paths.
    """
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)
