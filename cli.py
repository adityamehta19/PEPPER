"""
cli.py
======
Command-line interface definition for MetalMotifExtractor.

Examples
--------
Automatic EF-hand discovery mode::

    python metal_motif_extractor.py --metal CA --resolution 2.0 \\
        --distance 4 --motif_length 15 --mode uniprot

Excel input mode::

    python metal_motif_extractor.py --excel proteins.xlsx --metal CA

Text file input mode::

    python metal_motif_extractor.py --textfile ids.txt --metal MG \\
        --resolution 2.5

Multiple metals in one run (e.g. Ca2+ and Mg2+ sites from the same
structures)::

    python metal_motif_extractor.py --textfile ids.txt --metal CA MG \\
        --resolution 2.0
"""

from __future__ import annotations

import argparse
from typing import Any


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for MetalMotifExtractor.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser instance.
    """
    parser = argparse.ArgumentParser(
        prog="metal_motif_extractor.py",
        description=(
            "MetalMotifExtractor — extract metal-ion-binding structural "
            "motifs (e.g. EF-hand Ca2+ loops) from RCSB PDB structures."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ---- Input mode selection ----------------------------------------------
    mode_group = parser.add_argument_group("Input mode")
    mode_group.add_argument(
        "--mode",
        choices=["excel", "text", "uniprot"],
        default=None,
        help=(
            "Input mode. "
            "'excel' reads PDB IDs from an Excel file. "
            "'text' reads PDB IDs from a plain text file. "
            "'uniprot' discovers protein structures automatically via the "
            "UniProt REST API using the configured uniprot_query. "
            "Inferred from --excel/--textfile if not given explicitly."
        ),
    )
    mode_group.add_argument(
        "--excel", dest="excel_file", metavar="FILE",
        help="Path to Excel (.xlsx) file containing a single column of PDB IDs.",
    )
    mode_group.add_argument(
        "--textfile", dest="text_file", metavar="FILE",
        help="Path to plain-text file containing one PDB ID per line.",
    )
    mode_group.add_argument(
        "--uniprot_query", nargs="+", default=None, metavar="QUERY",
        help=(
            "UniProt REST API query string for --mode uniprot. "
            "Controls which protein family is discovered. "
            'Default (EF-hand): \'(ft_domain:"EF-hand") AND (reviewed:true)\'. '
            'C2H2 zinc finger: \'(ft_domain:"Zinc finger, C2H2-type") AND (reviewed:true)\'. '
            'RING finger: \'(ft_domain:"RING-type zinc finger") AND (reviewed:true)\'. '
            "See config.yaml for more examples. "
            "On Windows, if the query contains spaces or commas, set it in "
            "config.yaml instead to avoid shell quoting issues."
        ),
    )

    # ---- Metal / structural parameters -------------------------------------
    struct_group = parser.add_argument_group("Structural parameters")
    struct_group.add_argument(
        "--metal", type=str, nargs="+", default=None,
        help=(
            'Metal ion residue name(s) (e.g. "CA", "ZN", "MG", "FE2"). '
            'Pass multiple values to extract several metals in one run: '
            '--metal CA MG. '
            'Pass "any" to extract ALL known metal types from every structure '
            'in the protein family (skips per-structure metal presence filter): '
            '--metal any. '
            'Use "any" with --mode uniprot to get all metals across all '
            'EF-hand crystal forms.'
        ),
    )
    struct_group.add_argument(
        "--distance", dest="distance_cutoff", type=float, default=None,
        help="Distance cutoff (Angstrom) for identifying binding residues.",
    )
    struct_group.add_argument(
        "--motif_length", type=int, default=None,
        help="Number of residues in the extracted motif window.",
    )
    struct_group.add_argument(
        "--biological_assembly", action="store_true", default=None,
        help="Download and use the biological assembly instead of the "
             "asymmetric unit.",
    )

    # ---- Filtering ----------------------------------------------------------
    filt_group = parser.add_argument_group("Structure filtering")
    filt_group.add_argument(
        "--resolution", dest="max_resolution", type=float, default=None,
        help="Maximum acceptable resolution (Angstrom).",
    )
    filt_group.add_argument(
        "--include_nmr", action="store_true", default=None,
        help="Include NMR structures (which lack numeric resolution).",
    )
    filt_group.add_argument(
        "--methods", dest="allowed_methods", nargs="+", default=None,
        metavar="METHOD",
        help='Allowed experimental methods, e.g. --methods "X-RAY DIFFRACTION" '
             '"ELECTRON MICROSCOPY".',
    )

    # ---- Performance ----------------------------------------------------------
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument(
        "--workers", dest="max_workers", type=int, default=None,
        help="Number of parallel download threads.",
    )
    perf_group.add_argument(
        "--proxy", default=None, metavar="URL",
        help=(
            "HTTP proxy URL, e.g. http://fwg.iitgn.ac.in:80. "
            "Use when Python cannot reach external APIs but your browser can "
            "(proxy not auto-detected on some institutional networks). "
            "Equivalent to setting proxy: in config.yaml."
        ),
    )
    perf_group.add_argument(
        "--no_ssl_verify", action="store_true", default=False,
        help=(
            "Disable SSL certificate verification. Use this if you are "
            "behind an institutional proxy, university VPN, or corporate "
            "firewall that performs SSL inspection and causes "
            "CERTIFICATE_VERIFY_FAILED errors. Equivalent to setting "
            "ssl_verify: false in config.yaml."
        ),
    )

    # ---- Paths ----------------------------------------------------------------
    path_group = parser.add_argument_group("Paths")
    path_group.add_argument(
        "--config", default="config.yaml",
        help="Path to YAML configuration file.",
    )
    path_group.add_argument(
        "--output_dir", default=None,
        help="Root output directory.",
    )
    path_group.add_argument(
        "--pdb_dir", default=None,
        help="Directory for downloaded PDB files.",
    )
    path_group.add_argument(
        "--motif_dir", default=None,
        help="Directory for extracted motif PDB files.",
    )
    path_group.add_argument(
        "--metadata_dir", default=None,
        help="Directory for metadata CSV files.",
    )
    path_group.add_argument(
        "--log_dir", default=None,
        help="Directory for log files.",
    )

    return parser


def parse_args_to_overrides(args: argparse.Namespace) -> dict[str, Any]:
    """Convert parsed CLI namespace into a config-override dictionary.

    Parameters
    ----------
    args : argparse.Namespace
        Result of ``parser.parse_args()``.

    Returns
    -------
    dict
        Dictionary suitable for passing to ``utils.load_config(overrides=...)``.
        Keys with value ``None`` are excluded so that unset flags don't
        override YAML defaults.
    """
    # Infer input_mode if not explicitly given
    mode = args.mode
    input_file = None
    if mode is None:
        if args.excel_file:
            mode = "excel"
            input_file = args.excel_file
        elif args.text_file:
            mode = "text"
            input_file = args.text_file
    else:
        if mode == "excel":
            input_file = args.excel_file
        elif mode == "text":
            input_file = args.text_file

    overrides: dict[str, Any] = {
        "input_mode": mode,
        "input_file": input_file,
        # --uniprot_query uses nargs='+' so the shell can split on commas/spaces
        # on Windows; rejoin all tokens into one string here.
        "uniprot_query": " ".join(args.uniprot_query) if args.uniprot_query else None,
        "metal": args.metal,
        "distance_cutoff": args.distance_cutoff,
        "motif_length": args.motif_length,
        "use_biological_assembly": args.biological_assembly,
        "max_resolution": args.max_resolution,
        "include_nmr": args.include_nmr,
        "allowed_methods": args.allowed_methods,
        "max_workers": args.max_workers,
        "proxy": args.proxy,
        # --no_ssl_verify sets ssl_verify to False
        "ssl_verify": False if args.no_ssl_verify else None,
        "output_dir": args.output_dir,
        "pdb_dir": args.pdb_dir,
        "motif_dir": args.motif_dir,
        "metadata_dir": args.metadata_dir,
        "log_dir": args.log_dir,
    }
    return {k: v for k, v in overrides.items() if v is not None}
