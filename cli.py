"""
cli.py
======
Command-line interface definition for MetalMotifExtractor.

Examples
--------
Automatic EF-hand discovery mode::

    python metal_motif_extractor.py --metal CA --resolution 2.0 \\
        --distance 3 --motif_length 15 --mode uniprot

Excel input mode::

    python metal_motif_extractor.py --excel proteins.xlsx --metal CA --distance 2.8

Text file input mode::

    python metal_motif_extractor.py --textfile ids.txt --metal MG \\
        --resolution 2.5 --distance 2.8

Multiple metals in one run (e.g. Ca2+ and Mg2+ sites from the same
structures)::

    python metal_motif_extractor.py --textfile ids.txt --metal CA MG \\
        --resolution 2.0 --distance 2.8
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
        "--distance", dest="distance_cutoff", type=float, default=None, metavar="ANGSTROM",
        help="Distance cut-off in Angstrom for protein donor atoms (side-chain O/N/S and "
             "backbone carbonyl O). Default: 2.8. Water and non-protein ligand cut-offs "
             "follow this value unless --water_cutoff / --ligand_cutoff are given.",
    )
    struct_group.add_argument(
        "--motif_length", type=int, default=None,
        help="Number of residues in the extracted motif window.",
    )
    struct_group.add_argument(
        "--biological_assembly", action="store_true", default=None,
        help="Download and use the biological assembly (assembly 1, mmCIF "
             "format) instead of the asymmetric unit.",
    )
    struct_group.add_argument(
        "--asu_format", choices=["pdb", "cif"], default=None,
        help="File format of the asymmetric unit (default pdb; 'cif' = mmCIF). "
             "Running the same entries in both formats should give identical "
             "motifs; use it to check format handling. Ignored with "
             "--biological_assembly, which is always mmCIF.",
    )
    struct_group.add_argument(
        "--metal_cutoff", nargs="+", default=None, metavar="METAL=VALUE",
        help="Per-metal protein-donor cutoff (Angstrom) that overrides --distance "
             "for that metal, e.g. --metal_cutoff ZN=2.8 CA=3.2. The cutoff used "
             "is recorded for every motif.",
    )
    struct_group.add_argument(
        "--water_cutoff", type=float, default=None, metavar="ANGSTROM",
        help="Maximum metal...O(water) distance for a water molecule to be "
             "counted as coordinated to the ion. Default: the same as --distance.",
    )
    struct_group.add_argument(
        "--no_water", action="store_true", default=False,
        help="Do not search for coordinated water (waters are included by "
             "default: counted in the CSV and written into the motif files).",
    )
    struct_group.add_argument(
        "--ligand_cutoff", type=float, default=None, metavar="ANGSTROM",
        help="Maximum metal...N/O/S distance for an atom of a non-protein, "
             "non-water ligand (nucleotide, cofactor, substrate, anion, "
             "nucleic-acid phosphate, ...) to be counted as coordinated. "
             "Default: the same as --distance.",
    )
    struct_group.add_argument(
        "--no_ligands", action="store_true", default=False,
        help="Do not search for non-protein ligands (included by default: "
             "counted in the CSV and written into the motif files).",
    )
    struct_group.add_argument(
        "--geometry_max_rmsd", type=float, default=None, metavar="ANGSTROM",
        help="Optional. Label a site 'Irregular' when even its best-fitting "
             "ideal polyhedron has an RMSD above this value (Angstrom). By "
             "default no threshold is applied: the best-fitting geometry is "
             "always reported, together with its RMSD and the runner-up.",
    )
    struct_group.add_argument(
        "--write_findgeo_input", action="store_true", default=None,
        help="Write, for every motif, a small PDB file (metal + exactly the "
             "donor atoms counted in the coordination number) to --findgeo_dir, "
             "for the FindGeo comparison (see findgeo_compare.py).",
    )
    struct_group.add_argument(
        "--findgeo_dir", default=None, metavar="DIR",
        help="Directory for the FindGeo input files (default findgeo_inputs).",
    )
    struct_group.add_argument(
        "--model", dest="model_policy", choices=["first", "all"], default=None,
        help="Which model(s) of a multi-model file to process. 'all' (default) "
             "processes every model, each independently, and adds the model "
             "number to motif IDs; 'first' uses the first model only.",
    )
    struct_group.add_argument(
        "--first_model_only", action="store_true", default=False,
        help="Shortcut for --model first: process only the first model of "
             "multi-model files.",
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

    # ---- Secondary structure ---------------------------------------------------
    ss_group = parser.add_argument_group("Secondary structure")
    ss_group.add_argument(
        "--dssp_exe", default=None, metavar="PATH",
        help="Full path to mkdssp (or dssp). Default: found on PATH or in a conda "
             "environment. When given, only this executable is tried.",
    )
    ss_group.add_argument(
        "--vmd_exe", default=None, metavar="PATH",
        help="Full path to VMD (used for the STRIDE method). Default: found on PATH "
             "or in C:\\Program Files\\VMD. When given, only this executable is tried.",
    )
    ss_group.add_argument(
        "--check_ss", action="store_true", default=False,
        help="Test which secondary-structure methods work on this machine "
             "(DSSP, VMD/STRIDE, dihedral fallback) on a built-in helix, then exit.",
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
        "asu_format": args.asu_format,
        "dssp_exe": args.dssp_exe,
        "vmd_exe": args.vmd_exe,
        "metal_cutoffs": args.metal_cutoff,
        "water_cutoff": args.water_cutoff,
        # --no_water sets include_water to False; otherwise leave the config value
        "include_water": False if args.no_water else None,
        "ligand_cutoff": args.ligand_cutoff,
        "include_ligands": False if args.no_ligands else None,
        "geometry_max_rmsd": args.geometry_max_rmsd,
        "write_findgeo_input": args.write_findgeo_input,
        "findgeo_dir": args.findgeo_dir,
        "model_policy": "first" if args.first_model_only else args.model_policy,
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
