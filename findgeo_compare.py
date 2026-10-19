#!/usr/bin/env python3
"""
findgeo_compare.py
==================
Command-line driver for the FindGeo comparison (see ``findgeo_tools.py``).

Comparison A - same donor atoms (tests the geometry classifier alone)
---------------------------------------------------------------------
1. Extract motifs and write the donor-atom files:
       python metal_motif_extractor.py ... --write_findgeo_input
   (creates ``findgeo_inputs/<motif_id>.pdb``: the metal and exactly the
   donor atoms counted for that motif)
2. Run FindGeo on them (threshold large enough to keep every atom):
       python findgeo_compare.py run-donors --inputs findgeo_inputs \
           --out findgeo_runs --findgeo "python /path/to/findgeo.py"
3. Compare:
       python findgeo_compare.py same-donors \
           --motifs metadata/motif_metadata.csv --runs findgeo_runs \
           --out comparison_same_donors

Comparison B - each program defines its own coordination sphere
---------------------------------------------------------------
1. Run FindGeo on the full structures (its own sphere; default 2.8 A):
       python findgeo_compare.py run-full --pdbs pdb_files \
           --out findgeo_full --findgeo "python /path/to/findgeo.py"
2. Compare:
       python findgeo_compare.py own-spheres \
           --motifs metadata/motif_metadata.csv --runs findgeo_full \
           --out comparison_own_spheres
   For a controlled comparison also run MetalMotifExtractor with the same
   distance definition (--distance 2.8 --water_cutoff 2.8 --ligand_cutoff 2.8)
   and compare again; the remaining differences are then due to which atom
   types count as donors, not to the cutoff.

FindGeo is not bundled.  The ``-p -w -o -t`` options used here are those of
the stand-alone ``findgeo.py`` documented on the FindGeo help page; check the
options of the version you installed (a Java version was released in 2024).
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

import findgeo_tools as ft


def _write_outputs(df: pd.DataFrame, out: Path, summaries: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "sites.csv", index=False)
    bad = df[~df["status"].isin(["agree"])]
    bad.to_csv(out / "disagreements.csv", index=False)
    summary_json = {}
    for label, (overall, per_cn, confusion) in summaries.items():
        summary_json[label] = overall
        per_cn.to_csv(out / f"by_cn_{label}.csv", index=False)
        confusion.to_csv(out / f"confusion_{label}.csv")
    with (out / "summary.json").open("w") as fh:
        json.dump(summary_json, fh, indent=2, default=str)
    for label, ov in summary_json.items():
        print(f"\n--- {label} ---")
        for k, v in ov.items():
            print(f"  {k:45s} {v}")
    print(f"\nWrote results to {out}/")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("run-donors", help="run FindGeo on the donor-atom PDB files")
    a.add_argument("--inputs", default="findgeo_inputs")
    a.add_argument("--out", default="findgeo_runs")
    a.add_argument("--findgeo", required=True, help='e.g. "python /path/findgeo.py"')
    a.add_argument("--threshold", type=float, default=4.5,
                   help="FindGeo -t (A); must be >= the largest donor distance (default 4.5)")
    a.add_argument("--workers", type=int, default=4)

    b = sub.add_parser("run-full", help="run FindGeo on the full structures")
    b.add_argument("--pdbs", default="pdb_files")
    b.add_argument("--out", default="findgeo_full")
    b.add_argument("--findgeo", required=True)
    b.add_argument("--threshold", type=float, default=2.8,
                   help="FindGeo -t (A); FindGeo's default is 2.8")
    b.add_argument("--workers", type=int, default=4)
    b.add_argument("--metal", default=None,
                   help="FindGeo -m: only this metal, e.g. CA (default: all metals)")

    c = sub.add_parser("same-donors", help="comparison A")
    d = sub.add_parser("own-spheres", help="comparison B")
    for p in (c, d):
        p.add_argument("--motifs", default="metadata/motif_metadata.csv")
        p.add_argument("--runs", required=True)
        p.add_argument("--out", required=True)

    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.cmd == "run-donors":
        res = ft.run_findgeo(args.findgeo, args.inputs, args.out, args.threshold, args.workers)
        print(f"{sum(res.values())}/{len(res)} FindGeo runs finished without error")
        return 0
    if args.cmd == "run-full":
        res = ft.run_findgeo(args.findgeo, args.pdbs, args.out, args.threshold,
                             args.workers, name_upper=True,
                             extra_args=["-m", args.metal] if args.metal else None)
        print(f"{sum(res.values())}/{len(res)} FindGeo runs finished without error")
        return 0

    motifs = pd.read_csv(args.motifs)
    sites = ft.collect_runs(args.runs)
    print(f"{len(motifs)} motifs, {len(sites)} FindGeo site results read from {args.runs}")
    out = Path(args.out)

    if args.cmd == "same-donors":
        df = ft.compare_same_donors(motifs, sites)
        n_cn = int((df.get("cn_match") == False).sum()) if "cn_match" in df else 0   # noqa: E712
        _write_outputs(df, out, {"all": ft.summarise_comparison(df)})
        if n_cn:
            print(f"\nWARNING: {n_cn} records have a different CN in FindGeo; check the "
                  f"threshold used for run-donors and the excluded-donor setting.")
        return 0

    df = ft.compare_own_spheres(motifs, sites)
    have = df[df["status"] != "no_findgeo_site"]
    extra = {}
    if len(have):
        extra = {
            "cn_equal_pct": round(100 * float(have["cn_match"].mean()), 1),
            "identical_sphere_pct": round(100 * float(have["same_sphere"].mean()), 1),
            "mean_donor_precision": round(float(have["donor_precision"].mean()), 3),
            "mean_donor_recall": round(float(have["donor_recall"].mean()), 3),
        }
    summaries = {
        "all_records": ft.summarise_comparison(df),
        "identical_spheres_only": ft.summarise_comparison(df, only_identical_spheres=True),
    }
    summaries["all_records"][0].update(extra)
    _write_outputs(df, out, summaries)
    return 0


if __name__ == "__main__":
    sys.exit(main())
