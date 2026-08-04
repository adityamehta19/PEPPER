# MetalMotifExtractor

A general-purpose tool for extracting experimentally
verified metal-ion-binding structural motifs from RCSB Protein Data Bank
structures. Designed primarily for EF-hand Ca²⁺-binding proteins but
fully configurable for any metal-binding protein family with curated
UniProt Swiss-Prot domain annotations.

---

## Features

- **Any metal ion** — specify `CA`, `ZN`, `MG`, `MN`, `FE`, `FE2`, `FE3`,
  `CU`, `NI`, `CO`, `CD`, `HG`, `NA`, `K`, `LI`, `SR`, `BA`, lanthanides
  (`LA`, `GD`, `TB`, `EU`), or any other PDB HETATM residue name. Pass
  `--metal any` to extract all known metal types simultaneously.
- **Three input modes** — Excel spreadsheet, plain-text PDB ID list, or
  fully automatic discovery via the UniProt REST API.
- **Any protein family** — change one parameter (`uniprot_query`) to target
  zinc finger proteins, RING domains, rubredoxins, ferredoxins, or any
  other family with Swiss-Prot domain annotations.
- **Span-anchored window extraction** — guarantees every binding residue
  is always inside the extracted window, including the critical position-12
  glutamate in EF-hand loops that median-centred approaches miss.
- **Coordination analysis** — reports coordination number (donor atoms
  within cutoff, water excluded), geometry classification (tetrahedral,
  octahedral, pentagonal bipyramidal, etc.), τ₄/τ₅ indices, and all
  pairwise L–M–L angles.
- **Secondary structure assignment** — DSSP (primary), VMD/STRIDE
  (secondary, no extra install needed if VMD is available), or backbone
  dihedral fallback.
- **Multi-model extraction** — processes every model in multi-model PDB
  entries independently, preserving conformational diversity.
- **Quality filtering** — resolution, experimental method, and metal
  presence verification before download.
- **Parallel downloads** — multi-threaded structure retrieval.
- **Institutional network support** — configurable SSL bypass and explicit
  proxy settings for university/corporate networks.
- **Three output CSVs** — `structure_metadata.csv`, `motif_metadata.csv`
  (28 columns including coordination and secondary structure), and
  `secondary_structure.csv` (per-residue detail).
- **Full logging** — `logs/log.txt`.

---

## Installation

```bash
git clone <repository-url> MetalMotifExtractor
cd MetalMotifExtractor
pip install -r requirements.txt
```

Requires Python 3.11+.

### Optional dependencies for secondary structure

| Software | Purpose | Install |
|---|---|---|
| `mkdssp` | DSSP assignment (most accurate) | `sudo apt install dssp` or `conda install -c salilab dssp` |
| VMD | STRIDE assignment via headless batch mode | https://www.ks.uiuc.edu/Research/vmd/ |

If neither is installed, backbone dihedral angles are used as a fallback.

---

## Quick Start

### Automatic EF-hand discovery — all metals (recommended default)

Discovers all EF-hand proteins from UniProt Swiss-Prot and extracts
motifs for every metal ion found in each structure:

```bash
python metal_motif_extractor.py --metal any --resolution 2.0 \
    --distance 4 --motif_length 12 --mode uniprot
```

### Specific metal — Ca²⁺ only

```bash
python metal_motif_extractor.py --metal CA --resolution 2.0 \
    --distance 4 --motif_length 12 --mode uniprot
```

### Multiple metals in one run

```bash
python metal_motif_extractor.py --metal CA MG --resolution 2.0 \
    --distance 4 --motif_length 12 --mode uniprot
```

### Excel input (Mode A)

```bash
python metal_motif_extractor.py --excel proteins.xlsx --metal CA \
    --resolution 2.0 --distance 4 --motif_length 12
```

`proteins.xlsx` must contain a single column of PDB IDs (header optional):

| PDB_ID |
|--------|
| 1CLL   |
| 1CDL   |
| 3CLN   |

### Text file input (Mode B)

```bash
python metal_motif_extractor.py --textfile ids.txt --metal CA MG \
    --resolution 2.5 --distance 4 --motif_length 12
```

`ids.txt` — one PDB ID per line, `#` lines are comments:

```
# Calmodulin structures
1CLL
1CDL
3CLN
```

---

## Targeting Other Protein Families

Change `uniprot_query` in `config.yaml` or pass `--uniprot_query` on the
CLI to target any protein family with curated Swiss-Prot domain annotations.
All other pipeline stages — extraction, coordination analysis, secondary
structure, output — remain unchanged.

| Protein family | Metal | Query |
|---|---|---|
| EF-hand (default) | CA, MG | `(ft_domain:"EF-hand") AND (reviewed:true)` |
| C2H2 zinc finger | ZN | `(ft_domain:"Zinc finger, C2H2-type") AND (reviewed:true)` |
| RING finger (E3 ligase) | ZN | `(ft_domain:"RING-type zinc finger") AND (reviewed:true)` |
| LIM domain | ZN | `(ft_domain:"LIM zinc-binding") AND (reviewed:true)` |
| PHD finger | ZN | `(ft_domain:"PHD-type zinc finger") AND (reviewed:true)` |
| Rubredoxin | FE, FE2 | `(ft_domain:"Rubredoxin-type") AND (reviewed:true)` |
| 2Fe-2S ferredoxin | FE, FE2 | `(ft_domain:"2Fe-2S ferredoxin-type") AND (reviewed:true)` |
| EGF-like Ca-binding | CA | `(ft_domain:"EGF-like calcium-binding") AND (reviewed:true)` |
| FYVE zinc finger | ZN | `(ft_domain:"Zinc finger, FYVE-type") AND (reviewed:true)` |

**Via config.yaml:**
```yaml
uniprot_query: '(ft_domain:"RING-type zinc finger") AND (reviewed:true)'
metal: "ZN"
```

**Via CLI** (multi-word queries with commas are automatically reassembled):
```bash
python metal_motif_extractor.py --mode uniprot --metal ZN \
    --uniprot_query "Zinc finger, C2H2-type"
```

---

## UniProt Discovery Workflow (Mode C — `--mode uniprot`)

The `--mode uniprot` option automatically discovers protein structures
using the UniProt REST API. Note: despite the historical name `interpro`
appearing in some internal references, **InterPro is never contacted** —
all discovery calls go to `rest.uniprot.org`.

The workflow has three steps:

Check out our publication for the complete workflow.

## Command Line Interface

| Flag | Description |
|---|---|
| `--mode` | `excel`, `text`, or `uniprot` |
| `--excel FILE` | Excel file with PDB IDs (Mode A) |
| `--textfile FILE` | Text file with PDB IDs (Mode B) |
| `--uniprot_query QUERY` | UniProt query for `--mode uniprot`; multi-word queries are reassembled automatically |
| `--metal METAL [...]` | Metal code(s); space-separated. Use `any` for all known metals |
| `--distance` | Distance cutoff r_c (Å) for binding residue and coordination number calculation (default 4.0) |
| `--motif_length` | Requested window length in residues (default 15; minimum 12 for EF-hand) |
| `--resolution` | Maximum crystallographic resolution (Å) |
| `--include_nmr` | Accept NMR structures (no numeric resolution) |
| `--methods METHOD [...]` | Allowed experimental methods |
| `--biological_assembly` | Use biological assembly instead of asymmetric unit |
| `--workers` | Parallel download threads (default 8) |
| `--proxy URL` | HTTP proxy, e.g. `http://proxy.institution.ac.in:8080` |
| `--no_ssl_verify` | Disable SSL certificate verification (institutional networks) |
| `--config` | Path to YAML config file (default `config.yaml`) |
| `--output_dir`, `--pdb_dir`, `--motif_dir`, `--metadata_dir`, `--log_dir` | Override output paths |

Run `python metal_motif_extractor.py --help` for the full list.

---

## Configuration File

All parameters can be set persistently in `config.yaml`. CLI flags always
take precedence over config file values.

```yaml
# Metal ion(s): single name, list, or "any" for all 27 known metals
metal: "any"

# Input mode: excel, text, or uniprot
input_mode: "uniprot"

# UniProt query — change this to target a different protein family
uniprot_query: '(ft_domain:"EF-hand") AND (reviewed:true)'

max_resolution: 2.0          # Å — smaller = higher quality
include_nmr: false
allowed_methods:
  - "X-RAY DIFFRACTION"

distance_cutoff: 4.0         # Å — affects binding residues AND coordination number
motif_length: 12             # residues in extracted window (actual length reported in CSV)
use_biological_assembly: false
max_workers: 8

ssl_verify: false            # set false for institutional SSL-inspecting proxies
proxy: ""                    # e.g. "http://proxy.iitgn.ac.in:8080"
```

---

## Output Files

```
MetalMotifExtractor/
├── pdb_files/                        # Downloaded full PDB coordinate files
│   └── 1cll.pdb
├── motif_pdbs/                       # Extracted motif windows + coordinating ion
│   ├── 1CLL_A_CA_motif1.pdb
│   ├── 1CLL_A_CA_motif2.pdb
│   └── ...
├── metadata/
│   ├── structure_metadata.csv        # One row per downloaded structure (8 columns)
│   ├── motif_metadata.csv            # One row per extracted motif (28 columns)
│   └── secondary_structure.csv       # One row per residue per motif (14 columns)
├── logs/
│   └── log.txt                       # Full run log
└── outputs/                          # Reserved for downstream exports
```

### `structure_metadata.csv` — 8 columns

`PDB ID` · `Resolution` · `Experimental Method` · `Release Date` ·
`Organism` · `Protein Name` · `Chain IDs` · `Metal Ion Found`

The `Metal Ion Found` column reports which metals were actually extracted,
not just what was requested.

### `motif_metadata.csv` — 28 columns

**Identity:** Motif ID · PDB ID · Chain · Metal

**Quality:** Resolution · Experimental Method

**Binding site:** Binding Residues · Residue Numbers · Donor Atoms ·
Coordination Type · Distance Cutoff · Motif Length *(actual extracted
length, not necessarily the requested value)* · Sequence · Output Filename

**Coordination analysis:** Coordination Number *(donor atoms within
cutoff, water excluded)* · Coordination Geometry · tau4 · tau5 ·
Mean LML Angle (deg) · Min LML Angle (deg) · Max LML Angle (deg) ·
Largest Angle 1 (deg) · Largest Angle 2 (deg)

**Secondary structure:** SS Method *(DSSP / VMD_STRIDE / BioPDB_dihedral)* ·
SS Sequence *(per-residue 8-state codes, length = Motif Length)* ·
Helix Fraction · Strand Fraction · Coil Fraction

### `secondary_structure.csv` — 14 columns

One row per residue per motif. For N motifs of length L this file has
up to N × L rows.

`Motif ID` · `PDB ID` · `Chain` · `Metal` · `Position` *(1-based within
window)* · `Residue Name` · `Residue Number` · `One Letter` ·
`SS 8-state` · `SS 3-state` · `SS Full Name` · `Phi (deg)` ·
`Psi (deg)` · `SS Method`

---

## Software Architecture

The tool is organised as eight Python modules plus the CLI layer:

| File | Responsibility |
|---|---|
| `utils.py` | Config loading, logging, HTTP retry, SSL/proxy control, `ALL_METALS` constant |
| `downloader.py` | All three input modes; UniProt pagination; parallel PDB download |
| `metadata.py` | RCSB Data API queries; resolution/method filter; metal presence check |
| `motif.py` | Ion detection; KD-tree contact search; span-anchored window extraction; multi-model orchestration |
| `coordination.py` | Donor atom enumeration; L–M–L angles; τ₄/τ₅ indices; geometry classification |
| `secondary_structure.py` | DSSP → VMD/STRIDE → dihedral fallback pipeline |
| `metal_motif_extractor.py` | Main orchestration; CSV writers; summary report |
| `cli.py` | `argparse` command-line interface |
| `config.yaml` | All runtime parameters with documentation |
| `requirements.txt` | Python package dependencies |

---

## Algorithm Notes

### Span-anchored window extraction

The previous median-centred window algorithm failed for the canonical
EF-hand binding pattern. With binding residue indices {1, 3, 5, 12},
the median is 4, producing a window [−2:10] that excludes the critical
position-12 glutamate (the bidentate Glu providing two Ca²⁺ ligands).

The corrected algorithm anchors the window on the actual span of binding
residues. The span [i_min, i_max] is computed first; padding is then
distributed symmetrically outward to reach the requested motif length;
and clamping at chain termini shifts the window without shrinking it.
Every binding residue is guaranteed to be inside the window. When the
binding residue span exceeds the requested motif length (unusual, but
possible for multi-site ions), the window expands to cover the full span
rather than truncating any binding residue.

### Coordination number

Counts individual donor atoms within r_c of the metal ion, not residues.
A bidentate Glu (OE1 + OE2 both within r_c) contributes 2 to the count.
Water molecules are excluded regardless of distance. The coordination
number depends directly on the configured r_c — a smaller cutoff gives
a smaller coordination number for the same site.

### Multi-model extraction

When a PDB entry contains multiple models (NMR or multi-model X-ray),
each model is processed independently. Motifs from different models of
the same entry receive unique Motif IDs and are treated as independent
data points, preserving conformational diversity in the output dataset.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `CERTIFICATE_VERIFY_FAILED` | Institutional SSL-inspecting proxy | Add `ssl_verify: false` to `config.yaml` or use `--no_ssl_verify` |
| HTML redirect to campus login page | Network session expired (captive portal) | Open a browser, log in to the campus network portal, then re-run |
| Tool reaches external sites in browser but not in Python | Proxy not auto-detected by Python | Set `proxy: "http://your.proxy:port"` in `config.yaml` or use `--proxy` |
| `No PDB IDs were acquired` with `--mode uniprot` | UniProt blocked on network | Try after campus login; or use `--mode excel`/`--mode text` with a pre-built PDB list |
| `no ions found matching [...]` in log | Structure is apo or metal-depleted | Expected — tool correctly skips; re-run with `--metal any` to catch all metals |
| Motif length in CSV differs from `--motif_length` | Binding residue span exceeded requested length, or chain terminus | Correct behaviour — actual extracted length is reported |
| Position-12 Glu missing from EF-hand motif | Old version using median-centred window | Update to current version — span-anchored algorithm fixes this |
| `Failed to parse <id>` | Corrupted download | Delete the file in `pdb_files/` and re-run |
| SS Sequence all `C` | BioPDB dihedral fallback running on linear synthetic structure | Install VMD or mkdssp for accurate assignment on real structures |
| Slow run | Serial metadata API calls for large candidate sets | Increase `--workers`; RCSB rate limits apply so keep ≤ 8 |

---

## Supported Metal Ions

| Ion | Code | | Ion | Code |
|---|---|---|---|---|
| Calcium(II) | `CA` | | Lanthanum(III) | `LA` |
| Magnesium(II) | `MG` | | Gadolinium(III) | `GD` |
| Zinc(II) | `ZN` | | Terbium(III) | `TB` |
| Manganese(II) | `MN` | | Europium(III) | `EU` |
| Iron(II) | `FE2` | | Ytterbium(III) | `YB` |
| Iron(III) | `FE3` | | Samarium(III) | `SM` |
| Iron (unspecified) | `FE` | | Lutetium(III) | `LU` |
| Copper(II) | `CU` | | Neodymium(III) | `ND` |
| Copper(I) | `CU1` | | Cerium(III) | `CE` |
| Nickel(II) | `NI` | | Sodium(I) | `NA` |
| Cobalt(II) | `CO` | | Potassium(I) | `K` |
| Cadmium(II) | `CD` | | Lithium(I) | `LI` |
| Mercury(II) | `HG` | | Strontium(II) | `SR` |
| | | | Barium(II) | `BA` |

For complete iron coverage use `--metal FE FE2 FE3`.
Use `--metal any` to search for all 27 codes simultaneously.

---

## APIs Used

| API | Endpoint | Purpose |
|---|---|---|
| UniProt REST | `rest.uniprot.org/uniprotkb/search` | Protein family discovery (`--mode uniprot` only) |
| RCSB Data (entry) | `data.rcsb.org/rest/v1/core/entry/{id}` | Resolution, method, organism, protein name |
| RCSB Data (nonpolymer entity) | `data.rcsb.org/rest/v1/core/nonpolymer_entity/{id}/{eid}` | Metal presence verification |
| RCSB File Download | `files.rcsb.org/download/{id}.pdb` | Asymmetric unit coordinate file |
| RCSB Assembly Download | `files.rcsb.org/download/{id}-assembly1.cif` | Biological assembly (when enabled) |

---

## Citation

If you use this tool in published research, please cite:

> Berman, H.M. et al. "The Protein Data Bank." *Nucleic Acids Research*, 2000.
> The UniProt Consortium. "UniProt: the Universal Protein Knowledgebase in 2023." *Nucleic Acids Research*, 2023.
> Kabsch, W. & Sander, C. "Dictionary of protein secondary structure." *Biopolymers*, 1983. (if DSSP is used)
> Frishman, D. & Argos, P. "Knowledge-based protein secondary structure assignment." *Proteins*, 1995. (if VMD/STRIDE is used)

---

## License

Specify your project's license here prior to public release.

## Dependencies and Computational Requirements. 

The tool requires biopython ≥ 1.81, pandas ≥ 2.0, numpy ≥ 1.24, requests ≥ 2.31, PyYAML ≥ 6.0 and openpyxl ≥ 3.1. The optional mkdssp binary enables DSSP secondary-structure assignment. VMD activates STRIDE-based secondary structure assignment. No GPU or HPC infrastructure is necessary. Runtime scales linearly with the number of accepted structures.
