# PEPPER - PEPtide Pocket Extractor of Residues around metals

Extracts metal-ion-binding structural motifs (e.g. EF-hand Ca²⁺ loops) from
RCSB PDB structures: discovers candidate structures, filters them by
experimental quality, downloads coordinate files, and extracts fixed-length
sequence windows around metal coordination sites, together with the water
molecules coordinated to each ion.  Output: per-motif PDB files and
analysis-ready CSV tables.

> **Terminology.** The extracted object is a *sequence window* around the
> protein ligands of one ion (see *Algorithm notes*). It is not necessarily
> the complete coordination site.

---

## Features

- **Any metal ion** (`CA`, `ZN`, `MG`, `FE`, … or `any`).
- **Coordinated water**: waters whose oxygen lies within a configurable
  cutoff of the ion are counted, listed with their distances, and written
  (with coordinates) into the motif PDB file.
- **Coordinated non-protein ligands** (nucleotides, cofactors, substrates,
  anions, nucleic-acid phosphates, ...): their N/O/S atoms within a cutoff
  of the ion are counted and listed, and the whole ligand residue is
  written into the motif PDB file.
- **Three input modes**: Excel, plain-text ID list, or UniProt discovery.
- **Model-aware**: each model of a multi-model file is processed on its own;
  an ion is only ever matched to residues and waters of its own model.
- **PDB and mmCIF input**, including biological assemblies.
- **Parallel metadata retrieval and downloads.**
- **Auditable accounting**: every filtering stage has its own counter;
  `run_accounting.json` and `structure_outcomes.csv` explain every number.
- **Secondary structure** per residue (DSSP → VMD/STRIDE → dihedral fallback).

## Installation

```bash
pip install -r requirements.txt
# SciPy is required (geometry fitting)
# optional, for secondary structure: mkdssp (conda install -c conda-forge dssp,
# or sudo apt install dssp) and/or VMD; check with: python metal_motif_extractor.py --check_ss
```

Python 3.10+ is recommended.

### Secondary structure programs

`python metal_motif_extractor.py --check_ss` runs every method on a test helix
and says which one the tool will use, or why a program it found does not work.

- **mkdssp** is found on PATH and in conda/mamba environments (e.g.
  `C:\Users\<you>\miniforge3\envs\<env>\Library\bin\mkdssp.exe`); the
  environment does not need to be activated. mkdssp 4.x needs its data files
  (`components.cif`, `mmcif_pdbx.dic`, in `<env>/share/libcifpp`); the tool
  finds them next to the executable, or set `LIBCIFPP_DATA_DIR`. Each mkdssp
  4.x run takes about 2-3 s, mostly to read `components.cif`.
- **VMD** is found on PATH and in `C:\Program Files\VMD` (and the older
  `...\University of Illinois\VMD`). Otherwise give `--vmd_exe`.
- Either can be given explicitly with `--dssp_exe` / `--vmd_exe` (or
  `dssp_exe:` / `vmd_exe:` in config.yaml). The method in use is written to the
  log at the start of each run.

## Quick start

```bash
# UniProt discovery (default query: reviewed EF-hand proteins), all metals
python metal_motif_extractor.py --mode uniprot --metal any --resolution 2.0 --distance 3.0

# Specific metals from your own ID list
python metal_motif_extractor.py --textfile ids.txt --metal CA MG --resolution 2.0 --distance 3.0
python metal_motif_extractor.py --excel proteins.xlsx --metal ZN --distance 2.8

# Coordinated water: on by default; change the cutoff or switch it off
python metal_motif_extractor.py --textfile ids.txt --metal CA --distance 3.0 --water_cutoff 2.8
python metal_motif_extractor.py --textfile ids.txt --metal CA --distance 3.0 --no_water

# Non-protein ligands: on by default; change the cutoff or switch them off
python metal_motif_extractor.py --textfile ids.txt --metal MG --distance 3.0 --ligand_cutoff 2.8
python metal_motif_extractor.py --textfile ids.txt --metal MG --distance 3.0 --no_ligands

# Multi-model (e.g. NMR) files: every model is used by default; use only the first model with
python metal_motif_extractor.py --textfile ids.txt --metal CA --distance 3.0 --include_nmr --first_model_only

# Biological assembly 1 (mmCIF)
python metal_motif_extractor.py --textfile ids.txt --metal CA --distance 3.0 --biological_assembly
```

Input modes are `excel`, `text` and `uniprot`; the mode is inferred from
`--excel` / `--textfile` when `--mode` is not given.

## Several metals, and what is NOT supported

`--metal CA MG` (or `--metal FE FE2`) extracts several metals in one run; each
ion is processed on its own and labelled with its own metal in the IDs, the
files and the `Metal` column, and the summary and the coordination-number check
are given per metal. `FE` and `FE2` are different PDB components (Fe³⁺ and Fe²⁺)
and are kept separate.

**Single-atom ions only.** Iron–sulfur clusters (`FES` 2Fe-2S, `SF4` 4Fe-4S,
`F3S` 3Fe-4S) are single multi-atom residues, not ions, so `--metal FE FE2` does
**not** find the irons of a ferredoxin, and `--metal any` does not either. The
tool says so: structures whose only metal is such a cluster are counted in the
summary ("of these, only a metal cluster") or as `metal_cluster_not_supported`,
and `--metal FES|SF4|F3S` is refused. Heme and other cofactor metals are not
handled either (their metal is an atom inside a larger residue). Use protein
families whose metal is a single ion (for example rubredoxin, EF-hand, zinc
finger).

## Choosing the distance cut-off

**The tool has a default distance cutoff of 2.8 Å; user can change this by `--distance <Å>`** (or write
`distance_cutoff:` in your own config file). Without it the tool stops with a
message. The reason: the right value depends on the metal and the donor atoms,
and a cut-off that is too high counts atoms that are not coordinated, which can
give a coordination number the metal cannot physically have.

What the tool does with your choice:

1. `--distance` applies to protein donor atoms (side-chain O/N/S and backbone
   carbonyl O). Water and non-protein-ligand cut-offs are the same unless you
   set `--water_cutoff` / `--ligand_cutoff`. Per-metal protein cut-offs:
   `--metal_cutoff ZN=2.8 CA=3.2`.
2. After extraction the summary prints the cut-offs used and **compares the
   coordination numbers with the highest value each metal can have**. If sites
   exceed it, it says so, for example:

```
Distance cut-offs used (A)   : protein 4.0, water 4.0, ligand 4.0
Coordination number check (is the cut-off too high for the metal?):
  MG   sites    210 | limit  6 | above limit 31 (14.8%) | highest CN  9 | TOO HIGH
       CN distribution: 4:3 5:20 6:156 7:24 8:5 9:2
  *** WARNING: the distance cut-off (4.0 A) is too high for the coordination
  physically possible for MG: 31 of 210 sites (14.8%) have more than 6 donor
  atoms (highest 9). Lower --distance (or set --metal_cutoff MG=<smaller value>;
  also check --water_cutoff and --ligand_cutoff) and run again.
```


## Coordinated water

For every extracted motif:

1. All water molecules of the same model are collected (from **all** chains;
   deposited waters often carry a different chain ID from the protein).
2. A water is *coordinated* if its **oxygen** is within `water_cutoff`
   (default: the same as `--distance`) of the metal. This is a separate parameter from
   `distance_cutoff` (protein donors) because metal–water bonds are
   typically 2.0–2.6 Å and a 4 Å cutoff would also collect second-shell
   waters.
3. The count, identifiers and distances go to `motif_metadata.csv`
   (`Water Count`, `Water IDs`, `Water Distances (A)`, `Water Cutoff (A)`),
   and every coordinated water (all atoms, original coordinates and chain
   ID) is written into the motif PDB file as `HETATM` records.


## Coordinated non-protein ligands

A non-protein ligand is any residue that is neither an amino acid nor water:
HETATM groups (ATP/ADP, heme, NAD, substrates, inhibitors, sulfate, acetate,
glycerol, ...) and nucleic-acid residues (so RNA/DNA phosphate oxygens
binding Mg²⁺ are found). For every motif:

1. All such residues of the same model are searched, in all chains.
2. Each **N, O or S atom** within `ligand_cutoff` (default: the same as `--distance`) of the
   metal is a coordinated ligand atom (a bidentate phosphate counts two).
   Carbon, phosphorus and hydrogen atoms are never counted.
   Single-atom metal ions (e.g. a neighbouring Na⁺) are never counted as
   ligands.
3. `motif_metadata.csv` gets `Ligand Cutoff (A)`, `Ligand Atom Count`,
   `Ligand Residues` (`ATP:A:501`), `Ligand Atoms` (`ATP:A:501:O2G`),
   `Ligand Distances (A)` and `Common Additive Residues`.
4. The whole ligand residue (all atoms, original coordinates and chain ID)
   is written into the motif PDB file.

Residues with names typical of crystallisation additives (glycerol, ethylene
glycol, PEG, sulfate, acetate, ...) are counted but **flagged** in
`Common Additive Residues`; nothing is excluded. Judge flagged cases
individually (sulfate or acetate can be a genuine ligand).

## Coordination number

The **coordination number (CN) is the number of donor atoms** coordinated to the
ion. There is one donor list per motif, and the CN, and the
geometry come from it:

| Donor atoms in the list | When |
|---|---|
| protein donor atoms (side-chain O/N/S, backbone carbonyl O) of the window's chain | always |
| coordinated water oxygens | unless `--no_water` |
| N/O/S atoms of non-protein ligands | unless `--no_ligands` |

A bidentate carboxylate counts two atoms. With `--no_water --no_ligands` the
CN is just the number of protein donor atoms.

Protein donors: side-chain donors (Asp OD1/OD2, Glu OE1/OE2, Asn OD1,
Gln OE1, Ser OG, Thr OG1, Tyr OH, His ND1/NE2, Cys SG, Met SD) within
`distance_cutoff`, and backbone carbonyl O within the same `distance_cutoff`.
Every qualifying atom is counted.

CSV columns: `Coordination Number` (total donor atoms), `CN (protein donors)`,
`Water Count`, `Ligand Atom Count`, `All Donor Atoms` (every donor atom
counted, so the number can be checked by eye), the geometry columns below, and `Donors in Other Chains`.

## Coordination geometry (template fitting)

The geometry is the **ideal polyhedron that best fits the donor atoms**; there
are no threshold rules.

1. Each donor becomes a unit vector from the metal.
2. For every ideal polyhedron with the same number of vertices (CN 2–9), the
   vectors are fitted to the vertices over all rotations and all
   donor-to-vertex assignments (Kabsch rotation, from
   100 fixed starting rotations, so results are reproducible).
3. `Geometry RMSD (A)` is the RMSD of the best fit: the root-mean-square
   displacement of the donor atoms, placed on the sphere of the mean
   metal–donor distance, from the ideal vertices (angular deviation only).
4. The lowest-RMSD template is `Coordination Geometry`. `Geometry Second
   Best`, `Geometry Second RMSD (A)` and `Geometry Margin (A)` show the
   runner-up and the gap, so marginal assignments can be spotted.

By default **no threshold** is applied: the best fit is always the label, and
the RMSD tells you how good it is. Optionally `--geometry_max_rmsd X` labels
sites whose best fit is worse than X Å "Irregular" (choose X from your
benchmark).

Templates (`ideal_polyhedra.py`) are the ChemEnv reference geometries
(Waroquiers et al., 2017): CN 2 linear, angular; CN 3 trigonal planar,
trigonal pyramidal, T-shaped; CN 4 tetrahedral, square planar, seesaw; CN 5
trigonal bipyramidal, square pyramidal; CN 6 octahedral, trigonal prismatic,
pentagonal pyramidal; CN 7 pentagonal bipyramidal, capped octahedral, two
capped trigonal prisms; CN 8 square antiprismatic, dodecahedral, cubic,
hexagonal bipyramidal, two bicapped trigonal prisms; CN 9 tricapped trigonal
prismatic, capped square antiprismatic. Sites with one donor or more than
nine are not fitted. A geometry that is not in the library is assigned to the
closest member with a correspondingly large RMSD. Edit
`tools/make_ideal_polyhedra.py` and rerun it (needs `pip install pymatgen`) to
change the library.

The template-fit RMSD, with the runner-up and the margin, is the **single
statistic for every coordination number** (CN 2–9). `Coordination Geometry (protein donors
only)` and `Geometry RMSD (protein donors only) (A)` give the fit without water
and ligands.

`Donors in Other Chains` > 0 marks an inter-chain site: those protein donors
are not part of this record's CN (records are per chain), so its CN is
underestimated.

A motif is created only for ions that have at least one protein donor; an
ion coordinated exclusively by water yields no motif.

## Long windows, gaps and inter-chain sites

The window runs from the first to the last protein donor residue of one chain
(padded to the requested length when shorter, equal to the span when longer;
the requested length is a minimum). Donors far apart in sequence therefore give
a long window in which most residues are not in contact with the metal.
Per window the output gives `Binding Span (residues)` (modelled residues),
`Window Exceeds Requested`, `Window Numbering Span` and `Window Missing
Residues` (residue numbers inside the window that are not modelled, i.e. gaps).
An ion with protein donors in several chains gives one record per chain with
the same `Ion Chain`/`Ion Residue Number`/`Ion Serial`; each record's coordination
number and geometry use its own chain's protein donors (plus the ion's water and
ligand atoms), and `Donors in Other Chains` / `Other-Chain Donor Atoms` count and
list the rest. Group records by ion to rebuild the full site. See
`Long_windows_and_interchain_sites.md`.

## Models

By default **every model** of a multi-model file is processed, each on its own
(an ion is only matched to residues and waters of its own model), and the model
number is added to the motif ID and file name. X-ray structures usually have a
single model, so this mainly affects NMR ensembles (`--include_nmr`), where
every model yields its own, nearly identical motif: treat them as redundant.
Use `--first_model_only` (or `--model first`) to analyse only the first model.
DSSP/STRIDE analyse the first model of a file, so windows from other models get
their secondary structure from the dihedral method (the `SS Method` column says
which).

## Motif length report

The run summary reports how many motifs have the requested length, how many
exceeded it (the window is expanded to hold all protein donor residues), how
many are shorter (the chain is shorter than the window), the minimum, median,
mean and maximum length, and a length distribution (one row per length up to the
requested length, one row for each of the next 8 lengths, then bins of doubling
width). The same numbers are saved in `metadata/run_accounting.json`
(`motif_length`) and the distribution in `metadata/motif_length_distribution.csv`.

```
Requested motif length       : 12
  extracted at requested len : 1050 (81.7%)
  exceeded requested length  : 223 (17.4%)   (longest 140)
  shorter than requested     : 12 (0.9%)   (chain shorter than the window)
Motif length (residues)      : min 10, median 12, mean 13.48, max 140
Length distribution:
    length    motifs       %
        12      1050    81.7%  ##############################
        13        60     4.7%  ##
       ...
```

## Run time

Donor atoms are found with a k-d tree built once per model, the motif files
are written from small copies of the needed residues, and DSSP/STRIDE are run
once per structure (not once per motif). On synthetic stress structures
(`tools/stress_test.py`) the original code needed 43 s for 108 motifs and
276 s for 216 motifs; the current code needs 0.6 s and 1.7 s, and 24 s for a
60,000-atom complex with 1,350 ions. `structure_outcomes.csv` has an
`Extraction Seconds` column and structures slower than 60 s are logged, so
slow entries are easy to find.

## Per-metal cutoffs and site provenance

`--metal_cutoff ZN=2.8 CA=3.2` (or `metal_cutoffs:` in `config.yaml`) overrides
`--distance` for those metals; the cutoff used is written to `Distance Cutoff`.
The tool takes the metal identity **as deposited** and does not validate it
(misidentified metals and water-as-metal exist in the PDB; see CheckMyMetal).
To help screen doubtful sites it reports `Ion Occupancy`, `Ion B-factor`,
`Min Donor Distance (A)`, `Max Donor Distance (A)`, `Mean Donor Distance (A)`,
the donor list, geometry RMSD and `Donors in Other Chains`.
`tools/span_report.py` summarises window spans and these columns.

## Comparing two runs (e.g. asymmetric unit vs biological assembly)

Use a different `--output_dir`, `--motif_dir` and `--metadata_dir` for each run
(otherwise the second run overwrites the first one's CSV files and mixes the
motif files), then:

```bash
python tools/compare_runs.py --a asu/metadata/motif_metadata.csv \
       --b assembly/metadata/motif_metadata.csv --out comparison.csv
```

It reports, per PDB ID, motifs, chains and ions in each run and which
structures have fewer, the same, or more motifs. Fewer motifs in the assembly
run usually means the asymmetric unit holds several copies of the assembly.

## Checking PDB vs mmCIF handling

Two different things get mixed up here, so test them separately:

| Comparison | Same chains? | Expected |
|---|---|---|
| asymmetric unit as `.pdb` vs the same asymmetric unit as `.cif` (`--asu_format pdb` / `cif`) | yes | **identical** motifs |
| asymmetric unit vs biological assembly (`--biological_assembly`) | **no** (assembly 1 can be fewer or more chains) | counts **differ** |

```bash
python metal_motif_extractor.py --metal CA --distance 3.0 ... --asu_format pdb --output_dir run_pdb ...
python metal_motif_extractor.py --metal CA --distance 3.0 ... --asu_format cif --output_dir run_cif ...
python tools/compare_runs.py --a run_pdb/metadata/motif_metadata.csv --b run_cif/metadata/motif_metadata.csv
```
Use separate output, motif and metadata folders for each run. Any structure
that differs between the two formats points to a format-handling problem.

## UniProt queries (`--mode uniprot`)

The query string is sent to UniProt **exactly as written**; the tool adds and
changes nothing. So what is searched depends entirely on the syntax you use:

| You pass | What UniProt searches |
|---|---|
| `Zinc finger, C2H2-type` (no field) | free text over **all** fields of **all** entries, reviewed or not: broad, not a domain annotation |
| `(ft_domain:"EF-hand") AND (reviewed:true)` | entries with a **DOMAIN** feature whose description contains the phrase "EF-hand" (also matches "EF-hand 1", "EF-hand 2"), Swiss-Prot only |
| `(ft_zn_fing:"C2H2-type") AND (reviewed:true)` | entries with a **Zinc finger** feature ("C2H2-type 1", ...), Swiss-Prot only |

Text in double quotes is matched as a phrase, not necessarily as the whole
description. Words typed without a field name, such as `--uniprot_query 2Fe-2S ferredoxin-type`, are a free-text search over all entries; use a field query
such as `(ft_domain:"Rubredoxin-type") AND (reviewed:true)` (check the hit count at uniprot.org). UniProt annotates some zinc-finger types as *Zinc finger* features
(`ft_zn_fing`) and not as *DOMAIN* features, so the `ft_domain` zinc-finger
examples in `config.yaml` may return few or no entries; check the hit count at
uniprot.org before a long run. On Windows, quotes inside `--uniprot_query` are
often stripped by the shell; put the query in `config.yaml` (single-quoted YAML
string) instead. The query actually used is written to `logs/log.txt`
(`UniProt discovery mode: query = ...`).

## Processing order and accounting

1. **Acquire IDs.** UniProt mode: query with cursor pagination, collect *all*
   PDB cross-references of every entry, upper-case and de-duplicate them
   (the number of unique IDs is not bounded by the number of entries),
   fetch RCSB metadata once per ID, filter by resolution and method.
2. **Remaining filters.** Metadata is reused (never fetched or filtered
   twice). For Excel/text input the resolution/method filter is applied
   here. The **metal-presence filter** (does the entry contain a requested
   metal?) is applied only when specific metals are requested; with
   `--metal any` it is skipped.
3. **Download** (PDB format; mmCIF for assemblies, or when no PDB file
   exists).
4. **Extraction**, model by model.

The summary reports each removal under its own label (metadata unavailable,
no resolution, method, resolution, requested metal absent), and each
downloaded structure ends in exactly one outcome: `sites_extracted`,
`no_ion_of_requested_type`, `ion_without_protein_ligand`, or `parse_failure`.

## Output files

```
motif_pdbs/                 <PDB>_<chain>_<metal>_motif<n>.pdb   (window residues + ion + coordinated waters + coordinated ligands)
                            (<PDB>_m<model>_… for multi-model files; .cif if a chain ID is longer than 1 character)
metadata/
  structure_metadata.csv    one row per downloaded structure (+ Extraction Status)
  motif_metadata.csv        one row per motif
  secondary_structure.csv   one row per residue per motif
  structure_outcomes.csv    outcome of every downloaded structure
  run_accounting.json       every stage count and the settings of the run
  uniprot_pdb_xrefs.csv     UniProt accession -> PDB IDs (UniProt mode)
logs/log.txt
```

`motif_metadata.csv` columns: Motif ID, PDB ID, Model, Chain, Metal, Ion Chain, Ion Residue Number, Ion Serial,
Resolution, Experimental Method, Binding Residues, Residue Numbers, Donor
Atoms, Coordination Type, Distance Cutoff, Motif Length, Requested Motif Length, Binding Span (residues), Window Exceeds Requested, Sequence, Output
Filename, Coordination Number, Mean Donor Distance (A), Min Donor Distance (A), Max Donor Distance (A), Ion Occupancy, Ion B-factor, Window Numbering Span, Window Missing Residues, Other-Chain Donor Atoms, CN (protein donors), Coordination Geometry,
Geometry RMSD (A), Geometry Second Best, Geometry Second RMSD (A), Geometry
Margin (A), Coordination Geometry (protein donors only), Geometry RMSD
(protein donors only) (A),
Mean/Min/Max LML Angle, Largest Angle 1/2, Water Cutoff (A), Water Count,
Water IDs, Water Distances (A), Ligand Cutoff (A), Ligand Atom Count,
Ligand Residues, Ligand Atoms, Ligand Distances (A), Common Additive
Residues, All Donor Atoms, Donors in Other Chains, SS Method,
SS Sequence, Helix/Strand/Coil Fraction.

## Command line

| Flag | Description |
|---|---|
| `--mode` | `excel`, `text` or `uniprot` |
| `--excel FILE`, `--textfile FILE` | PDB ID input |
| `--uniprot_query …` | UniProt query for `--mode uniprot` |
| `--metal M [M …]` | Metal residue name(s), or `any` |
| `--distance` | **Required.** Protein-donor cutoff (Å); there is no default |
| `--metal_cutoff METAL=VALUE ...` | Per-metal protein-donor cutoff overriding `--distance` |
| `--water_cutoff` | Metal…O(water) cutoff (Å); default: same as `--distance` |
| `--max_cn METAL=N ...` | Highest coordination number of a metal, for the cut-off check (built in: MG 6, ZN 6, CA 9) |
| `--no_water` | Do not search for coordinated water |
| `--ligand_cutoff` | Metal…N/O/S cutoff (Å) for non-protein ligands; default: same as `--distance` |
| `--no_ligands` | Do not search for non-protein ligands |
| `--write_findgeo_input`, `--findgeo_dir` | Write the metal + counted donor atoms of every motif as a PDB for FindGeo (default dir `findgeo_inputs`) |
| `--geometry_max_rmsd` | Optional: label sites whose best geometry fit has RMSD above this (Å) "Irregular" |
| `--model {all,first}` | Which model(s) of multi-model files to use (default `all`) |
| `--first_model_only` | Shortcut for `--model first` |
| `--motif_length` | Window length (residues) |
| `--biological_assembly` | Use assembly 1 (mmCIF) |
| `--asu_format {pdb,cif}` | File format of the asymmetric unit (default `pdb`) |
| `--resolution`, `--include_nmr`, `--methods` | Structure filtering |
| `--workers`, `--proxy`, `--no_ssl_verify` | Performance / network |
| `--config`, `--output_dir`, `--pdb_dir`, `--motif_dir`, `--metadata_dir`, `--log_dir` | Paths |
| `--check_ss`| Test which secondary-structure methods work on this machine  |

CLI flags override `config.yaml`.

## Algorithm notes

For each ion, protein ligands are residues of one chain with a valid donor
atom (side-chain O/N/S, or backbone carbonyl O) within `distance_cutoff`
(the same cutoff for side-chain and backbone donors). The window spans the first to the last ligand
residue and is padded to `motif_length` when the span is shorter (extra
residue to the C-terminal side); when the span is longer, the window is
exactly the span and can exceed `motif_length`. Near chain ends the window
is shifted, not shortened. Chains are processed separately, so a site with
ligands in several chains gives one record per chain; each record counts
only that chain's protein donors (the water count is per ion).

## Limitations

- Only the asymmetric unit is searched unless `--biological_assembly` is
  used; ligands from symmetry mates are otherwise missed.
- Waters and ligands that are not modelled in the deposited structure cannot
  be counted. Ligand donor atoms are limited to N, O and S (halides and
  other donors are not counted).
- The tool operates on deposited coordinates and therefore inherits the source metal assignments and associated issues such as misidentified or poorly modeled metal ions in the PDB.
- Geometry is limited to the template library. 
- Records are per protein chain: protein donors of the same ion in other
  chains are not in the CN (see `Donors in Other Chains`).

## Troubleshooting

| Issue | Likely cause | Resolution |
|---|---|---|
| SSL / `CERTIFICATE_VERIFY_FAILED` | Institutional SSL inspection | `ssl_verify: false` or `--no_ssl_verify` |
| Empty API responses / HTML instead of JSON | Proxy or firewall | `--proxy`; check domains are reachable |
| `No PDB IDs were acquired` | Empty input or UniProt unreachable | Check input file / network |
| Many structures with `no_ion_of_requested_type` | Apo forms or other metals | Expected; see `structure_outcomes.csv` |
| `SS Method` is `BioPDB_dihedral` | mkdssp/VMD not found or not working, or windows from a non-first model | Run `--check_ss` (it prints the reason); see log |
| `--check_ss`: DSSP `FOUND BUT NOT WORKING`, `mmcif_pdbx.dic` / `components.cif` | mkdssp 4 cannot find its data files | Set `LIBCIFPP_DATA_DIR` to `<conda env>/share/libcifpp` |
| `--check_ss`: DSSP `0xC0000135` | A DLL of mkdssp is missing (Windows) | Use the `mkdssp.exe` inside its conda env (`--dssp_exe`), or reinstall `dssp` there |
| VMD/STRIDE gives only coil, or `unable to read output file from Stride` | Space in the temporary-folder path | Handled automatically; if no space-free folder exists, set `TMP` to e.g. `C:\Temp` |

## API references

RCSB Search/Data/File APIs (search.rcsb.org, data.rcsb.org, files.rcsb.org),
UniProt REST API (rest.uniprot.org), Biopython `Bio.PDB`.

## Citation

Berman, H.M. et al. "The Protein Data Bank." *Nucleic Acids Research*, 2000.
The UniProt Consortium. "UniProt: the Universal Protein Knowledgebase."
*Nucleic Acids Research*, 2023.

## License

Specify your project's license here.
