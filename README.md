# toolkit-for-pdb

Scripts for analyzing structures from the PDB.

## Current Functionality
Takes an input pdb, parses for amino acids (canonical and non canonical), generates SMILES for each amino acid, computes RDKit molecular descriptors, and projects a chosen feature onto B-factors for visualization. 

## Up Next
In progress to-dos (for better/more accurate functionality):
1. Properly cap the scraped amino acids 
2. Handle multi-chain PDBs with missing residues, and have convenient reporting
3. Handle "duplicate" amino acids that allegedly map to unique SMILES

In progress to-dos (nice-to-haves):
1. Try out other featurizations (other cheminformatics tools, chemical language model embeddings, etc)
2. Flag interesting ncAAs (N-methylation, beta amino acids,etc)

So far have only tried on straight-forward cyclosporine example (PDBID: 1CYB). 

## Overview

The pipeline has three stages, all configured by a single `config.toml` and run together via `run_pipeline.sh`:

```
scrape.py  →  smiles.csv  →  featurize.py  →  features.csv  →  replace_bfactors.py  →  *_bfactor_<feature>.pdb
```

| Script | Input | Output |
|---|---|---|
| `scrape.py` | PDB file | `smiles.csv` — one SMILES per unique residue type |
| `featurize.py` | `smiles.csv` | `features.csv` — 217 RDKit descriptors per residue |
| `replace_bfactors.py` | `features.csv` + PDB | PDB with B-factors replaced by a chosen descriptor |

Each run creates a sub-folder under `outputs/<project_name>/` so results from different configs don't overwrite each other.

---

## Setup

The project uses [pixi](https://pixi.sh) to manage the environment (conda + PyPI dependencies).

```bash
# Install pixi if you don't have it
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repo and then install all dependencies with

pixi install
```

Dependencies include: Python 3.14, RDKit, BioPython, pandas, typer.

---

## Configuration

Copy and edit `config.toml` to set up a new run:

```toml
[project]
name       = "cyclosporin_a"   # sub-folder name under outputs/
pdb_file   = "data/1CYB.pdb"  # path to the input PDB (relative to repo root)
output_dir = "outputs"

[smiles]
# WIP for capping
cap_termini  = false   # if true: add ACE to N-terminus, NME to C-terminus
zwitterionic = false   # if true: apply NH3+ / COO- formal charges (mutually exclusive with cap_termini)

[bfactor]
feature_column = "MolLogP"   # which column from features.csv to write into B-factors
```

---

## Running the pipeline

```bash
# Run all three scripts in sequence
./run_pipeline.sh config.toml

# Or run each script individually
pixi run python scripts/scrape.py          config.toml
pixi run python scripts/featurize.py       config.toml
pixi run python scripts/replace_bfactors.py config.toml
```

Each script also supports `--help`.

---

## Output structure

```
outputs/
└── cyclosporin_a/
    ├── smiles.csv                    # resname, smiles
    ├── features.csv                  # resname, smiles, MolWt, MolLogP, ... (217 descriptors)
    └── 1CYB_bfactor_MolLogP.pdb     # original PDB with B-factors = MolLogP per residue type
```

To visualize in PyMOL/ChimeraX/VMD, color by B-factor after loading the output PDB.

---

## What scrape.py reports

For each chain it prints:

- Chain type: `potential protein or peptide`, `DNA`, `RNA`, or `small molecule`. NOTE: current assumption is that a "potential protein or peptide" chain contains at least one canonical amino acid.
- SEQRES residue count vs. residues with coordinates (and which specific residues are missing)
- Count of canonical vs. non-canonical amino acids
- Any conflicting SMILES — i.e., two instances of the same residue name that produce different SMILES (flags stereochemistry issues or atom-level differences)

SMILES are generated from the atoms present in the structure (not a lookup table), so they reflect what is actually resolved.

---

## Repo structure

```
toolkit-for-pdb/
├── toolkit_pdb/           # installable package with shared helpers
│   ├── constants.py       # residue-type sets, bond distance thresholds
│   └── pdb_utils.py       # classify_chain, is_n_methylated, residue_to_smiles,
│                          #   residue_to_rdkit_mol, grab_backbone_atom_indices,
│                          #   parse_seqres_names
├── scripts/
│   ├── scrape.py          # stage 1: PDB → smiles.csv
│   ├── featurize.py       # stage 2: smiles.csv → features.csv
│   └── replace_bfactors.py# stage 3: features.csv + PDB → bfactor PDB
├── data/                  # input PDB files
├── outputs/               # one sub-folder per config run (git-ignored)
├── config.toml            # example configuration
├── run_pipeline.sh        # runs all three stages in sequence
├── pyproject.toml         # makes toolkit_pdb pip-installable (editable)
└── pixi.toml              # environment definition
```
