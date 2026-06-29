 #This script reports useful information about a loaded PDB.
# It also builds a per-residue SMILES dictionary for downstream RDKit featurization.

import csv
import tomllib
from datetime import date
from pathlib import Path
from typing import Annotated

import typer
from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa

from toolkit_pdb.constants import NON_AA
from toolkit_pdb.pdb_utils import classify_chain, residue_to_smiles, parse_seqres_names

app = typer.Typer()


@app.command()
def main(config_path: Annotated[Path, typer.Argument(help="Path to config.toml")]):
    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    repo_root  = Path(__file__).parent.parent
    pdb_file   = repo_root / cfg["project"]["pdb_file"]
    datestamp  = date.today().strftime("%Y_%m_%d")
    output_dir = repo_root / cfg["project"]["output_dir"] / f"{cfg['project']['name']}_{datestamp}"

    # Cap stuff is WIP
    # cap_termini  = cfg["smiles"]["cap_termini"]
    # zwitterionic = cfg["smiles"]["zwitterionic"]

    output_dir.mkdir(parents=True, exist_ok=True)

    # I want the parser to throw warnings if multiple residues have the same identifier
    # as per https://biopython.org/docs/dev/Tutorial/chapter_pdb.html#sec-problem-structures
    parser     = PDBParser(PERMISSIVE=0)
    loaded_pdb = parser.get_structure("loaded_pdb", pdb_file)

    n_models = len(loaded_pdb)
    print(f"PDB file: {pdb_file} has {n_models} model(s)")
    print("Selecting first one")  # add utility for conformer selection later if we think other models will have variations in sequence
    selected_structure = loaded_pdb[0]

    # SEQRES length per chain — reflects the full deposited sequence
    seqres_lengths: dict[str, int] = {}
    for record in SeqIO.parse(str(pdb_file), "pdb-seqres"):
        chain_id = record.id.split(":")[-1]
        seqres_lengths[chain_id] = len(record.seq)

    # SEQRES 3-letter residue names in order — used to identify missing residues by name
    seqres_names = parse_seqres_names(pdb_file)

    # resname -> canonical SMILES from first occurrence
    aa_smiles: dict[str, str] = {}
    # resname -> all distinct SMILES seen (populated only on conflict)
    aa_smiles_conflicts: dict[str, set[str]] = {}

    for chain in selected_structure:
        chain_id   = chain.get_id()
        chain_type = classify_chain(chain)  # Update later to not assume that a protein/peptide chain has at least one canonical AA
        print(f"Chain {chain_id}: {chain_type}")

        if chain_type == "potential protein or peptide":
            valid_AAs = [r for r in chain.get_residues() if r.get_resname().strip() not in NON_AA]
            n_coords  = len(valid_AAs)
            n_seqres  = seqres_lengths.get(chain_id, "?")
            missing   = (n_seqres - n_coords) if isinstance(n_seqres, int) else "?"
            print(f"  SEQRES residues (full sequence):  {n_seqres}")
            print(f"  Residues with coordinates:        {n_coords}")
            print(f"  Missing coordinates:              {missing}")

            canonical     = [r for r in valid_AAs if is_aa(r, standard=True)]
            non_canonical = [r for r in valid_AAs if not is_aa(r, standard=True)]
            print(f"  Canonical AAs:                    {len(canonical)}")

            if non_canonical:
                nc_info = []
                for name in sorted({r.get_resname().strip() for r in non_canonical}):
                    rep = next(r for r in non_canonical if r.get_resname().strip() == name)
                    # N-methylations tagging WIP
                    # tag = " (N-methylated)" if is_n_methylated(rep) else ""
                    tag = ""
                    nc_info.append(f"{name}{tag}")
                print(f"  Non-canonical residues:           {', '.join(nc_info)}")

            for residue in valid_AAs:
                resname = residue.get_resname().strip()
                smiles  = residue_to_smiles(residue) #, cap=cap_termini, zwitterionic=zwitterionic) WIP
                if smiles is None:
                    print(f"  WARNING: could not generate SMILES for {resname} {residue.get_id()[1]}")
                    continue
                if resname not in aa_smiles:
                    aa_smiles[resname] = smiles
                elif smiles != aa_smiles[resname]:
                    aa_smiles_conflicts.setdefault(resname, {aa_smiles[resname]}).add(smiles)

    print("\nAA SMILES dictionary:")
    for resname, smiles in sorted(aa_smiles.items()):
        print(f"  {resname}: {smiles}")

    if aa_smiles_conflicts:
        print("\nFLAG — conflicting SMILES for the same residue name:")
        for resname, smiles_set in sorted(aa_smiles_conflicts.items()):
            print(f"  {resname}:")
            for s in smiles_set:
                print(f"    {s}")

    smiles_csv = output_dir / "smiles.csv"
    with open(smiles_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["resname", "smiles"])
        for resname, smiles in sorted(aa_smiles.items()):
            writer.writerow([resname, smiles])
    print(f"\nSaved: {smiles_csv}")


if __name__ == "__main__":
    app()
