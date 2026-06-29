# Reads smiles.csv produced by scrape.py, computes RDKit molecular descriptors
# for each residue, and saves features.csv with one column per descriptor.

import tomllib
from datetime import date
from pathlib import Path
from typing import Annotated

import typer
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

app = typer.Typer()


@app.command()
def main(config_path: Annotated[Path, typer.Argument(help="Path to config.toml")]):
    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    repo_root  = Path(__file__).parent.parent
    datestamp  = date.today().strftime("%Y_%m_%d")
    output_dir = repo_root / cfg["project"]["output_dir"] / f"{cfg['project']['name']}_{datestamp}"

    smiles_csv   = output_dir / "smiles.csv"
    features_csv = output_dir / "features.csv"

    df = pd.read_csv(smiles_csv)

    rows = []
    for _, row in df.iterrows():
        temp_mol = Chem.MolFromSmiles(row["smiles"])
        desc = {name: fn(temp_mol) for name, fn in Descriptors.descList}
        if not desc:
            print(f"WARNING: could not parse SMILES for {row['resname']} — skipping descriptors.")
        rows.append({"resname": row["resname"], "smiles": row["smiles"], **desc})

    out = pd.DataFrame(rows)
    out.to_csv(features_csv, index=False)
    print(f"Computed {len(out.columns) - 2} descriptors for {len(out)} residues.")
    print(f"Saved: {features_csv}")


if __name__ == "__main__":
    app()
