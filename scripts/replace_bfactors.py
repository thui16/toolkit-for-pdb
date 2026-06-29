# Reads features.csv produced by featurize.py and writes a new PDB file with
# B-factors replaced by the selected feature column value, per residue type.

import tomllib
from datetime import date
from pathlib import Path
from typing import Annotated

import typer
import pandas as pd
from Bio.PDB import PDBParser, PDBIO

app = typer.Typer()


@app.command()
def main(config_path: Annotated[Path, typer.Argument(help="Path to config.toml")]):
    with open(config_path, "rb") as f:
        cfg = tomllib.load(f)

    repo_root      = Path(__file__).parent.parent
    pdb_file       = repo_root / cfg["project"]["pdb_file"]
    datestamp      = date.today().strftime("%Y_%m_%d")
    output_dir     = repo_root / cfg["project"]["output_dir"] / f"{cfg['project']['name']}_{datestamp}"
    feature_column = cfg["bfactor"]["feature_column"]

    features_csv = output_dir / "features.csv"
    df = pd.read_csv(features_csv, index_col="resname")

    if feature_column not in df.columns:
        raise typer.BadParameter(
            f"Feature '{feature_column}' not found in {features_csv}.\n"
            f"Available: {', '.join(df.columns)}"
        )

    feature_map: dict[str, float] = df[feature_column].to_dict()

    parser     = PDBParser(PERMISSIVE=0)
    structure = parser.get_structure("loaded_pdb", pdb_file)

    missing_resnames: set[str] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                val = feature_map.get(resname)
                if val is None:
                    missing_resnames.add(resname)
                    val = 0.0
                for atom in residue:
                    atom.set_bfactor(float(val))

    if missing_resnames:
        print(f"WARNING: no feature value for {sorted(missing_resnames)} — B-factor set to 0.0")

    out_pdb = output_dir / f"{pdb_file.stem}_bfactor_{feature_column}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_pdb))
    print(f"Replaced B-factors with '{feature_column}'.")
    print(f"Saved: {out_pdb}")


if __name__ == "__main__":
    app()
