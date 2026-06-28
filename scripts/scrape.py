#This script is used to report useful information to the user about the loaded PDB.
#Information includes
#It also builds AA level dictionary with smiles values.
#This dictionary will be used downstream for RDKit calculations, B-factor replacements, visualizations.

# Assumptions: no AAs are named "HOH" nor are only 2 letters that conventionally are for nucleic acids.
# Also selects one first conformer

from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa

PDB_FILE = Path(__file__).parent.parent / "data" / "1CYB.pdb"
# PDB_FILE = Path(__file__).parent.parent / "data" / "3VD0.pdb"

DNA_RESIDUES = {"DA", "DT", "DC", "DG", "DI"}
RNA_RESIDUES = {"A", "U", "C", "G", "I"}
WATER = "HOH"

# I want the parser to throw warnings if multiple residues have the same identifier
# as per https://biopython.org/docs/dev/Tutorial/chapter_pdb.html#sec-problem-structures
parser = PDBParser(PERMISSIVE=0)
loaded_pdb = parser.get_structure("loaded_pdb", PDB_FILE)

n_models = len(loaded_pdb)
print(f"PDB file: {PDB_FILE} has {n_models} model(s)")
print(f"Selecting first one") # add utility for conformer selection later if we think other models will have variations in sequence

selected_structure = loaded_pdb[0]

protein_residues = []
dna_residues = []


# SEQRES length per chain — reflects the full deposited sequence
seqres_lengths = {}
for record in SeqIO.parse(str(PDB_FILE), "pdb-seqres"):
    chain_id = record.id.split(":")[-1]
    seqres_lengths[chain_id] = len(record.seq)

def classify_chain(chain):
    # Assumes that to be a potential peptide or protein, you have at least one canonical amino acid. 
    # I should update to classify as presence of at least one peptide (amide) bond
    residues = [r for r in chain.get_residues() if r.get_resname().strip() != WATER]
    if not residues:
        return "empty"
    if any(is_aa(r, standard=True) for r in residues):
        return "potential protein or peptide"
    names = {r.get_resname().strip() for r in residues}
    if names <= DNA_RESIDUES:
        return "DNA"
    if names <= RNA_RESIDUES:
        return "RNA"
    return "small molecule"

for chain in selected_structure:
    chain_id = chain.get_id()
    chain_type = classify_chain(chain)
    print(f"Chain {chain_id}: {chain_type}")

    if chain_type == "potential protein or peptide":
        n_coords  = sum(
            1 for r in chain.get_residues()
            if r.get_resname().strip() not in DNA_RESIDUES | RNA_RESIDUES | {WATER}
        )
        n_seqres  = seqres_lengths.get(chain_id, "?")
        missing   = (n_seqres - n_coords) if isinstance(n_seqres, int) else "?"
        print(f"  SEQRES residues (full sequence):  {n_seqres}")
        print(f"  Residues with coordinates:        {n_coords}")
        print(f"  Missing coordinates:              {missing}")

        non_nucleic = [
            r for r in chain.get_residues()
            if r.get_resname().strip() not in DNA_RESIDUES | RNA_RESIDUES | {WATER}
        ]
        canonical     = [r for r in non_nucleic if is_aa(r, standard=True)]
        non_canonical = [r for r in non_nucleic if not is_aa(r, standard=True)]
        print(f"  Canonical AAs:                    {len(canonical)}")
        if non_canonical:
            names = sorted({r.get_resname().strip() for r in non_canonical})
            print(f"  Non-canonical residues:           {', '.join(names)}")

# print(selected_structure.header["has_missing_residues"])