from pathlib import Path
from rdkit import Chem
from rdkit.Chem import RWMol
from Bio.PDB import is_aa

from .constants import DNA_RESIDUES, RNA_RESIDUES, WATER


def parse_seqres_names(pdb_file: Path) -> dict[str, list[str]]:
    """Parse SEQRES records directly, returning chain_id -> ordered list of 3-letter residue names."""
    names: dict[str, list[str]] = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("SEQRES"):
                chain = line[11]          # column 12 (0-indexed 11) = chain ID
                residues = line[19:70].split()
                names.setdefault(chain, []).extend(residues)
    return names


def classify_chain(chain) -> str:
    """Classify a BioPython Chain as protein, DNA, RNA, small molecule, or empty.

    A chain is a potential protein/peptide if it contains at least one canonical
    amino acid. DNA and RNA are identified by exclusive residue name sets.
    Everything else (ligands, ions) is 'small molecule'.
    """
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


def residue_to_rdkit_mol(residue) -> Chem.Mol | None:
    """Build an RDKit mol from a BioPython residue's atoms via a HETATM PDB block. Bonds inferred by proximityBonding"""
    lines = []
    for atom in residue.get_atoms():
        aname   = atom.get_name()
        element = (atom.element or aname.strip()[0]).strip()
        x, y, z = atom.get_vector()
        resname = residue.get_resname().strip()
        resseq  = residue.get_id()[1]
        serial  = atom.get_serial_number()
        lines.append(
            f"HETATM{serial:5d} {aname:<4s} {resname:<3s} A{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    lines.append("END")
    return Chem.MolFromPDBBlock("\n".join(lines), sanitize=True, removeHs=True)


def grab_backbone_atom_indices(mol) -> tuple[int, int, int, int | None]:
    """Identify backbone N, C, and carbonyl O atom indices from an RDKit mol.

    Algorithm:
      1. Find CA by PDB atom name — required anchor for all backbone detection.
      2. Find backbone N: search for nitrogen within 2 bonds of CA, preferring
         a direct (1-bond) neighbor. Raises ValueError if none found — the residue
         is likely not an alpha or beta amino acid.
      3. Find backbone C: carbon directly bonded to CA that has at least one oxygen
         neighbor (carbonyl character). This distinguishes backbone C from CB and
         avoids misidentifying sidechain carbonyls (e.g. Asp CG), which are never
         directly bonded to CA. Raises ValueError if none found.
      4. Find carbonyl O: oxygen neighbor of backbone C with PDB name 'O'.
    """
    # 1. Locate CA
    ca_atom = next(
        (a for a in mol.GetAtoms()
         if (info := a.GetPDBResidueInfo()) and info.GetName().strip() == "CA"),
        None,
    )
    if ca_atom is None:
        raise ValueError("Cannot find CA atom — backbone detection requires an alpha carbon.")
    ca_idx = ca_atom.GetIdx()

    # 2. Find backbone N: prefer direct bond, fall back to 2-bond neighbour
    n_atom = next((nb for nb in ca_atom.GetNeighbors() if nb.GetAtomicNum() == 7), None)
    if n_atom is None:
        for nb in ca_atom.GetNeighbors():
            n_atom = next(
                (nb2 for nb2 in nb.GetNeighbors()
                 if nb2.GetAtomicNum() == 7 and nb2.GetIdx() != ca_atom.GetIdx()),
                None,
            )
            if n_atom is not None:
                break
    if n_atom is None:
        raise ValueError(
            "Cannot find backbone N within 2 bonds of CA — "
            "residue may not be an alpha or beta amino acid."
        )
    n_idx = n_atom.GetIdx()

    # 3. Find backbone C: carbon bonded directly to CA with an oxygen neighbour (carbonyl)
    c_atom = next(
        (nb for nb in ca_atom.GetNeighbors()
         if nb.GetAtomicNum() == 6
         and any(nb2.GetAtomicNum() == 8 for nb2 in nb.GetNeighbors())),
        None,
    )
    if c_atom is None:
        raise ValueError(
            "Cannot find backbone C — no carbon directly bonded to CA has an oxygen neighbour."
        )
    c_idx = c_atom.GetIdx()

    # 4. Find carbonyl O: oxygen neighbour of backbone C named 'O'
    o_idx = next(
        (nb.GetIdx() for nb in c_atom.GetNeighbors()
         if nb.GetAtomicNum() == 8
         and (info := nb.GetPDBResidueInfo()) and info.GetName().strip() == "O"),
        None,
    )

    return n_idx, ca_idx, c_idx, o_idx


def residue_to_smiles(residue) -> str | None:
# def residue_to_smiles(residue, cap: bool = False, zwitterionic: bool = False) -> str | None:
    """Convert a BioPython residue to a canonical SMILES string via RDKit.

    WIP on cap stuff
    cap=True:          ACE on N-terminus, NME on C-terminus.
    zwitterionic=True: NH3+ on N, COO- on C (formal charges, adds OXT).
    """
    mol = residue_to_rdkit_mol(residue)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)
    # if not cap and not zwitterionic:
    #     return Chem.MolToSmiles(mol)

    # n_idx, ca_idx, c_idx, o_idx = grab_backbone_atom_indices(mol)

    # rw = RWMol(mol)

    # if cap:
    #     if n_idx is not None:
    #         ace_c  = rw.AddAtom(Chem.Atom(6))
    #         ace_o  = rw.AddAtom(Chem.Atom(8))
    #         ace_me = rw.AddAtom(Chem.Atom(6))
    #         rw.AddBond(n_idx, ace_c,  Chem.rdchem.BondType.SINGLE)
    #         rw.AddBond(ace_c, ace_o,  Chem.rdchem.BondType.DOUBLE)
    #         rw.AddBond(ace_c, ace_me, Chem.rdchem.BondType.SINGLE)
    #     if c_idx is not None:
    #         nme_n  = rw.AddAtom(Chem.Atom(7))
    #         nme_me = rw.AddAtom(Chem.Atom(6))
    #         rw.AddBond(c_idx,  nme_n,  Chem.rdchem.BondType.SINGLE)
    #         rw.AddBond(nme_n,  nme_me, Chem.rdchem.BondType.SINGLE)

    # if zwitterionic:
    #     if n_idx is not None:
    #         rw.GetAtomWithIdx(n_idx).SetFormalCharge(1)
    #     if c_idx is not None:
    #         oxt = rw.AddAtom(Chem.Atom(8))
    #         rw.GetAtomWithIdx(oxt).SetFormalCharge(-1)
    #         rw.AddBond(c_idx, oxt, Chem.rdchem.BondType.SINGLE)
    #         if o_idx is not None:
    #             bond = rw.GetBondBetweenAtoms(c_idx, o_idx)
    #             if bond:
    #                 bond.SetBondType(Chem.rdchem.BondType.DOUBLE)

    # try:
    #     out_mol = rw.GetMol()
    #     Chem.SanitizeMol(out_mol)
    #     return Chem.MolToSmiles(out_mol)
    # except Exception:
    #     return None
