# python3 deduplicate.py all.smi > unique.smi

#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdchem
import sys

def invert_chiral_tag(tag):
    if tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
    elif tag == rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
        return rdchem.ChiralType.CHI_TETRAHEDRAL_CW
    return tag

def mirror_mol(mol):
    """
    Return a copy with all tetrahedral stereocenters inverted.
    Double-bond E/Z stereochemistry is left unchanged.
    """
    m = Chem.Mol(mol)
    for atom in m.GetAtoms():
        tag = atom.GetChiralTag()
        atom.SetChiralTag(invert_chiral_tag(tag))
    return m

def enantiomer_invariant_key(smiles):
    """
    Same key for an enantiomeric pair, different keys for diastereomers.
    E/Z stereochemistry is preserved.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Ensure stereochemistry is fully assigned from the SMILES
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    smi1 = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    mir = mirror_mol(mol)
    Chem.AssignStereochemistry(mir, force=True, cleanIt=True)
    smi2 = Chem.MolToSmiles(mir, canonical=True, isomericSmiles=True)

    return min(smi1, smi2)

def deduplicate_enantiomers(infile):
    seen = set()
    with open(infile) as f:
        for line in f:
            line = line.strip()

            parts = line.split()
            if len(parts) < 2:
                continue

            smi = parts[0]
            name = parts[1]

            key = enantiomer_invariant_key(smi)
            if key is None:
                print(f"Skipping invalid SMILES: {smi}", file=sys.stderr)
                continue

            if key not in seen:
                seen.add(key)
                print(f"{smi}\t{name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} mol.smi")
        sys.exit(1)

    deduplicate_enantiomers(sys.argv[1])
