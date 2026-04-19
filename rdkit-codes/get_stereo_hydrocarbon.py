from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions
)

# Example: 2-butene
smiles = "CC=CC"
mol_id = "bigQM7w_000098"

# Stereo enumeration options
opts = StereoEnumerationOptions(
    tryEmbedding=False,
    unique=True,
    maxIsomers=1024
)

def is_hydrocarbon_mol(mol: Chem.Mol) -> bool:
    """True if molecule contains only C and H atoms."""
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ("C", "H"):
            return False
    return True

mol = Chem.MolFromSmiles(smiles)

if mol is None:
    raise ValueError(f"Invalid SMILES: {smiles}")

if not is_hydrocarbon_mol(mol):
    raise ValueError(f"Not a hydrocarbon: {smiles}")

# Write no-stereo version
mol_no_stereo = Chem.Mol(mol)
Chem.RemoveStereochemistry(mol_no_stereo)

canonical_no_stereo = Chem.MolToSmiles(
    mol_no_stereo,
    isomericSmiles=False
)

print("No-stereo SMILES:")
print(f"{canonical_no_stereo}\t{mol_id}")

# Enumerate stereoisomers
isomers = list(EnumerateStereoisomers(mol, options=opts))

seen = set()
stereo_idx = 1

print("\nStereoisomers:")
for iso in isomers:
    iso_smi = Chem.MolToSmiles(iso, isomericSmiles=True)
    if iso_smi in seen:
        continue
    seen.add(iso_smi)

    stereo_name = f"{mol_id}_stereo_{stereo_idx:05d}"
    print(f"{iso_smi}\t{stereo_name}")
    stereo_idx += 1

print(f"\nTotal stereoisomers found: {len(seen)}")
