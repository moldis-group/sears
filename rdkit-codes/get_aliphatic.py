from rdkit import Chem
  
input_file = "full.smi"            # contails both aliphatic and aromatic hydrocarbons
output_aliphatic = "aliphatic.smi" # aliphatic subset
output_aromatic = "atomatic.smi"   # aromatic subset

n_total = 0
n_aliphatic = 0
n_aromatic = 0

with open(input_file) as infile, \
     open(output_aliphatic, "w") as f_aliph, \
     open(output_aromatic, "w") as f_arom:

    for line in infile:
        line = line.strip()

        parts = line.split()
        if len(parts) < 2:
            continue

        smiles, mol_id = parts[0], parts[1]
        n_total += 1

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        # any aromatic ring?
        is_aromatic = any(bond.GetIsAromatic() for bond in mol.GetBonds())

        if is_aromatic:
            f_arom.write(f"{smiles}\t{mol_id}\n")
            n_aromatic += 1
        else:
            f_aliph.write(f"{smiles}\t{mol_id}\n")
            n_aliphatic += 1

print("Total hydrocarbons:", n_total)
print("Aliphatic (used in main pipeline):", n_aliphatic)
print("Aromatic (benzene/toluene etc):", n_aromatic)
print("Aliphatic file:", output_aliphatic)
print("Aromatics file:", output_aromatic)
