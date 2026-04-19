from rdkit import Chem

input_file = "RH.smi" # SMILES of hydrocarbons 
output_file = "R_rad.smi" # SMILES of R radicals

def generate_C_radicals(mol):
    radicals = []

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue

        # Find attached hydrogens
        neighbors = atom.GetNeighbors()
        H_neighbors = [n for n in neighbors if n.GetSymbol() == "H"]

        if not H_neighbors:
            continue

        for H in H_neighbors:
            mol_copy = Chem.RWMol(mol)

            c_idx = atom.GetIdx()
            h_idx = H.GetIdx()

            # Remove that hydrogen atom
            mol_copy.RemoveAtom(h_idx)

            # After removal, atom indices may shift, but C index is still valid
            c_atom = mol_copy.GetAtomWithIdx(c_idx)
            c_atom.SetNumRadicalElectrons(1)

            Chem.SanitizeMol(mol_copy)

            radicals.append(Chem.RemoveHs(mol_copy))

    return radicals


n_written = 0

with open(input_file) as infile, open(output_file, "w") as outfile:

    for line in infile:
        line = line.strip()
        if not line:
            continue

        # skip header or any non-data line
        if line.startswith("SMILES"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        smiles, mol_id = parts[0], parts[1]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        radicals = generate_C_radicals(mol)

        seen = set()
        rad_idx = 1

        for rad in radicals:
            rad_smi = Chem.MolToSmiles(rad, isomericSmiles=True)
            if rad_smi in seen:
                continue
            seen.add(rad_smi)

            new_id = f"{mol_id}_C_rad_{rad_idx:05d}"
            outfile.write(f"{rad_smi}\t{new_id}\n")
            rad_idx += 1
            n_written += 1

print("Total C-centered radicals written:", n_written)
print("Output:", output_file)
