from rdkit import Chem
  
input_file = "full.smi" # contains systems with any number of rings or acyclic
output_single = "single_or_acyclic.smi" # subset with one or no ring
output_multiring = "multiring.smi" # subset with multiple rings

n_total = 0
n_single_or_acyclic = 0
n_multiring = 0

with open(input_file) as infile, \
     open(output_single, "w") as f_single, \
     open(output_multiring, "w") as f_multi:


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

        n_rings1 = mol.GetRingInfo().NumRings()
        ring_info = mol.GetRingInfo()
        n_rings = ring_info.NumRings()

        print(n_rings1,n_rings)

        if n_rings <= 1:
            f_single.write(f"{smiles}\t{mol_id}\n")
            n_single_or_acyclic += 1
        else:
            f_multi.write(f"{smiles}\t{mol_id}\n")
            n_multiring += 1

print("Total hydrocarbons:", n_total)
print("Acyclic or single-ring:", n_single_or_acyclic)
print("Multi-ring:", n_multiring)
print("Kept in:", output_single)
print("Excluded (multi-ring) in:", output_multiring)
