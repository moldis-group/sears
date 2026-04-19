import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToXYZBlock

input_file = "file.smi"

output_root = "geoms"      # in the folder 'geoms' multiple folders will be created with same name as molecule_id in the second column of "file.smi". Each folder  will contain rdkit_guess.xyz

os.makedirs(output_root, exist_ok=True)

with open(input_file, "r") as f:
    first_line = True
    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header like: SMILES  ID
        if first_line and line.upper().startswith("SMILES"):
            first_line = False
            continue
        first_line = False

        parts = line.split()
        if len(parts) < 2:
            continue

        smiles = parts[0]
        name = parts[1]  # folder and file prefix

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Failed to parse SMILES: {smiles}")
            continue

        # explicit hydrogens
        mol = Chem.AddHs(mol)

        # ----- 3D embedding -----
        params = AllChem.ETKDGv3()
        params.randomSeed = 42

        conf_id = AllChem.EmbedMolecule(mol, params)

        if conf_id == -1:
            print(f"ETKDG failed for: {smiles}, trying random coords.")
            params.useRandomCoords = True
            conf_id = AllChem.EmbedMolecule(mol, params)

        if conf_id == -1:
            print(f"Embedding really failed for: {smiles}")
            continue

        # geometry optimization (UFF)
        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)

        # ----- write XYZ into its own folder -----
        xyz_block = MolToXYZBlock(mol)

        mol_folder = os.path.join(output_root, name)
        os.makedirs(mol_folder, exist_ok=True)

        out_path = os.path.join(mol_folder, "rdkit_guess.xyz")
        with open(out_path, "w") as out:
            out.write(xyz_block)

        print(f"Wrote {out_path}")

print("Done.")
