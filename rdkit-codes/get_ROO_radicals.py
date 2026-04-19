#!/usr/bin/env python3

from rdkit import Chem

input_file = "R_rad.smi" # SMILES of R radicals
output_file = "ROO.smi"  # SMILES of ROO

from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)

opts = StereoEnumerationOptions(
    tryEmbedding=False,
    unique=True,
    maxIsomers=16  # enough for 1–2 centers
)


input_file = "../step_4_R_radical/bigQM7w_hydrocarbons_C_radicals_deduplicated.smi"
output_roo = "bigQM7w_ROO_radicals_stereo_maxIsomers_16.smi"


def get_radical_carbons(mol):
    """Return indices of C atoms with one radical electron."""
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "C" and a.GetNumRadicalElectrons() == 1]

def generate_ROO_from_R(mol):
    mol_H = Chem.AddHs(mol)
    results = []

    for c_idx in get_radical_carbons(mol_H):
        rw = Chem.RWMol(mol_H)
        c_atom = rw.GetAtomWithIdx(c_idx)

        # add O–O
        o1_idx = rw.AddAtom(Chem.Atom("O"))
        o2_idx = rw.AddAtom(Chem.Atom("O"))
        rw.AddBond(c_idx, o1_idx, Chem.BondType.SINGLE)
        rw.AddBond(o1_idx, o2_idx, Chem.BondType.SINGLE)

        # move radical C· -> O2·
        c_atom.SetNumRadicalElectrons(0)
        o2_atom = rw.GetAtomWithIdx(o2_idx)
        o2_atom.SetNumRadicalElectrons(1)

        try:
            base = rw.GetMol()
            Chem.SanitizeMol(base)
            base = Chem.RemoveHs(base)

            # unset chirality at attachment carbon so RDKit will enumerate both
            attach_atom = base.GetAtomWithIdx(c_idx)
            attach_atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)

            # enumerate stereoisomers
            isos = list(EnumerateStereoisomers(base, options=opts))

            seen_local = set()
            for iso in isos:
                smi = Chem.MolToSmiles(iso, isomericSmiles=True)
                if smi in seen_local:
                    continue
                seen_local.add(smi)
                results.append(iso)

        except Exception:
            continue

    return results


global_seen_roo = set()

n_in = 0
n_roo = 0

with open(input_file) as infile, \
     open(output_roo, "w") as f_roo:

    for line in infile:
        line = line.strip()

        parts = line.split()
        if len(parts) < 2:
            continue

        smiles, mol_id = parts[0], parts[1]
        n_in += 1

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        roo_mols = generate_ROO_from_R(mol)

        roo_idx = 1
        for rm in roo_mols:
            roo_smi = Chem.MolToSmiles(rm, isomericSmiles=True)
            if roo_smi in global_seen_roo:
                continue
            global_seen_roo.add(roo_smi)

            new_id = f"{mol_id}_ROO_{roo_idx:05d}"
            f_roo.write(f"{roo_smi}\t{new_id}\n")
            roo_idx += 1
            n_roo += 1

print("Input C-radicals read:", n_in)
print("Unique ROO· radicals written:", n_roo)
print("ROO· file:", output_roo)

