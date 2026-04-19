#!/usr/bin/env python3

from rdkit import Chem
from collections import deque

input_file = "ROO.smi" # SMILES of ROO radicals
hyb='sp2'
output_file = "QOOH_"+hyb+".smi"  # SMILES of QOOH

MAX_SHELL = 8  # 0=alpha, 1=beta, 2=gamma, 3=delta, 4=epsilon, ...


def find_peroxy_pattern(mol):
    """
    Given a ROO· radical, identify:
    - idx_O_rad: O atom with 1 radical electron (O·)
    - idx_O_partner: the other O in the O–O bond
    - idx_C_anchor: carbon attached to the peroxy group (C–O–O·)

    Returns:
        (idx_O_rad, idx_O_partner, idx_C_anchor)
    or None if pattern not found.
    """
    o_rad_idx = None
    for a in mol.GetAtoms():
        if a.GetSymbol() == "O" and a.GetNumRadicalElectrons() == 1:
            o_rad_idx = a.GetIdx()
            break
    if o_rad_idx is None:
        return None

    o_rad_atom = mol.GetAtomWithIdx(o_rad_idx)

    # find O–O partner
    o_partner_idx = None
    for nei in o_rad_atom.GetNeighbors():
        if nei.GetSymbol() == "O":
            o_partner_idx = nei.GetIdx()
            break
    if o_partner_idx is None:
        return None

    o_partner_atom = mol.GetAtomWithIdx(o_partner_idx)

    # find attached carbon (anchor)
    c_anchor_idx = None
    for nei in o_partner_atom.GetNeighbors():
        if nei.GetSymbol() == "C":
            c_anchor_idx = nei.GetIdx()
            break
    if c_anchor_idx is None:
        return None

    return o_rad_idx, o_partner_idx, c_anchor_idx


def carbon_distances(mol, anchor_idx):
    """
    BFS on carbon-only graph to get distances (in C–C bonds) from anchor carbon.

    Returns:
        dict {C_idx: distance}
    """
    dist = {anchor_idx: 0}
    q = deque([anchor_idx])

    while q:
        i = q.popleft()
        atom_i = mol.GetAtomWithIdx(i)
        for nei in atom_i.GetNeighbors():
            j = nei.GetIdx()
            if mol.GetAtomWithIdx(j).GetSymbol() != "C":
                continue
            if j not in dist:
                dist[j] = dist[i] + 1
                q.append(j)
    return dist


def is_sp_carbon(atom):
    """
    True only for carbon atoms assigned sp hybridization by RDKit.
    """
    return (
        atom.GetSymbol() == "C"
        and atom.GetHybridization() == Chem.HybridizationType.SP
    )

def is_sp2_carbon(atom):
    """
    True only for carbon atoms assigned sp2 hybridization by RDKit.
    """
    return (
        atom.GetSymbol() == "C"
        and atom.GetHybridization() == Chem.HybridizationType.SP2
    )

def is_sp3_carbon(atom):
    """
    True only for carbon atoms assigned sp3 hybridization by RDKit.
    """
    return (
        atom.GetSymbol() == "C"
        and atom.GetHybridization() == Chem.HybridizationType.SP3
    )


def generate_QOOH_from_ROO(mol, max_shell=4):
    """
    Given a ROO· radical mol, generate all QOOH radicals via intramolecular
    H abstraction from carbons within distance 0..max_shell from the anchor
    carbon along the carbon skeleton.

    Only donor carbons that are sp3 are allowed.

    Returns:
        list of RDKit Mol objects (QOOH radicals)
    """
    # Work with explicit Hs
    mol_H = Chem.AddHs(mol)

    pattern = find_peroxy_pattern(mol_H)
    if pattern is None:
        return []

    o_rad_idx, o_partner_idx, c_anchor_idx = pattern

    # Distances on carbon-only graph
    dist = carbon_distances(mol_H, c_anchor_idx)

    qooh_mols = []

    for c_idx, d in dist.items():
        if d > max_shell:
            continue

        c_atom = mol_H.GetAtomWithIdx(c_idx)

        # donor carbon must be sp/sp2/sp3
        if hyb == 'sp': 
            if not is_sp_carbon(c_atom):
                continue
        if hyb == 'sp2': 
            if not is_sp2_carbon(c_atom):
                continue
        if hyb == 'sp3': 
            if not is_sp3_carbon(c_atom):
                continue

        # donor carbon must have at least one H
        H_neighbors = [n for n in c_atom.GetNeighbors() if n.GetSymbol() == "H"]
        if not H_neighbors:
            continue

        # use only one H per donor carbon position
        donor_H = H_neighbors[0]

        rw = Chem.RWMol(mol_H)

        # remove donor H (H abstraction)
        donor_h_idx = donor_H.GetIdx()
        rw.RemoveAtom(donor_h_idx)

        # heavy atom indices remain unchanged since only H is removed

        # add H to the radical O, converting ROO· -> ROOH
        o_rad_atom = rw.GetAtomWithIdx(o_rad_idx)
        new_H_idx = rw.AddAtom(Chem.Atom("H"))
        rw.AddBond(o_rad_idx, new_H_idx, Chem.BondType.SINGLE)
        o_rad_atom.SetNumRadicalElectrons(0)

        # donor carbon becomes new radical center
        c_atom_new = rw.GetAtomWithIdx(c_idx)
        c_atom_new.SetNumRadicalElectrons(1)

        try:
            m_new = rw.GetMol()
            Chem.SanitizeMol(m_new)
            m_new = Chem.RemoveHs(m_new)
            qooh_mols.append(m_new)
        except Exception:
            # skip invalid constructs
            continue

    return qooh_mols


def main():
    global_seen_qooh = set()
    n_in = 0
    n_qooh = 0

    with open(input_file) as infile, open(output_file, "w") as out:
        out.write("SMILES\tID\n")

        for line in infile:
            line = line.strip()
            if not line or line.startswith("SMILES"):
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            smiles, mol_id = parts[0], parts[1]
            n_in += 1

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            qooh_mols = generate_QOOH_from_ROO(mol, max_shell=MAX_SHELL)

            idx = 1
            for qm in qooh_mols:
                smi = Chem.MolToSmiles(qm, isomericSmiles=True)
                if smi in global_seen_qooh:
                    continue

                global_seen_qooh.add(smi)
                new_id = f"{mol_id}_QOOH_{idx:05d}"
                out.write(f"{smi}\t{new_id}\n")
                idx += 1
                n_qooh += 1

    print("ROO· radicals read:", n_in)
    print("Unique QOOH radicals written:", n_qooh)
    print("Output QOOH file:", output_file)


if __name__ == "__main__":
    main()
