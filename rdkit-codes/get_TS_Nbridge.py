#
# Due to the following line, [N] is not always appearing...rdkit uses [N], [NH], [NH2+], or [N.]
# 
# n_atom.SetAtomicNum(7)
# n_atom.SetFormalCharge(0)
# n_atom.SetNumExplicitHs(0)
# n_atom.SetNoImplicit(True)
#
# So, use the code with Sbridge and replace S with N
#

#!/usr/bin/env python3

from rdkit import Chem
from rdkit.Chem import rdchem

input_file = "QOOH_sp3.smi"
output_file = "TS_Nbridge.smi"


def find_o_peroxy_and_h_ooh(molH):
    """
    In a QOOH, find the hydroperoxide O atom:
      - O bonded to another O and to H
    Return (o_idx, h_idx).
    """
    for atom in molH.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        neigh = atom.GetNeighbors()
        h_nei = [n for n in neigh if n.GetSymbol() == "H"]
        o_nei = [n for n in neigh if n.GetSymbol() == "O"]
        if len(h_nei) == 1 and len(o_nei) >= 1:
            return atom.GetIdx(), h_nei[0].GetIdx()
    return None, None


def find_donor_carbon(molH):
    """
    Find the unique carbon with one radical electron.
    """
    cand = [
        a.GetIdx() for a in molH.GetAtoms()
        if a.GetSymbol() == "C" and a.GetNumRadicalElectrons() == 1
    ]
    if len(cand) != 1:
        return None
    return cand[0]


def get_n_bridge_ring_size(mol, n_idx, o_idx, c_idx):
    """
    Find ring(s) containing N, peroxy O, and donor C.
    Return the smallest such ring size.
    """
    ring_info = mol.GetRingInfo()
    sizes = []
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        if n_idx in ring_set and o_idx in ring_set and c_idx in ring_set:
            sizes.append(len(ring))
    return min(sizes) if sizes else None


def get_atom_idx_by_mapnum(mol, mapnum):
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == mapnum:
            return atom.GetIdx()
    return None


def clear_all_atom_maps(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)


def clear_all_atom_chiral_tags(mol):
    for atom in mol.GetAtoms():
        atom.SetChiralTag(rdchem.ChiralType.CHI_UNSPECIFIED)


def assign_and_label_cip(mol):
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    return mol


def smiles_without_maps(mol):
    m = Chem.Mol(mol)
    clear_all_atom_maps(m)
    return Chem.MolToSmiles(m, isomericSmiles=True)


def make_n_bridge_template_with_maps(mol):
    """
    Build the N-bridge template from the QOOH radical.

    Steps:
      - add explicit Hs
      - convert OOH hydrogen -> N
      - connect N to donor radical carbon
      - keep all existing stereochemistry fixed
      - mark donor carbon, N, and peroxy O with atom-map numbers

    Returns:
      (mol_noH, donor_mapnum, n_mapnum, o_mapnum, ring_size)
    or
      (None, None, None, None, None) on failure.
    """
    molH = Chem.AddHs(mol)

    o_idx, h_idx = find_o_peroxy_and_h_ooh(molH)
    if o_idx is None:
        return None, None, None, None, None

    c_idx = find_donor_carbon(molH)
    if c_idx is None:
        return None, None, None, None, None

    rw = Chem.RWMol(molH)

    # Convert hydroperoxide H into N
    n_atom = rw.GetAtomWithIdx(h_idx)
    n_atom.SetAtomicNum(7)
    n_atom.SetFormalCharge(0)
    n_atom.SetNumExplicitHs(0)
    n_atom.SetNoImplicit(True)

    # Ensure N-O and N-C bonds exist
    if rw.GetBondBetweenAtoms(h_idx, o_idx) is None:
        rw.AddBond(h_idx, o_idx, rdchem.BondType.SINGLE)
    if rw.GetBondBetweenAtoms(h_idx, c_idx) is None:
        rw.AddBond(h_idx, c_idx, rdchem.BondType.SINGLE)

    try:
        m_newH = rw.GetMol()
        Chem.SanitizeMol(m_newH)
        Chem.AssignStereochemistry(m_newH, force=True, cleanIt=True)
    except Exception:
        return None, None, None, None, None

    # Atom-map numbers so indices can be recovered safely after RemoveHs()
    donor_mapnum = 901
    n_mapnum = 902
    o_mapnum = 903
    m_newH.GetAtomWithIdx(c_idx).SetAtomMapNum(donor_mapnum)
    m_newH.GetAtomWithIdx(h_idx).SetAtomMapNum(n_mapnum)
    m_newH.GetAtomWithIdx(o_idx).SetAtomMapNum(o_mapnum)

    ring_size = get_n_bridge_ring_size(m_newH, h_idx, o_idx, c_idx)

    try:
        m_new = Chem.RemoveHs(m_newH)
        Chem.AssignStereochemistry(m_new, force=True, cleanIt=True)
    except Exception:
        return None, None, None, None, None

    return m_new, donor_mapnum, n_mapnum, o_mapnum, ring_size


def clone_with_donor_chirality(base_mol, donor_idx, donor_tag):
    """
    Make one copy with donor carbon explicitly set to one tetrahedral sense.
    Existing substrate stereochemistry remains fixed.
    """
    m = Chem.Mol(base_mol)
    atom = m.GetAtomWithIdx(donor_idx)
    atom.SetChiralTag(donor_tag)

    try:
        assign_and_label_cip(m)
        Chem.SanitizeMol(m)
        assign_and_label_cip(m)
    except Exception:
        return None

    return m


def canonical_smiles_same_stereo(mol):
    """
    Canonical isomeric SMILES with atom maps removed.
    """
    m = Chem.Mol(mol)
    clear_all_atom_maps(m)
    assign_and_label_cip(m)
    return Chem.MolToSmiles(m, isomericSmiles=True)


def canonical_smiles_no_stereo(mol):
    """
    Canonical connectivity-only SMILES with atom maps removed.
    """
    m = Chem.Mol(mol)
    clear_all_atom_maps(m)
    clear_all_atom_chiral_tags(m)
    return Chem.MolToSmiles(m, isomericSmiles=False)


def nbridge_to_qooh(molN):
    """
    Convert N-bridge scaffold back to QOOH by:
      - changing N -> H
      - removing N-C bond
      - keeping only the O-H bond

    Returns a sanitized heavy-atom molecule without explicit Hs.
    """
    m = Chem.Mol(molN)
    rw = Chem.RWMol(m)

    n_idxs = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "N"]
    if len(n_idxs) != 1:
        return None
    n_idx = n_idxs[0]

    n_atom = rw.GetAtomWithIdx(n_idx)
    nbrs = [n.GetIdx() for n in n_atom.GetNeighbors()]
    o_nbrs = [i for i in nbrs if rw.GetAtomWithIdx(i).GetSymbol() == "O"]
    c_nbrs = [i for i in nbrs if rw.GetAtomWithIdx(i).GetSymbol() == "C"]

    if len(o_nbrs) != 1 or len(c_nbrs) != 1:
        return None

    o_idx = o_nbrs[0]
    c_idx = c_nbrs[0]

    bond_nc = rw.GetBondBetweenAtoms(n_idx, c_idx)
    if bond_nc is None:
        return None
    rw.RemoveBond(n_idx, c_idx)

    n_atom = rw.GetAtomWithIdx(n_idx)
    n_atom.SetAtomicNum(1)
    n_atom.SetFormalCharge(0)
    n_atom.SetNumExplicitHs(0)
    n_atom.SetNoImplicit(True)

    try:
        mH = rw.GetMol()
        Chem.SanitizeMol(mH)
        m = Chem.RemoveHs(Chem.AddHs(mH))
        Chem.AssignStereochemistry(m, force=True, cleanIt=True)
        return m
    except Exception:
        return None


def nbridge_to_roo(molN):
    """
    Convert N-bridge scaffold to ROO radical by:
      - removing the N atom entirely
      - removing both N-O and N-C bonds
      - adding one radical electron to the donor carbon

    Returns a sanitized heavy-atom molecule without explicit Hs.
    """
    m = Chem.Mol(molN)
    rw = Chem.RWMol(m)

    n_idxs = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == "N"]
    if len(n_idxs) != 1:
        return None
    n_idx = n_idxs[0]

    n_atom = rw.GetAtomWithIdx(n_idx)
    nbrs = [n.GetIdx() for n in n_atom.GetNeighbors()]
    o_nbrs = [i for i in nbrs if rw.GetAtomWithIdx(i).GetSymbol() == "O"]
    c_nbrs = [i for i in nbrs if rw.GetAtomWithIdx(i).GetSymbol() == "C"]

    if len(o_nbrs) != 1 or len(c_nbrs) != 1:
        return None

    c_idx = c_nbrs[0]

    # Remove N atom after noting donor carbon
    rw.RemoveAtom(n_idx)

    # Adjust donor carbon index if needed after deletion
    if c_idx > n_idx:
        c_idx -= 1

    c_atom = rw.GetAtomWithIdx(c_idx)
    c_atom.SetNumRadicalElectrons(1)
    c_atom.SetNoImplicit(True)

    try:
        m0 = rw.GetMol()
        Chem.SanitizeMol(m0)
        m = Chem.RemoveHs(Chem.AddHs(m0))
        Chem.AssignStereochemistry(m, force=True, cleanIt=True)
        return m
    except Exception:
        return None


def endpoints_match(iso1, iso2):
    """
    Keep only TS pairs for which BOTH endpoint checks pass:
      - same QOOH endpoint
      - same ROO endpoint

    Uses stereochemistry-aware canonical isomeric SMILES.
    """
    q1 = nbridge_to_qooh(iso1)
    q2 = nbridge_to_qooh(iso2)
    if q1 is None or q2 is None:
        return False, None, None, None, None

    r1 = nbridge_to_roo(iso1)
    r2 = nbridge_to_roo(iso2)
    if r1 is None or r2 is None:
        return False, None, None, None, None

    q1_smi = canonical_smiles_same_stereo(q1)
    q2_smi = canonical_smiles_same_stereo(q2)
    r1_smi = canonical_smiles_same_stereo(r1)
    r2_smi = canonical_smiles_same_stereo(r2)

    ok = (q1_smi == q2_smi) and (r1_smi == r2_smi)
    return ok, q1_smi, q2_smi, r1_smi, r2_smi


n_total = 0
n_written = 0
n_skipped = 0
n_not_exact2 = 0
n_endpoint_fail = 0

with open(input_file) as infile, open(output_file, "w") as out:

    for line in infile:
        line = line.strip()

        parts = line.split()
        if len(parts) < 2:
            continue

        smiles, mol_id = parts[0], parts[1]
        n_total += 1

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            n_skipped += 1
            continue

        base_N, donor_mapnum, n_mapnum, o_mapnum, ring_size = make_n_bridge_template_with_maps(mol)
        if base_N is None:
            n_skipped += 1
            continue

        donor_idx = get_atom_idx_by_mapnum(base_N, donor_mapnum)
        n_idx = get_atom_idx_by_mapnum(base_N, n_mapnum)

        if donor_idx is None or n_idx is None:
            n_skipped += 1
            continue

        donor_atom = base_N.GetAtomWithIdx(donor_idx)

        # Require donor carbon to be a potential tetrahedral stereocenter
        if donor_atom.GetDegree() < 3:
            n_skipped += 1
            continue

        # Build exactly two versions by explicitly setting donor chirality
        iso1 = clone_with_donor_chirality(
            base_N, donor_idx, rdchem.ChiralType.CHI_TETRAHEDRAL_CW
        )
        iso2 = clone_with_donor_chirality(
            base_N, donor_idx, rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
        )

        if iso1 is None or iso2 is None:
            n_skipped += 1
            continue

        smi1 = smiles_without_maps(iso1)
        smi2 = smiles_without_maps(iso2)

        # Collapse check: not actually two unique N-bridge stereomers
        if smi1 == smi2:
            n_not_exact2 += 1
            print(f"{mol_id}: donor stereoflip collapsed to one unique isomer; skipped")
            continue

        # Endpoint-consistency filter:
        # keep only pairs that map back to the SAME ROO and SAME QOOH endpoints
        ok, q1_smi, q2_smi, r1_smi, r2_smi = endpoints_match(iso1, iso2)
        if not ok:
            n_endpoint_fail += 1
            print(f"{mol_id}: endpoint mismatch; skipped")
            if q1_smi is not None and q2_smi is not None:
                print(f"  QOOH1: {q1_smi}")
                print(f"  QOOH2: {q2_smi}")
            if r1_smi is not None and r2_smi is not None:
                print(f"  ROO1 : {r1_smi}")
                print(f"  ROO2 : {r2_smi}")
            continue

        # Canonicalize ordering for reproducible _01/_02 labels
        smis = sorted([smi1, smi2])

        for i, smi_out in enumerate(smis, start=1):
            new_id = f"{mol_id}_Nbridge_{i:02d}"
            out.write(f"{smi_out},{new_id},{ring_size}-mem\n")
            if ring_size is not None:
                print(f"{new_id}: new N-ring is {ring_size}-membered")
            else:
                print(f"{new_id}: could not determine N-ring size")
            n_written += 1

print("Total QOOH entries read:", n_total)
print("Exact two N-bridge stereomers written:", n_written)
print("Skipped (pattern fail / sanitize / invalid):", n_skipped)
print("Collapsed to <2 unique donor stereomers:", n_not_exact2)
print("Failed same-endpoint filter:", n_endpoint_fail)
print("Output file:", output_file)
