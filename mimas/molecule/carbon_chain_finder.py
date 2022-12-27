import numpy as np
import math
import itertools
import copy

import rdkit
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.rdchem import Mol, RWMol, BondType
import rdkit.Chem.rdmolops


def find_carbon_chain(mol=None, smiles=None, max_break_bond=2):
    """
    Find carbon chain in a molecule with the following criteria:

    1. Connect with rest part with one or two bonds.
    2. It contains at least one atom.
    3. All the atoms in the fragments are carbon.

    :param mol: RDKit molecule
    :param smiles: SMILES string
    :return: A list of carbon chain information. The list is sorted by the length of carbon chain. \
    Every carbon chain information contains two lists: the first list is the index of atoms in the carbon chain, like: (C1, C2, ..., CN). \
    The second list is the breaked bonds to generate the carbon chain, like: [(A1, B1), (A2, B2), ..., (AN, BN)]. The second atom index BN is always in the carbon chain.
    """
    if mol is None:
        mol = Chem.MolFromSmiles(smiles)

    # Get carbon atom index
    atom_c_idx = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C"}

    # Only consider atoms connected with simple or double bonds.
    all_not_single_or_double_bonds = [bond for bond in mol.GetBonds()
                                      if bond.GetBondType() not in {Chem.BondType.SINGLE, Chem.BondType.DOUBLE}]
    exclude_atom_idx = {bond.GetBeginAtomIdx() for bond in all_not_single_or_double_bonds} | \
        {bond.GetEndAtomIdx() for bond in all_not_single_or_double_bonds}
    atom_c_idx = atom_c_idx-exclude_atom_idx

    # Only consider atoms not in the ring
    rdkit.Chem.rdmolops.FastFindRings(mol)
    ring = mol.GetRingInfo()
    exclude_atom_idx = {y for x in ring.AtomRings() for y in x}
    atom_c_idx = atom_c_idx-exclude_atom_idx

    # Get all bond list
    c_c_bond_list = _find_c_c_bond_list(mol, atom_c_idx)

    # Generate all carbon chain list
    carbon_chain_all = _generate_all_carbon_chain_list(mol, atom_c_idx, c_c_bond_list, max_break_bond=max_break_bond)

    # Refine the carbon chain list.
    carbon_chain_refined = _refine_carbon_chain_list(carbon_chain_all)

    return carbon_chain_refined


def _refine_carbon_chain_list(carbon_chain_all):
    carbon_chain_all.sort(key=lambda x: len(x[0]), reverse=True)
    carbon_chain_refined = []
    for cur_cc in carbon_chain_all:
        carbon_chain, _ = cur_cc
        for known_carbon_chain, _ in carbon_chain_refined:
            if carbon_chain.issubset(known_carbon_chain):
                break
        else:
            carbon_chain_refined.append(cur_cc)
    return carbon_chain_refined


def _generate_all_carbon_chain_list(mol, atom_c_idx, c_c_bond_list, max_break_bond):
    carbon_chain_list = []
    for break_num in range(1, max_break_bond+1):
        for break_bonds in itertools.combinations(c_c_bond_list, break_num):
            mol_tmp = RWMol(mol)
            for bond in break_bonds:
                mol_tmp.RemoveBond(*bond)
            all_fragments = list(Chem.GetMolFrags(mol_tmp))
            if len(all_fragments) == 1:
                continue
            elif len(all_fragments) > 2:
                all_fragments.sort(key=lambda x: -len(x))
            for frag in all_fragments:
                frag = set(frag)
                if frag == (frag & atom_c_idx):
                    final_bonds = []
                    is_wanted = True
                    # This is wanted fragmentation.
                    for bond in break_bonds:
                        if bond[0] in frag:
                            final_bonds.append((bond[1], bond[0]))
                        elif bond[1] in frag:
                            final_bonds.append((bond[0], bond[1]))
                        else:
                            is_wanted = False
                            break
                    if is_wanted:
                        # Sort the bond list to make the fist bond is in the largest fragment.
                        if len(final_bonds) == 2:
                            if len(frag) <= 1:
                                continue
                            if final_bonds[0][0] not in all_fragments[0]:
                                final_bonds = [final_bonds[1], final_bonds[0]]

                        # Record this carbon chain
                        carbon_chain_list.append([frag, tuple(final_bonds)])
    return carbon_chain_list


def _find_c_c_bond_list(mol, atom_c_idx):
    c_c_bond_list = []
    for bond in mol.GetBonds():
        # 1. Only continue with single bond
        if not(bond.GetBondType() == Chem.BondType.SINGLE):
            continue

        atom_1_idx, atom_2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if atom_1_idx not in atom_c_idx or atom_2_idx not in atom_c_idx:
            continue

        c_c_bond_list.append((atom_1_idx, atom_2_idx))
    return c_c_bond_list


def test():
    mol_smiles = "OP(O)(=O)OCCC(NC(=O)CCC)C(O)\C=C\C"
    print(find_carbon_chain(smiles=mol_smiles))
    return
    peaks = np.array(
        [[226.083, 100], [73.064, 200], [79.0, 999]],
        dtype=np.float32)
    lipid_spec.add_spectrum(peaks, '[M+NH4]+', ms2_tolerance_da=0.05)
    matches = lipid_spec.match_spectrum(
        327.1679,
        np.array(
            [[226.083, 100], [87.080, 200], [79.0, 999]],
            dtype=np.float32),
        '[M+NH4]+', ms1_tolerance_ppm=10., ms2_tolerance_da=0.05)

    print(matches)
    return


if __name__ == '__main__':
    test()
