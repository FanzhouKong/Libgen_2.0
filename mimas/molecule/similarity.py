from rdkit import Chem
from rdkit.Chem import rdFMCS


def __remove_bonds_in_smarts(mol, smarts):
    removed_bond_number = 0
    sub_atom = mol.GetSubstructMatch(Chem.MolFromSmarts(smarts))
    # print(sub_atom)
    for i in sub_atom:
        all_bonds = mol.GetAtomWithIdx(i).GetBonds()
        for bond in all_bonds:
            atom_1, atom_2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            mol.RemoveBond(atom_1, atom_2)
            if atom_1 in sub_atom and atom_2 in sub_atom:
                removed_bond_number += 1
    return mol, removed_bond_number


def calculate_similarity(mol1, mol2):
    mol1 = Chem.rdchem.RWMol(mol1)
    mol2 = Chem.rdchem.RWMol(mol2)
    bond_number_common = 0
    bond_number_mol1 = len(mol1.GetBonds())
    bond_number_mol2 = len(mol2.GetBonds())

    while True:
        res = rdFMCS.FindMCS([mol1, mol2], timeout=60, threshold=1,
                             ringMatchesRingOnly=True, completeRingsOnly=True,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                             ringCompare=rdFMCS.RingCompare.StrictRingFusion,
                             maximizeBonds=True, matchValences=True)
        if res.numBonds == 0:
            break

        common_s = res.smartsString
        # print(common_s)

        mol1, _ = __remove_bonds_in_smarts(mol1, common_s)
        mol2, _ = __remove_bonds_in_smarts(mol2, common_s)
        bond_number_common += res.numBonds
        # print(bond_num)

    return {
        "mol1_bond_number": bond_number_mol1,
        "mol2_bond_number": bond_number_mol2,
        "common_bond_number": bond_number_common,
        "similarity": (2 * bond_number_common) / (bond_number_mol1 + bond_number_mol2)
    }


if __name__ == '__main__':
    import pprint

    s1 = "InChI=1S/C10H7NO2/c12-10(13)9-6-5-7-3-1-2-4-8(7)11-9/h1-6H,(H,12,13)"
    s2 = "InChI=1S/C10H7NO2/c12-10(13)9-5-7-3-1-2-4-8(7)6-11-9/h1-6H,(H,12,13)"
    mol_a = Chem.MolFromInchi(s1)
    mol_b = Chem.MolFromInchi(s2)
    result = calculate_similarity(mol_a, mol_b)
    pprint.pprint(result)

    s1 = "CC(CCC)c1ccccc1"
    s2 = "CCC(CC)c1ccccc1"
    mol_a = Chem.rdchem.RWMol(Chem.MolFromSmiles(s1))
    mol_b = Chem.rdchem.RWMol(Chem.MolFromSmiles(s2))
    result = calculate_similarity(mol_a, mol_b)
    pprint.pprint(result)
