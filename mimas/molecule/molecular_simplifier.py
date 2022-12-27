#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from .carbon_chain_finder import find_carbon_chain


class MolecularSimplifier:
    """
    Simplify a molecular structure by removing carbon chains.
    """

    def __init__(self, mol) -> None:
        self.mol = mol

    def get_simplified_mol(self) -> Chem.Mol:
        """
        Get the simplified molecular structure.
        """
        # 1. Find the carbon chains.
        carbon_chain_list = find_carbon_chain(self.mol, max_break_bond=2)

        # 2. Replace the carbon chains with a single carbon atom.
        atoms_to_remove = []
        mol_new = Chem.RWMol(self.mol)
        for carbon_chain_old, break_bond_list in carbon_chain_list:
            if len(carbon_chain_old) == len(break_bond_list):
                continue
            # Generate a new carbon atom.
            new_carbon_chain = []
            atoms_to_remove += list(carbon_chain_old)
            for atom_origin_a, atom_origin_b in break_bond_list:
                atom_new = mol_new.AddAtom(Chem.Atom('C'))
                if len(new_carbon_chain) > 0:
                    # Add a bond between the new carbon atom and the previous carbon atom.
                    bond_new = mol_new.AddBond(new_carbon_chain[-1], atom_new, order=Chem.rdchem.BondType.SINGLE)
                new_carbon_chain.append(atom_new)
                # Attach the new carbon atom to the molecule.
                mol_new.RemoveBond(atom_origin_a, atom_origin_b)
                mol_new.AddBond(atom_origin_a, atom_new, order=Chem.rdchem.BondType.SINGLE)
        atoms_to_remove = [mol_new.GetAtomWithIdx(i) for i in atoms_to_remove]
        [mol_new.RemoveAtom(atom.GetIdx()) for atom in atoms_to_remove]

        # 3. Return the simplified molecular structure.
        return mol_new.GetMol()


if __name__ == '__main__':
    mol_smiles = "NCCOP(O)(=O)OCC(COC=CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCC=CCC=C(O)CCCCC"
    mol = Chem.MolFromSmiles(mol_smiles)
    ms = MolecularSimplifier(mol)
    simplified_mol = ms.get_simplified_mol()
    print(Chem.MolToSmiles(simplified_mol))
