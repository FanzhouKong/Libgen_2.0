from rdkit import Chem
from rdkit.Chem.rdchem import RWMol
import re
import itertools


class Substructure:
    def __init__(self, mol=None, smiles=None):
        if mol is not None:
            self.mol = mol
        else:
            self.mol = Chem.MolFromSmiles(smiles)
        self._smiles = ""
        self._inchikey = ""

    def break_bond_num(self):
        n = len([m for m in re.finditer("\[13CH*\d*\]", self.smiles)])
        return n

    def heavy_atom_num(self):
        return self.mol.GetNumHeavyAtoms() - self.break_bond_num()

    @property
    def smiles(self):
        if not self._smiles:
            self._smiles = Chem.MolToSmiles(self.mol)
        return self._smiles

    @property
    def smiles_for_output(self):
        smiles = Chem.MolToSmiles(self.mol)
        smiles = re.sub(r"\[13CH*\d*\]", "[*]", smiles)
        return smiles

    @property
    def inchikey(self):
        if not self._inchikey:
            self._inchikey = Chem.MolToInchiKey(self.mol)
        return self._inchikey

    def __getstate__(self):
        return self.smiles

    def __setstate__(self, state):
        mol = Chem.MolFromSmiles(state)
        self.__init__(mol)

    def __hash__(self):
        return hash(self.inchikey)

    def __str__(self):
        return self.smiles

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Substructure):
            return self.inchikey == __o.inchikey
        else:
            return False

    def get_mol_for_matching(self):
        smiles = re.sub(r"\[13CH*\d*\]", "[*]", self.smiles)
        mol = Chem.MolFromSmiles(smiles)
        return mol

    def reset_mol(self, mol=None):
        if mol is None:
            mol = self.mol
        self.__init__(mol)

    def is_all_c(self):
        smiles_no_c = re.sub("\[13CH\d\]", "", self.smiles)
        smiles_no_c = \
            smiles_no_c.replace("C", ""). \
            replace("(", "").replace(")", ""). \
            replace("=", "")
        return len(smiles_no_c) == 0

    def neutralize_atoms(self):
        # Remove charge
        if self.smiles.find('+') >= 0 or self.smiles.find('-') >= 0:
            pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
            at_matches = self.mol.GetSubstructMatches(pattern)
            at_matches_list = [y[0] for y in at_matches]
            if len(at_matches_list) > 0:
                for at_idx in at_matches_list:
                    atom = self.mol.GetAtomWithIdx(at_idx)
                    chg = atom.GetFormalCharge()
                    hcount = atom.GetTotalNumHs()
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(hcount - chg)
                    atom.UpdatePropertyCache()
            self.reset_mol()

    def find_substructures(self, break_num_list):
        finder = SubstructureFinder(self.mol)
        finder.find_possible_bond_break()
        all_ss = {}
        for b in break_num_list:
            ss = finder.find_substructures_with_bond_break(b)
            all_ss.update(ss)
        return all_ss


class SubstructureFinder:
    def __init__(self, mol):
        self.mol = mol
        self.critical_bonds = []

    def _is_neighbor_atom_c(self, atom_idx):
        atom = self.mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != "C":
            return False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() != "C":
                return False
        return True

    def find_possible_bond_break(self):
        # find all bonds
        all_bonds = []
        for bond in self.mol.GetBonds():
            bond_type = bond.GetBondType()
            if bond_type != Chem.rdchem.BondType.SINGLE:
                continue
            atom_1, atom_2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            all_bonds.append([atom_1, atom_2, bond_type])

        # Remove bond one by one
        for bond in all_bonds:
            mol = RWMol(self.mol)
            atom_1, atom_2, bond_type = bond
            mol.RemoveBond(atom_1, atom_2)
            mol_frags = Chem.GetMolFrags(mol)
            if len(mol_frags) == 2:
                # Record this bond
                self.critical_bonds.append(bond)

    def find_substructures_with_bond_break(self, break_num):
        """
        Find all possible substructures with `break_num` bond break. 
        Output the substructures, the corresponding atom index, and the breaked bond.
        The first atom in the bond will be in the original molecule, the second atom in the bond will be in the substructure, the third entry will be the breaked bond type.

        :param break_num: int
        :return: dict: {inchikey: Substructure,{atom index},[breaked bond]}
        """
        result = {}
        for cur_break_bonds in itertools.combinations(self.critical_bonds, break_num):
            mol = RWMol(self.mol)
            all_atoms_number = [x.GetIdx() for x in mol.GetAtoms()]
            all_atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
            for bond in cur_break_bonds:
                atom_1, atom_2, bond_type = bond
                mol.RemoveBond(atom_1, atom_2)
            all_mol_frags = Chem.GetMolFrags(mol)

            for frag in all_mol_frags:
                for bond in cur_break_bonds:
                    atom_1, atom_2, bond_type = bond
                    if atom_1 not in frag and atom_2 not in frag:
                        break
                else:
                    if mol is None:
                        mol = RWMol(self.mol)
                        all_atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
                        for bond in cur_break_bonds:
                            atom_1, atom_2, bond_type = bond
                            mol.RemoveBond(atom_1, atom_2)

                    bonds_output = []
                    for bond in cur_break_bonds:
                        atom_1, atom_2, bond_type = bond
                        if atom_1 in frag:
                            atom_stay = atom_1
                            bond = atom_2, atom_1, bond_type
                            bonds_output.append(bond)
                        elif atom_2 in frag:
                            atom_stay = atom_2
                            bond = atom_1, atom_2, bond_type
                            bonds_output.append(bond)
                        else:
                            raise RuntimeError()

                        if bond_type == Chem.BondType.SINGLE:
                            atom_add = mol.AddAtom(Chem.AtomFromSmiles("[13CH3]"))
                        elif bond_type == Chem.BondType.DOUBLE:
                            atom_add = mol.AddAtom(Chem.AtomFromSmiles("[13CH2]"))
                        elif bond_type == Chem.BondType.TRIPLE:
                            atom_add = mol.AddAtom(Chem.AtomFromSmiles("[13CH]"))
                        else:
                            atom_add = mol.AddAtom(Chem.AtomFromSmiles("[13C]"))
                        # print(Chem.MolToSmiles(mol))
                        mol.AddBond(atom_add, atom_stay, bond_type)
                        # print("#", Chem.MolToSmiles(mol))
                    for a in all_atoms_number:
                        if a not in frag:
                            mol.RemoveAtom(all_atoms[a].GetIdx())

                    if mol is not None:
                        ss = Substructure(mol)
                        ss.neutralize_atoms()
                        result[ss.inchikey] = ss, set(frag), bonds_output
                    mol = None

        return result
