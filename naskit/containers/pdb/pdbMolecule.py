from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom
from ...exceptions import InvalidPDB



class PdbMolecule:
    __slots__ = ("__atoms", "__name_idx_map")
    
    def __init__(self):
        self.__atoms = []
        self.__name_idx_map = {}
        
    def _remap(self):
        self.__name_idx_map.clear()
        self.__name_idx_map = {atom.name:i for i, atom in enumerate(self.__atoms)}
        
        
    def __getitem__(self, i: int):
        return self.__atoms[i]
    
    def __len__(self):
        return len(self.__atoms)
    
    def __iter__(self):
        return iter(self.__atoms)
    
    def __str__(self):
        return "\n".join([str(a) for a in self.__atoms])
    
        
    def add_atom(self, atom: PdbAtom):
        if len(self.__atoms):
            if self.__name_idx_map.get(atom.name) is not None:
                raise InvalidPDB(f"Atom ({atom.atomn}) with name {atom.name} "
                                 f"already exists in molecule.")
            
            a = self.__atoms[0]
            if atom.mol_name!=a.mol_name:
                raise InvalidPDB(f"All atoms of a molecule (number {a.moln}) "
                                 f"must have the same molecule name ({a.mol_name}), "
                                 f"got {atom.mol_name}.")

            if atom.chain!=a.chain:
                raise InvalidPDB(f"All atoms of a molecule (number {a.moln}) "
                                 f"must have the same chain name ({a.chain}), "
                                 f"got {atom.chain}.")
            
        self.__atoms.append(atom)
        self.__name_idx_map[atom.name] = len(self.__atoms)
        
    def get_atom_idx(self, name: str):
        return self.__name_idx_map.get(name)
    
    def get_atom(self, i: Union[int, str]):
        if isinstance(i, int):
            return self.__atoms[i]
        elif isinstance(i, str):
            return self.__atoms[self.__name_idx_map[i]]
        else:
            raise IndexError(f"Invalid argument of type {type(i)}, accepted: int index or str name.")
            
    def delete_atom(self, i: Union[int, str]):
        if isinstance(i, int):
            return self.__atoms.pop(i)
        elif isinstance(i, str):
            return self.__atoms.pop(self.__name_idx_map[i])
        else:
            raise IndexError(f"Invalid argument of type {type(i)}, accepted: int index or str name.")
        self._remap()
        
        
    def renum_atoms(self, initn: int = 1):
        for i, a in enumerate(self.__atoms):
            a.atomn = initn + i
        
    
    @property
    def natoms(self):
        return len(self)
    
    @property
    def moln(self):
        return self.__atoms[0].moln
    
    @moln.setter
    def moln(self, moln: int):
        for a in self.__atoms:
            a.moln = moln
            
    @property
    def name(self):
        return self.__atoms[0].mol_name
    
    @name.setter
    def name(self, name: str):
        for a in self.__atoms:
            a.mol_name = name
    
    @property
    def chain(self):
        return self.__atoms[0].chain
    
    @chain.setter
    def chain(self, chain_name: str):
        if (not chain_name.isascii()) or (not chain_name.isupper() or len(chain_name)!=1):
            raise ValueError(f"Chain name must be a single upper ascii character, got {chain_name}.")
            
        for a in self.__atoms:
            a.chain = chain_name
    
    @property
    def coords(self):
        return np.stack([a.coords for a in self.__atoms])
    
    @coords.setter
    def coords(self, coords: np.ndarray):
        if coords.shape[0]!=len(self) or coords.shape[1]!=3:
            raise ValueError(f"Coords matrix must have shape: ({len(self)}, 3), got {coords.shape}")
            
        for i in range(coords.shape[0]):
            self.__atoms[i].coords = coords[i]
            
            
class PdbResidue(PdbMolecule):
    def __init__(self):
        super().__init__()
        
        
class NucleicAcidResidue(PdbResidue):
    def __init__(self):
        super().__init__()
        
        
class AminoacidResidue(PdbResidue):
    def __init__(self):
        super().__init__()