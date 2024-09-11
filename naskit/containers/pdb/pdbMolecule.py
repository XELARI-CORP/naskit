from typing import Union, List, Iterable, Tuple
import numpy as np
from .pdbAtom import PdbAtom
from .pdbDraw import PDBDraw
from ...exceptions import InvalidPDB
from ...utils.math3d import align



class PdbMolecule(PDBDraw):
    __slots__ = ("__atoms", "__name_idx_map")
    
    def __init__(self):
        self.__atoms = []
        self.__name_idx_map = {}
        
    def _remap(self):
        self.__name_idx_map.clear()
        self.__name_idx_map = {atom.aname:i for i, atom in enumerate(self.__atoms)}
        
        
    def __getitem__(self, i: Union[int, str, list, tuple]):
        if isinstance(i, int):
            return self.__atoms[i]
            
        elif isinstance(i, str):
            return self.__atoms[self.__name_idx_map[i]]

        elif isinstance(i, (list, tuple)):
            if not all([isinstance(j, (int, str)) for j in i]):
                raise IndexError(f"Expected all indices to be int index or string name.")

            submol = self.__class__()
            for j in i:
                PdbMolecule.add_atom(submol, self[j])
            
            return submol
            
        else:
            raise IndexError(f"Expected int index, string name or lits/tuple of these types got {type(i)}.")

    
    def __len__(self):
        return len(self.__atoms)

    def __contains__(self, atom_name: str):
        return self.__name_idx_map.get(atom_name) is not None
    
    def __iter__(self):
        return iter(self.__atoms)
    
    def __repr__(self):
        return f"{self.mname} {self.__class__.__name__} with {len(self)} atoms at {hex(id(self))}"
    
    def __str__(self):
        return "\n".join([str(a) for a in self.__atoms])
    
        
    def copy(self):
        copied_mol = self.__class__()
        for a in self.__atoms:
            copied_mol.add_atom(a.copy(), skip_validation=True)
        
        return copied_mol

    def dist(self, a: Union["PdbAtom", "PdbMolecule"]):
        if len(a.coords.shape)==1:
            return np.linalg.norm((self.coords - a.coords), axis=1)
        return np.linalg.norm((self.coords[:, np.newaxis, :] - a.coords), axis=2)

    def translate(self, lang: str = "amber"):
        for a in self.__atoms:
            a.translate(lang)
        
    def add_atom(self, atom: PdbAtom, skip_validation: bool = False):
        if len(self.__atoms) and (not skip_validation):
            if self.__name_idx_map.get(atom.aname) is not None:
                raise InvalidPDB(f"Atom (number {atom.anum}) with name {atom.aname} "
                                 f"is already in molecule {self.mname} (number {self.mnum}).")
            
            a = self.__atoms[0]
            if atom.mname!=a.mname:
                raise InvalidPDB(f"All atoms of a molecule (number {a.mnum}) "
                                 f"must have the same molecule name ({a.mname}), "
                                 f"got {atom.mname}.")

            if atom.chain!=a.chain:
                raise InvalidPDB(f"All atoms of a molecule (number {a.mnum}) "
                                 f"must have the same chain name ({a.chain}), "
                                 f"got {atom.chain}.")
            
        self.__atoms.append(atom)
        self.__name_idx_map[atom.aname] = len(self.__atoms) - 1
        
    def get_atom_idx(self, name: str):
        return self.__name_idx_map.get(name)
    
    def delete_atom(self, i: Union[int, str]):
        if isinstance(i, int):
            self.__atoms.pop(i)
        elif isinstance(i, str):
            self.__atoms.pop(self.__name_idx_map[i])
        else:
            raise IndexError(f"Invalid argument of type {type(i)}, accepted: int index or str name.")
        self._remap()
        
    def renum_atoms(self, initn: int = 1):
        for i, a in enumerate(self.__atoms):
            a.anum = initn + i
        
    def atoms(self):
        for a in self.__atoms:
            yield a
            
            
    @property
    def natoms(self):
        return len(self)
    
    @property
    def mnum(self):
        return self.__atoms[0].mnum
    
    @mnum.setter
    def mnum(self, mnum: int):
        for a in self.__atoms:
            a.mnum = mnum
            
    @property
    def mname(self):
        return self.__atoms[0].mname
    
    @mname.setter
    def mname(self, mname: str):
        for a in self.__atoms:
            a.mname = mname
    
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
            
            
    def embed_molecule_fragment(self, 
                                other: Union["PdbMolecule", "AminoacidResidue", "NucleicAcidResidue"],
                                source_atoms: Iterable[str],
                                embed_atoms: Iterable[str],
                                correspondence: Iterable[Tuple[str, str]],
                                source_origin_atom: str = None,
                                embed_origin_atom: str = None
                               ):
        """
        Substitutes source atoms of molecule to embed atoms of other molecule.
        Other molecule is aligned by pairs of atoms provided in correspondence.
        
        :param other: aligned structure
        :param source_atoms: atoms to be deleted
        :param embed_atoms: atoms to be embedded
        :param correspondence: pairs of indices (ai, bi) for aligning, ai - changed molecule, bi - embedded
        :param source_origin_atom: name of source molecule's atom which coordinates will be used for origin shift
        :param embed_origin_atom: name of embed molecule's atom which coordinates will be used for origin shift
        """
        other = other.copy()
        
        aidx, bidx = [], []
        for a1, a2 in correspondence:
            if (i:=self.get_atom_idx(a1)) is None:
                raise KeyError(f"Target molecule does not have {a1} atom.")
            if (j:=other.get_atom_idx(a2)) is None:
                raise KeyError(f"Embedded fragment does not have {a2} atom.")
            
            aidx.append(i)
            bidx.append(j)
            
        source_origin_idx = None if (source_origin_atom is None) else self.get_atom_idx(source_origin_atom)
        embed_origin_idx = None if (embed_origin_atom is None) else other.get_atom_idx(embed_origin_atom)
        
        other.coords = align(self.coords, other.coords, 
                             aidx, bidx,
                             source_origin_idx, embed_origin_idx)
        other.mnum = self.mnum
        other.mname = self.mname
        other.chain = self.chain
        
        for aname in source_atoms:
            if aname not in self:
                if aname[0]=="H":
                    continue
                else:
                    raise ValueError(f"Source molecule {self.mname} {self.mnum} does not have {aname} atom to delete.")
                
            self.delete_atom(aname)
        
        for aname in embed_atoms:
            a = other[aname]
            self.add_atom(a)











