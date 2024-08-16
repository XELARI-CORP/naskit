from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule
from .pdbResidue import PdbResidue, NucleicAcidResidue, AminoacidResidue
from .pdbDraw import PDBDraw
from .pdb_ss_parsing import SSParsing
from ...exceptions import InvalidPDB


    
class PDBCompounds(PDBDraw):
    __slots__ = ("__comps", )
    
    def __init__(self):
        self.__comps = []
        
        
    def __getitem__(self, i: Union[int, slice]):
        if isinstance(i, slice):
            slice_comp = self.__class__()
            for c in self.__comps[i]:
                PDBCompounds.add(slice_comp, c)

            return slice_comp
            
        return self.__comps[i]
    
    def __len__(self):
        return len(self.__comps)
    
    def __iter__(self):
        return iter(self.__comps)
    
    def __str__(self):
        return "\n".join([str(c) for c in self])
        
        
    def copy(self):
        copied_comp = self.__class__()
        for c in self.__comps:
            PDBCompounds.add(copied_comp, c.copy())
        
        return copied_comp

    def dist(self, a: Union["PdbAtom", "PdbMolecule", "PDBCompounds"]):
        if len(a.coords.shape)==1:
            return np.linalg.norm((self.coords - a.coords), axis=1)
        return np.linalg.norm((self.coords[:, np.newaxis, :] - a.coords), axis=2)
    
    def add(self, compound: Union[PdbMolecule, 
                                  NucleicAcidResidue, AminoacidResidue, 
                                  "NucleicAcidChain", "ProteinChain"]):   
        self.__comps.append(compound)    
        
    def renum_atoms(self, initn: int = 1):
        offset = 0
        for c in self.__comps:
            c.renum_atoms(initn + offset)
            offset += c.natoms
            
    def renum_mols(self, initn: int = 1):
        offset = 0
        for c in self.__comps:
            if isinstance(c, (PdbMolecule, NucleicAcidResidue, AminoacidResidue)):
                c.moln = (initn + offset)
                offset += 1
            else: # chain
                c.renum_mols(initn + offset)
                offset += len(c)
                
    def rename_chains(self, chain_name: str = "A"):
        ch = ord(chain_name)
        for i, c in enumerate(self.__comps):
            if isinstance(c, (PdbMolecule, NucleicAcidResidue, AminoacidResidue)):
                c.chain = chr(ch)
            else:
                c.rename_chains(chr(ch))
                ch = ch+1 if ch<90 else 65
                
        
    def atoms(self):
        for c in self.__comps:
            for a in c.atoms():
                yield a
            
            
    @property
    def natoms(self):
        return sum([c.natoms for c in self.__comps])
    
    @property
    def coords(self):
        return np.concatenate([c.coords for c in self.__comps], axis=0)
    
    @coords.setter
    def coords(self, coords: np.ndarray):
        if coords.shape[0]!=self.natoms or coords.shape[1]!=3:
            raise ValueError(f"Coords matrix must have shape: ({self.natoms}, 3), got {coords.shape}")
            
        offset = 0
        for c in self.__comps:
            c.coords = coords[offset:(offset + c.natoms)]
            offset += c.natoms


# CHAIN
    
class PDBChain(PDBCompounds):
    def __init__(self):
        super().__init__()
        
    def __repr__(self):
        return f"{self.__class__.__name__} with {len(self)} {self[0].__class__.__name__} at {hex(id(self))}"
    
    def __str__(self):
        return "\n".join([str(c) for c in self]) + "\nTER"
                
    def add(self, residue: Union[NucleicAcidResidue, AminoacidResidue]):
        if len(self):
            res = self[0]
            if residue.chain!=res.chain:
                raise InvalidPDB(f"All residues of a chain must have the same chain name ({res.chain}), "
                                 f"got {residue.chain} in residue with first atom - "
                                 f"{residue[0].atomn} {residue[0].name}.")
                
            if residue.moln==res.moln:
                raise InvalidPDB(f"Residue with molecule number {residue.moln} "
                                 f"already exists in chain. Tried to add residue with first atom - "
                                 f"{residue[0].atomn} {residue[0].name}.")
            
        super().add(residue)
        
        
class NucleicAcidChain(PDBChain, SSParsing):
    def __init__(self):
        super().__init__()
        
    def add(self, residue: NucleicAcidResidue):
        if not isinstance(residue, NucleicAcidResidue):
            last_residue_message = ""
            if len(self):
                last_residue_message = f"Last residue starts with atom - {residue[0].atomn} {residue[0].name}"
            raise InvalidPDB(f"Expected NucleicAcidResidue, got {type(residue)}. {last_residue_message}")
            
        super().add(residue)
        
        
    def to_rna(self):
        for c in self: c.to_rna()
        
    def to_dna(self):
        for c in self: c.to_dna()
        
    @property
    def natype(self):
        if all([c.natype=='rna' for c in self]):
            return 'rna'
        elif all([c.natype=='dna' for c in self]):
            return 'dna'
        else:
            return None

    @property
    def seq(self):
        return "".join([res.name.lstrip('D') for res in self])
            
        
class ProteinChain(PDBChain):
    def __init__(self):
        super().__init__()
        
    def add(self, residue: AminoacidResidue):
        if not isinstance(residue, AminoacidResidue):
            last_residue_message = ""
            if len(self):
                last_residue_message = f"Last residue starts with atom - {residue[0].atomn} {residue[0].name}"
            raise InvalidPDB(f"Expected AminoacidResidue, got {type(residue)}. {last_residue_message}")
            
        super().add(residue)
        
        
# PDB container
        
class PDB(PDBCompounds):
    def __init__(self):
        super().__init__()
        
    def __repr__(self):
        s = [f"PDB at {hex(id(self))}"]
        nmols = 0
        moli = None
        molt = None
        mol_name = None
        for i, c in enumerate(self):
            if isinstance(c, (NucleicAcidChain, ProteinChain)):
                if nmols:
                    s.append(f"[{moli:>3}] - {nmols} {molt} with name {mol_name}")
                    nmols = 0
                    mol_name = None
                s.append(f"[{i:>3}] - {repr(c)}")
            
            else:
                if mol_name is not None and c[0].mol_name!=mol_name:
                    s.append(f"[{moli:>3}] - {nmols} {molt} with name {mol_name}")
                    nmols = 0
                    
                if nmols==0: 
                    moli = i
                    molt = c.__class__.__name__
                    mol_name = c.name
                nmols+=1
        
        if nmols: s.append(f"[{moli:>3}] - {nmols} {molt} with name {mol_name}")
        return "\n    ".join(s)
        
        
class PDBModels:
    __slots__ = ("__models", "header")
    
    def __init__(self, models: List[PDB], header: List[str] = []):
        self.__models = models
        self.header = header
        
    def __getitem__(self, i: int):
        return self.__models[i]
    
    def __len__(self):
        return len(self.__models)
    
    def __iter__(self):
        return iter(self.__models)
    
    def __str__(self):
        s = []
        for i, m in enumerate(self.__models):
            s.append(f"MODEL        {i+1}".ljust(80))
            s.append(str(m))
            s.append(f"ENDMDL")
        
        return "\n".join(s)