from typing import Union, List
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule, PdbResidue, NucleicAcidResidue, AminoacidResidue
from ...exceptions import InvalidPDB


    
class PDBCompounds:
    __slots__ = ("__comps", )
    
    def __init__(self):
        self.__comps = []
        
        
    def __getitem__(self, i: int):
        return self.__comps[i]
    
    def __len__(self):
        return len(self.__comps)
    
    def __iter__(self):
        return iter(self.__comps)
        
        
    def add(self, compound: Union[PdbMolecule, 
                                  NucleicAcidResidue, AminoacidResidue, 
                                  "NucleicAcidChain", "ProteinChain"]):   
        self.__comps.append(compound)


# CHAIN
    
class PDBChain(PDBCompounds):
    def __init__(self):
        super().__init__()
        
    def __repr__(self):
        return f"{self.__class__.__name__} with {len(self)} {self[0].__class__.__name__} residues at {hex(id(self))}"
                
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
        
        
class NucleicAcidChain(PDBChain):
    def __init__(self):
        super().__init__()
        
    def add(self, residue: NucleicAcidResidue):
        if not isinstance(residue, NucleicAcidResidue):
            last_residue_message = ""
            if len(self):
                last_residue_message = f"Last residue starts with atom - {residue[0].atomn} {residue[0].name}"
            raise InvalidPDB(f"Expected NucleicAcidResidue, got {type(residue)}. {last_residue_message}")
            
        super().add(residue)
        
        
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
                s.append(f"[{i:>3}] - {repr(c)}")
            
            else:
                if mol_name is not None and c[0].mol_name!=mol_name:
                    s.append(f"[{moli:>3}] - {nmols} {molt} with name {mol_name}")
                    nmols = 0
                    
                if nmols==0: 
                    moli = i
                    molt = c.__class__.__name__
                    mol_name = c[0].mol_name
                nmols+=1
        
        if nmols: s.append(f"[{moli:>3}] - {nmols} {molt} with name {mol_name}")
        return "\n    ".join(s)
        
        
class PDBModels:
    __slots__ = ("__models", "header")
    
    def __init__(self, models: List[PDB], header: str = ""):
        self.__models = models
        self.header = header
        
    def __getitem__(self, i: int):
        return self.__models[i]
    
    def __len__(self):
        return len(self.__models)