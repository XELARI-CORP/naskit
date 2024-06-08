from typing import Union, List
import importlib.resources
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule
from ...exceptions import InvalidPDB



class PdbResidue(PdbMolecule):
    def __init__(self):
        super().__init__()
        
        
class AminoacidResidue(PdbResidue):
    def __init__(self):
        super().__init__()
    
        
class NucleicAcidResidue(PdbResidue):
    def __init__(self):
        super().__init__()
        
    
    @property
    def natype(self):
        return "dna" if "D" in self.name else "rna"
    
    
    def to_rna(self):
        if self.natype=='rna':
            return
        
        self.change_sugar('ribose')
        if self.name=="DT":
            self.change_nucleobase("U")
            self.name = "U"
        else:
            self.name = self.name[1] # DC -> C
        
    def to_dna(self):
        if self.natype=='dna':
            return
        
        self.change_sugar('deoxyribose')
        if self.name=="U":
            self.change_nucleobase("T")
            self.name = "DT"
        else:
            self.name = "D" + self.name
            
            
    def change_sugar(self, sugar: str):
        if sugar=='ribose' or sugar=='rna':
            self.embed_molecule_fragment(RIBOSE_CORE,
                                         source_atoms=DEOXYRIBOSE_SOURSE_ATOMS,
                                         embed_atoms=RIBOSE_SOURSE_ATOMS,
                                         correspondence=RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS)
            
        elif sugar=='deoxyribose' or sugar=='dna':
            self.embed_molecule_fragment(DEOXYRIBOSE_CORE, 
                                         source_atoms=RIBOSE_SOURSE_ATOMS, 
                                         embed_atoms=DEOXYRIBOSE_SOURSE_ATOMS, 
                                         correspondence=RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS)
            
        else:
            raise ValueError(f"Sugar must be 'ribose', 'rna' or 'deoxyribose', 'dna'. got {sugar}")
            
    def change_nucleobase(self, base: str):
        print(f'Change from {self.name} to {base}')
        
        
# TEMPLATES

PACKAGE_PATH = importlib.resources.files('naskit')

RIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "H4'", "H3'", "H2'", "HO2'", "H1'")
DEOXYRIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "H4'", "H3'", "H2'", "H2''", "H1'")
RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS = (("C1'", "C1'"), ("C4'", "C4'"), ("O4'", "O4'"))

DEOXYRIBOSE_CORE = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"deoxyribose.pdb") as f:
    for l in f:
        DEOXYRIBOSE_CORE.add_atom(PdbAtom.from_pdb_line(l[:-1]))
    
RIBOSE_CORE = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"ribose.pdb") as f:
    for l in f:
        RIBOSE_CORE.add_atom(PdbAtom.from_pdb_line(l))

# Pyrimidine

THYMINE_NB = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"thymine.pdb") as f:
    for l in f:
        THYMINE_NB.add_atom(PdbAtom.from_pdb_line(l))

URACIL_NB = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"uracil.pdb") as f:
    for l in f:
        URACIL_NB.add_atom(PdbAtom.from_pdb_line(l))
        
CYTOSINE_NB = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"cytosine.pdb") as f:
    for l in f:
        CYTOSINE_NB.add_atom(PdbAtom.from_pdb_line(l))
        
# Purine

ADENINE_NB = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"adenine.pdb") as f:
    for l in f:
        ADENINE_NB.add_atom(PdbAtom.from_pdb_line(l))
        
GUANINE_NB = NucleicAcidResidue()
with open(PACKAGE_PATH/"resources"/"pdb"/"guanine.pdb") as f:
    for l in f:
        GUANINE_NB.add_atom(PdbAtom.from_pdb_line(l))