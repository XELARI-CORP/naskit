from typing import Union, List, Tuple, Iterable
import numpy as np
from .pdbAtom import PdbAtom
from .pdbMolecule import PdbMolecule
from ...exceptions import InvalidPDB
from ...utils.package_resources import get_package_path



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
    
    
    def is_protonated(self):
        return any([a.element=='H' for a in self.atoms()])
        
    
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
            self._embed_fragment_with_hydrogen_check(RIBOSE_CORE,
                                                     source_atoms=DEOXYRIBOSE_SOURSE_ATOMS,
                                                     embed_atoms=RIBOSE_SOURSE_ATOMS,
                                                     correspondence=RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS)
            
        elif sugar=='deoxyribose' or sugar=='dna':
            self._embed_fragment_with_hydrogen_check(DEOXYRIBOSE_CORE, 
                                                     source_atoms=RIBOSE_SOURSE_ATOMS, 
                                                     embed_atoms=DEOXYRIBOSE_SOURSE_ATOMS, 
                                                     correspondence=RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS)
            
        else:
            raise ValueError(f"Sugar must be 'ribose', 'rna' or 'deoxyribose', 'dna'. got {sugar}")
            
    def change_nucleobase(self, base: str):
        if base==self.name.lstrip("D"):
            return
        
        new_mol = BASE_NAME_CONTAINER_MAP.get(base)
        if new_mol is None:
            raise ValueError(f"Expected base name A, G, C, U, T, got {base}.")
            
        source_atoms = BASE_NAME_SOURCE_ATOMS_MAP.get(self.name[-1])
        embed_atoms = BASE_NAME_SOURCE_ATOMS_MAP.get(base)
        
        self_type = NAME_BASE_TYPE_MAP[self.name[-1]]
        other_type = NAME_BASE_TYPE_MAP[base]
        correspondence = BASE_CORESPONDANCE_MAP[f"{self_type}-{other_type}"]
        source_origin_atom = BASE_ORIGIN_ATOM_MAP[self.name[-1]]
        embed_origin_atom = BASE_ORIGIN_ATOM_MAP[base]
        
        self._embed_fragment_with_hydrogen_check(new_mol,
                                                 source_atoms=source_atoms, 
                                                 embed_atoms=embed_atoms, 
                                                 correspondence=correspondence,
                                                 source_origin_atom=source_origin_atom,
                                                 embed_origin_atom=embed_origin_atom
                                                )
        
    def _embed_fragment_with_hydrogen_check(self, 
                                            other: Union["PdbMolecule", "AminoacidResidue", "NucleicAcidResidue"],
                                            source_atoms: Iterable[str],
                                            embed_atoms: Iterable[str],
                                            correspondence: Iterable[Tuple[str, str]],
                                            source_origin_atom: str = None,
                                            embed_origin_atom: str = None
                                           ):
        
        was_protonated = self.is_protonated()
        self.embed_molecule_fragment(other, source_atoms, embed_atoms, correspondence, source_origin_atom, embed_origin_atom)
        if not was_protonated:
            h_atoms = [a.name for a in self.atoms() if a.element=='H']
            for h_atom in h_atoms:
                self.delete_atom(h_atom)
        
        
# TEMPLATES
PACKAGE_PATH = get_package_path()

## SUGARE

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
        
## BASE

BASE_NAME_SOURCE_ATOMS_MAP = {
    "A":("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "H8", "H61", "H62", "H2"), 
    "G":("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "H8", "H1", "H21", "H22"), 
    "C":("N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "H41", "H42", "H5", "H6"), 
    "U":("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "H3", "H5", "H6"), 
    "T":("N1", "C2", "O2", "N3", "C4", "O4", "C5", "C7", "C6", "H3", "H71", "H72", "H73", "H6")
                             }

PYRIMIDINE_CORE_ATOMS = ("N1", "C2", "N3", "C4", "C5", "C6")
PURINE_CORE_ATOMS = ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4")
PYRIMIDINE_TO_PURINE_CORE_CORESPONDANCE = (("N1", "N9"), ("C2", "C4"), ("C6", "C8"))

NAME_BASE_TYPE_MAP = {
    "C":"Pyrimidine", "U":"Pyrimidine", "T":"Pyrimidine", 
    "A":"Purine", "G":"Purine"
                     }

BASE_CORESPONDANCE_MAP = {
    "Pyrimidine-Pyrimidine": tuple([(a, a) for a in PYRIMIDINE_CORE_ATOMS]),
    "Purine-Purine": tuple([(a, a) for a in PURINE_CORE_ATOMS]),
    "Pyrimidine-Purine": PYRIMIDINE_TO_PURINE_CORE_CORESPONDANCE,
    "Purine-Pyrimidine": tuple([(b, a) for a, b in PYRIMIDINE_TO_PURINE_CORE_CORESPONDANCE])
                         }

BASE_ORIGIN_ATOM_MAP = {"C":"N1", "U":"N1", "T":"N1", "A":"N9", "G":"N9"}

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
        
BASE_NAME_CONTAINER_MAP = {"A":ADENINE_NB, "G":GUANINE_NB, "C":CYTOSINE_NB, "U":URACIL_NB, "T":THYMINE_NB}