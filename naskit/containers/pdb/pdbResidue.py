from typing import Union, List
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
        
        
# Molecule fragment embedding

RIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "H4'", "H3'", "H2'", "HO2'", "H1'")
DEOXYRIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "H4'", "H3'", "H2'", "H2''", "H1'")
RIBOSE_DEOXYRIBOSE_ALIGN_CORRESPONDENCE_ATOMS = (("C1'", "C1'"), ("C4'", "C4'"), ("O4'", "O4'"))

deoxyribose_lines = """\
ATOM      1  C4'  DG A   1      -0.209   0.437   1.353  1.00  0.00           C  
ATOM      2  O4'  DG A   1      -0.680   1.478   0.486  1.00  0.00           O  
ATOM      3  C3'  DG A   1       0.422  -0.647   0.457  1.00  0.00           C  
ATOM      4  O3'  DG A   1      -0.460  -1.731   0.221  1.00  0.00           O  
ATOM      5  C2'  DG A   1       0.604   0.079  -0.864  1.00  0.00           C  
ATOM      6  C1'  DG A   1      -0.638   0.976  -0.844  1.00  0.00           C  
ATOM      7  H4'  DG A   1      -1.039  -0.002   1.906  1.00  0.00           H  
ATOM      8  H3'  DG A   1       1.383  -0.994   0.845  1.00  0.00           H  
ATOM      9  H2'  DG A   1       1.525   0.662  -0.828  1.00  0.00           H  
ATOM     10 H2''  DG A   1       0.606  -0.608  -1.713  1.00  0.00           H  
ATOM     11  H1'  DG A   1      -1.516   0.347  -1.019  1.00  0.00           H  \
"""

DEOXYRIBOSE_CORE = NucleicAcidResidue()
for l in deoxyribose_lines.split("\n"):
    DEOXYRIBOSE_CORE.add_atom(PdbAtom.from_pdb_line(l))
    
    
ribose_lines = """\
ATOM      1  C4'   G A   1       0.780  -0.535  -1.189  1.00  0.00           C  
ATOM      2  O4'   G A   1       1.528   0.615  -0.740  1.00  0.00           O  
ATOM      3  C3'   G A   1      -0.659  -0.242  -0.834  1.00  0.00           C  
ATOM      4  O3'   G A   1      -1.434  -1.436  -0.666  1.00  0.00           O  
ATOM      5  C2'   G A   1      -0.535   0.533   0.450  1.00  0.00           C  
ATOM      6  O2'   G A   1      -0.586  -0.340   1.585  1.00  0.00           O  
ATOM      7  C1'   G A   1       0.818   1.198   0.380  1.00  0.00           C  
ATOM      8  H4'   G A   1       1.111  -1.402  -0.621  1.00  0.00           H  
ATOM      9  H3'   G A   1      -1.101   0.397  -1.603  1.00  0.00           H  
ATOM     10  H2'   G A   1      -1.307   1.293   0.520  1.00  0.00           H  
ATOM     11 HO2'   G A   1       0.022  -1.063   1.416  1.00  0.00           H  
ATOM     12  H1'   G A   1       1.361   0.985   1.302  1.00  0.00           H  \
"""

RIBOSE_CORE = NucleicAcidResidue()
for l in ribose_lines.split("\n"):
    RIBOSE_CORE.add_atom(PdbAtom.from_pdb_line(l))