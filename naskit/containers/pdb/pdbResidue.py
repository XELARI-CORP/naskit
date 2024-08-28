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


    def is_purine(self) -> bool:
        return self.name.lstrip("D") in "AG" and all([a in self for a in PURINE_CORE_ATOMS])

    def is_pyrimidine(self) -> bool:
        return self.name.lstrip("D") in "CUT" and all([a in self for a in PYRIMIDINE_CORE_ATOMS])

    def base_normal_vec(self):
        anames = ("N9", "C4", "C8") if self.is_purine() else ("N1", "C2", "C6")
        o, a, b = [self[aname] for aname in anames]
        
        v1 = a.coords-o.coords 
        v2 = b.coords-o.coords
        n = np.cross(v1, v2)
        n /= np.linalg.norm(n)
        return n

    
    def _add_atom_on_axis(self, dira1: str, dira2: str, aname: str, bond_len: float):
        ac = self[dira2].copy()
        
        bondv = self[dira2].coords - self[dira1].coords
        bondv = bondv/np.linalg.norm(bondv)
        
        ac.coords = self[dira2].coords + bond_len*bondv
        ac.element = aname[0]
        ac.name = aname
        ac.atomn = 1 + max([a.atomn for a in self])
        self.add_atom(ac)
        
    
    def fill3end(self):
        if "HO3'" not in self:
            self._add_atom_on_axis("C3'", "O3'", "HO3'", 0.96)

    
    def fill5end(self):
        if 'P' in self:
            self._add_atom_on_axis("O5'", "P", "OP3", 1.6097)
            self._add_atom_on_axis("P", "OP3", "HOP3", 0.96)
        elif "HO5'" not in self:
            self._add_atom_on_axis("C5'", "O5'", "HO5'", 0.96)
        
    
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
        
        new_mol = NT_TEMPLATE_MAP.get(base)
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

RIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "C2'", "O2'", "C1'", "H2'", "HO2'")
DEOXYRIBOSE_SOURSE_ATOMS = ("C4'", "O4'", "C3'", "C2'", "C1'",
                            # alias
                            ("H2'", "H2'1"), 
                            ("H2''", "H2'2")
                           )
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

PYRIMIDINE_CORE_ATOMS = ("N1", "C2", "N3", "C4", "C5", "C6")
PURINE_CORE_ATOMS = ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4")
PYRIMIDINE_TO_PURINE_CORE_CORESPONDANCE = (("N1", "N9"), ("C2", "C4"), ("C6", "C8"))

BASE_NAME_SOURCE_ATOMS_MAP = {
    "A": PURINE_CORE_ATOMS + ("N6", "H8", "H61", "H62", "H2"), 
    "G": PURINE_CORE_ATOMS + ("O6", "N2", "H8", "H1", "H21", "H22"), 
    "C": PYRIMIDINE_CORE_ATOMS + ("O2", "N4", "H41", "H42", "H5", "H6"), 
    "U": PYRIMIDINE_CORE_ATOMS + ("O2", "O4", "H3", "H5", "H6"), 
    "T": PYRIMIDINE_CORE_ATOMS + ("O2", "O4", "C7", "H3", "H71", "H72", "H73", "H6")
                             }

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
        
NT_TEMPLATE_MAP = {}
for nt_name, pdb_name in zip(["A", "G", "C", "U", "T"], 
                             ["adenine.pdb", "guanine.pdb", # Purine
                              "cytosine.pdb", "uracil.pdb", "thymine.pdb"] # Pyrimidine
                            ):
    
    with open(PACKAGE_PATH/"resources"/"pdb"/pdb_name) as f:
        NT = NucleicAcidResidue()
        for l in f:
            NT.add_atom(PdbAtom.from_pdb_line(l))

        NT_TEMPLATE_MAP[nt_name] = NT







