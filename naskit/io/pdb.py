from typing import Union, List, Dict, Optional
from pathlib import Path
from io import TextIOWrapper
from tempfile import _TemporaryFileWrapper
import urllib
import tempfile
import numpy as np

from ..containers.pdb.pdbAtom import PdbAtom
from ..containers.pdb.pdbMolecule import PdbMolecule
from ..containers.pdb.pdbResidue import PdbResidue, NucleicAcidResidue, AminoacidResidue
from ..containers.pdb.pdbContainer import PDB, PDBModels, NucleicAcidChain, ProteinChain
from ..exceptions import InvalidPDB



NA_NAMES = {"A", "U", "G", "C", "I", "T",
            "DA", "DT", "DG", "DC"}

AMINOACID_NAMES = {'ALA', 'CYS', 'ASP', 'GLU', 
                   'PHE', 'GLY', 'ILE', 'LYS', 
                   'LEU', 'MET', 'PRO', 'GLN', 
                   'ARG', 'SER', 'THR', 'VAL', 
                   'TRP', 'TYR', 'ASH', 'ASN', 
                   'HIS', 'HID', 'HIE', 'HIP',
                          'HSD', 'HSE', 'HSP'}


class pdbRead:
    def __init__(self, file: Union[str, Path, TextIOWrapper, _TemporaryFileWrapper]):
        if isinstance(file, (str, Path)):
            self._file = open(file)
        elif isinstance(file, (TextIOWrapper, _TemporaryFileWrapper)):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper. Got {type(file)}")
        
    def __enter__(self):
        return self
    
    def close(self):
        self._file.close()

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
        
    def read(self, 
             derive_element: bool = False, 
             element_derive_func = None
            ):
        
        lines = self._file.readlines()
        lines = list(map(lambda l: l.rstrip(), lines))
        
        # Header
        for i, l in enumerate(lines):
            if l.startswith("ATOM") or \
                l.startswith("HETATM") or \
                l.startswith("MODEL"):
                break
                
        header = "" if i==0 else "\n".join(lines[:i])
        lines = lines[i:]
        tokens = self.parse_atoms(lines, 
                                  derive_element=derive_element, 
                                  element_derive_func=element_derive_func)
        models = self.split_models(tokens)
        
        for i, model in enumerate(models):
            model = self.parse_mols(model)
            model = self.parse_chains(model)
            pdb = PDB()
            for component in model: pdb.add(component)
            models[i] = pdb
            
        return PDBModels(models, header)
    
        
    def parse_atoms(self, 
                    lines, 
                    derive_element: bool = False, 
                    element_derive_func = None
                   ):
        
        tokens = []
        
        for l in lines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                atom = PdbAtom.from_pdb_line(l, 
                                             derive_element=derive_element, 
                                             element_derive_func=element_derive_func)
                if atom.altloc not in (' ', 'A'): # anisotropic temperature factors
                    continue
                tokens.append(atom)
                
            elif l.startswith("TER") or l.startswith("MODEL") or l.startswith("ENDMDL"):
                tokens.append(l.strip())
                
        return tokens
    
    
    def split_models(self, tokens):
        models = []
        model_state = 0 # 0 - undefined, 1 - opened (after MODEL), 2 - closed (ENDMDL)
        model_components = []
        
        for a in tokens:
            if isinstance(a, PdbAtom) or a.startswith("TER"):
                model_components.append(a)
                
            elif a.startswith("MODEL"):
                if model_state==1:
                    raise InvalidPDB(f"PDB contains second MODEL line ({a}) without closing ENDMDL before.")
                    
                if len(model_components):
                    raise InvalidPDB(f"PDB contains atom or TER lines before opening ({a}) line.")
                    
                model_state = 1
                
            else: # ENDMDL
                if model_state!=1:
                    raise InvalidPDB(f"PDB contains closing ({a}) line without opening MODEL line before.")
                    
                if len(model_components)==0:
                    raise InvalidPDB(f"PDB contains empty model. Closed at line ({a}).")
                    
                models.append(model_components)
                model_components = []
                model_state = 2
                
        if len(model_components):
            if model_state!=0:
                raise InvalidPDB(f"The last PDB model is not closed by ENDMDL.")
            else: # signle model without MODEL and ENDMDL in pdb file
                models.append(model_components)
        
        return models
    
    
    def parse_mols(self, atom_tokens):
        mol_tokens = []
        mol_atoms = []
        reset_mol_idx = True
        
        for at in atom_tokens:
            if isinstance(at, PdbAtom):
                if reset_mol_idx:
                    cur_mol_idx = at.mnum
                    reset_mol_idx = False
                    
                if at.mnum!=cur_mol_idx:
                    if len(mol_atoms)>0:
                        mol_tokens.append(self.make_mol(mol_atoms))
                    mol_atoms = []
                    cur_mol_idx = at.mnum
                    
                mol_atoms.append(at)
                
            else:
                if len(mol_atoms):
                    mol_tokens.append(self.make_mol(mol_atoms))
                mol_atoms = []
                reset_mol_idx = True
                mol_tokens.append(at)
                
        if len(mol_atoms):
            mol_tokens.append(self.make_mol(mol_atoms))
                
        return mol_tokens
                
        
    def parse_chains(self, mol_tokens):
        chains = []
        chain_mols = []
        reset_chain = True
        
        for i, m in enumerate(mol_tokens):
            if isinstance(m, (NucleicAcidResidue, AminoacidResidue)):
                if reset_chain:
                    cur_chain = m.chain
                    reset_chain = False
                    
                if m.chain!=cur_chain:
                    chains.append(self.make_chain(chain_mols))
                    chain_mols = []
                    cur_chain = m.chain
                    
                chain_mols.append(m)
                
            elif isinstance(m, PdbMolecule):
                if len(chain_mols):
                    chains.append(self.make_chain(chain_mols))
                    chain_mols = []
                    reset_chain = True
                chains.append(m)
                
            else: # TER
                if i==0:
                    raise InvalidPDB(f"TER on first line.")
                elif isinstance(mol_tokens[i-1], str) and mol_tokens[i-1].startswith("TER"):
                    raise InvalidPDB(f"Two TERs in a row.")
                
                if len(chain_mols)!=0:
                    chains.append(self.make_chain(chain_mols))
                    chain_mols = []
                reset_chain = True
                
        if len(chain_mols):
            chains.append(self.make_chain(chain_mols))
            
        return chains
        
        
    def make_mol(self, atoms):
        if atoms[0].mname in NA_NAMES:
            m = NucleicAcidResidue()
        elif atoms[0].mname in AMINOACID_NAMES:
            m = AminoacidResidue()
        else:
            m = PdbMolecule()
            
        for atom in atoms:
            m.add_atom(atom)
        return m
    
    
    def make_chain(self, mols):
        if isinstance(mols[0], NucleicAcidResidue):
            chain = NucleicAcidChain()
        elif isinstance(mols[0], AminoacidResidue):
            chain = ProteinChain()
        else:
            raise ValueError(f"Expected NucleicAcidResidue or AminoacidResidue, "
                             f"got first molecule of type {type(mols[0])}.")
        
        for m in mols:
            chain.add(m)
        
        return chain
    
    
class pdbWrite:
    def __init__(self, file: Union[str, Path, TextIOWrapper, _TemporaryFileWrapper]):  
        if isinstance(file, (str, Path)):
            self._file = open(file, 'w')
        elif isinstance(file, (TextIOWrapper, _TemporaryFileWrapper)):
            self._file = file
        else:
            raise TypeError(f"Invalid file type. Accepted - string, Path, TextIOWrapper. Got {type(file)}")
        
    def __enter__(self):
        return self
    
    def close(self):
        self._file.close()

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


    def write(self, 
              data: Union[PdbMolecule, NucleicAcidResidue, AminoacidResidue,
                          NucleicAcidChain, ProteinChain, 
                          PDB, PDBModels],
              write_header: bool = False
             ):
        
        if not isinstance(data, (PDBModels, PDB, 
                               NucleicAcidChain, ProteinChain, 
                               PdbMolecule, NucleicAcidResidue, AminoacidResidue)):
            raise TypeError(f"Pdb writer can not write object of type {type(data)}. "
                             f"Expected - PdbMolecule, NucleicAcidResidue, AminoacidResidue, "
                             f"NucleicAcidChain, ProteinChain, PDB, PDBModels.")
        
        if isinstance(data, PDBModels) and write_header and len(data.header):
            self._file.write(data.header + "\n")
        
        self._file.write(str(data))


def request_pdb(pdb_id: str, 
                save_path: Optional[str] = None, 
                **kwargs
               ) -> PDB:
    
    req = urllib.request.Request(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    with urllib.request.urlopen(req) as response:
       txt = response.read().decode()

    if save_path is not None:
        with open(save_path, 'w') as f:
            f.write(txt)
        
    fp = tempfile.TemporaryFile('w+')
    fp.write(txt)
    fp.seek(0)

    with pdbRead(fp) as f:
        pdb = f.read(**kwargs)
    
    return pdb
    