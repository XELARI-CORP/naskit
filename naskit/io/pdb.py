from typing import Union, List, Dict
from pathlib import Path
from io import TextIOWrapper
from tempfile import _TemporaryFileWrapper
import numpy as np

from ..containers.pdb.pdbAtom import PdbAtom
from ..containers.pdb.pdbMolecule import PdbMolecule, PdbResidue, NucleicAcidResidue, AminoacidResidue
from ..containers.pdb.pdbContainer import PDB, PDBModels, NucleicAcidChain, ProteinChain
from ..exceptions import InvalidPDB



NA_NAMES = {"A", "U", "G", "C", "I", 
            "DA", "DT", "DG", "DC"}

AMINOACID_NAMES = {'ALA', 'CYS', 'ASP', 'GLU', 
                   'PHE', 'GLY', 'ILE', 'LYS', 
                   'LEU', 'MET', 'PRO', 'GLN', 
                   'ARG', 'SER', 'THR', 'VAL', 
                   'TRP', 'TYR', 'ASH',  'ASN', 
                   'HID', 'HIE', 'HIP', 'HIS'}


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
        
        
    def read(self):
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
        tokens = self.parse_atoms(lines)
        models = self.split_models(tokens)
        
        for i, model in enumerate(models):
            model = self.parse_mols(model)
            model = self.parse_chains(model)
            pdb = PDB()
            for component in model: pdb.add(component)
            models[i] = pdb
            
        return PDBModels(models, header)
    
        
    def parse_atoms(self, lines):
        tokens = []
        
        for l in lines:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                atom = PdbAtom.from_pdb_line(l)
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
                    cur_mol_idx = at.moln
                    reset_mol_idx = False
                    
                if at.moln!=cur_mol_idx:
                    if len(mol_atoms)>0:
                        mol_tokens.append(self.make_mol(mol_atoms))
                    mol_atoms = []
                    cur_mol_idx = at.moln
                    
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
        
        for m in mol_tokens:
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
                if len(chain_mols)==0:
                    raise InvalidPDB(f"PDB contains empty chain at ({m}).")
                    
                chains.append(self.make_chain(chain_mols))
                chain_mols = []
                reset_chain = True
                
        if len(chain_mols):
            chains.append(self.make_chain(chain_mols))
            
        return chains
        
        
    def make_mol(self, atoms):
        if atoms[0].mol_name in NA_NAMES:
            m = NucleicAcidResidue()
        elif atoms[0].mol_name in AMINOACID_NAMES:
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
              reenum_atoms: bool = False,
              reenum_mols: bool = False,
              rename_chains: bool = False,
              write_header: bool = False
             ):
        
        ## re enum
        
        if isinstance(data, PDBModels):
            if write_header and len(data.header):
                self._file.write(data.header + "\n")
            self.write_models(data)
            
        elif isinstance(data, PDB):
            self.write_pdb(data)
            
        elif isinstance(data, (NucleicAcidChain, ProteinChain)):
            self.write_chain(data)
            
        elif isinstance(data, (PdbMolecule, NucleicAcidResidue, AminoacidResidue)):
            self.write_molecule(data)
        
        else:
            raise TypeError(f"Pdb writer can not write object of type {type(data)}. "
                             f"Expected - PdbMolecule, NucleicAcidResidue, AminoacidResidue, "
                             f"NucleicAcidChain, ProteinChain, PDB, PDBModels.")
            
            
    def write_models(self, data):
        for i, m in enumerate(data):
            self._file.write(f"MODEL        {i+1}".ljust(80) + "\n")
            self.write_pdb(m)
            self._file.write(f"ENDMDL\n")
        
    def write_pdb(self, data):
        for i, c in enumerate(data):
            if isinstance(c, (NucleicAcidChain, ProteinChain)):
                self.write_chain(c)
            else:
                self.write_molecule(c)
        
    def write_chain(self, data):
        for m in data:
            self.write_molecule(m)
        self._file.write(f"TER\n")
        
    def write_molecule(self, data):
        for a in data:
            self._file.write(f"{str(a)}\n")