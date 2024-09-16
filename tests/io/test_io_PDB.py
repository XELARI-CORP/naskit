import pytest
import tempfile
from naskit.containers.pdb import PdbAtom, NucleicAcidResidue, AminoacidResidue, NucleicAcidChain, ProteinChain
from naskit import pdbRead
from naskit.exceptions import InvalidPDB



class TestReadPdbMakeMolecule:

    def test_no_element(self):
        pdb_text = """\
ATOM      1  O5'  DG     1      25.240  25.210  31.260  1.00  0.00            
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_mol_with_same_atom(self):
        pdb_text = """\
ATOM     29  N   ASN A   6      36.804  15.140  26.697  1.00 28.50           N  
ATOM     30  N   ASN A   6      36.804  15.140  26.697  1.00 28.50           N  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_dif_mol_name(self):
        pdb_text = """\
HETATM 2609  P   PO4 A 301      -4.401   0.932 108.199  1.00 45.91           P  
HETATM 2610  O1  L   A 301      -4.921  -0.287 107.368  1.00 56.99           O  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_dif_chain_in_mol(self):
        pdb_text = """\
HETATM 2609  P   PO4 A 301      -4.401   0.932 108.199  1.00 45.91           P  
HETATM 2610  O1  PO4 B 301      -4.921  -0.287 107.368  1.00 56.99           O 
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

                
class TestReadPdbSplitModels:
    
    def test_open_model(self):
        pdb_text = """\
MODEL
ATOM   2594  N   GLN B 172      -4.507 -21.576  70.408  1.00 19.99           N  
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

                
    def test_two_open(self):
        pdb_text = """\
MODEL
MODEL
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()
                
                
    def test_last_open_model(self):
        pdb_text = """\
MODEL
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
ENDMDL
MODEL
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()
                
                
    def test_close_model(self):
        pdb_text = """\
ATOM   2594  N   GLN B 172      -4.507 -21.576  70.408  1.00 19.99           N  
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
ENDMDL
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()
                
                
    def test_two_close(self):
        pdb_text = """\
MODEL
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
ENDMDL
ENDMDL
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_empty_model(self):
        pdb_text = """\
MODEL  
ENDMDL
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_atoms_outside_model(self):
        pdb_text = """\
MODEL
ATOM   2594  N   GLN B 172      -4.507 -21.576  70.408  1.00 19.99           N  
ENDMDL
ATOM   2595  CA  GLN B 172      -4.273 -22.613  69.383  1.00 22.77           C  
MODEL
ATOM   2599  CG  GLN B 172      -3.101 -24.913  69.186  1.00 34.07           C  
ENDMDL
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()

    
    def test_TER_outside_model(self):
        pdb_text = """\
MODEL
ATOM   2599  CG  GLN B 172      -3.101 -24.913  69.186  1.00 34.07           C  
ENDMDL
TER
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()


class TestReadPdbChain:
    
    def test_mix_na_amin(self):
        pdb_text = """\
ATOM    572  P     C A  19      67.015  33.891 -97.657  1.00  0.00           P  
ATOM    573  OP1   C A  19      67.343  32.450 -97.581  1.00  0.00           O  
ATOM    603  P     G A  20      71.325  33.829-102.603  1.00  0.00           P  
ATOM    604  OP1   G A  20      71.579  32.512-103.231  1.00  0.00           O  
ATOM    305  N   THR A  41      50.378  -1.061  42.314  1.00 40.42           N  
ATOM    306  CA  THR A  41      50.831  -1.724  43.529  1.00 37.70           C  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pytest.raises(InvalidPDB):
            with pdbRead(fp) as f:
                pdb = f.read()


    def test_add_residue_with_dif_chain_name(self):
        chain = NucleicAcidChain()
        
        m1 = NucleicAcidResidue()
        m1.add_atom(PdbAtom.from_pdb_line("ATOM    665  P     U A  22      78.912  34.354-105.928  1.00  0.00           P  "))
        chain.add(m1)
        
        m2 = NucleicAcidResidue()
        m2.add_atom(PdbAtom.from_pdb_line("ATOM    666  OP1   U B  23      79.027  35.512-106.840  1.00  0.00           O  "))
        
        with pytest.raises(InvalidPDB):
            chain.add(m2)
            
            
    def test_add_residue_with_same_mol_num(self):
        chain = NucleicAcidChain()
        
        m1 = NucleicAcidResidue()
        m1.add_atom(PdbAtom.from_pdb_line("ATOM    665  P     U A  22      78.912  34.354-105.928  1.00  0.00           P  "))
        chain.add(m1)
        
        m2 = NucleicAcidResidue()
        m2.add_atom(PdbAtom.from_pdb_line("ATOM    666  OP1   U A  22      79.027  35.512-106.840  1.00  0.00           O  "))
        
        with pytest.raises(InvalidPDB):
            chain.add(m2)   