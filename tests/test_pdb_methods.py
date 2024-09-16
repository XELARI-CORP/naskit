import pytest
import tempfile
import numpy as np
from naskit import pdbRead
from naskit.exceptions import InvalidPDB



class TestPdbContainersEditing:

    def test_coords_assign(self):
        pdb_text = """\
HETATM 2622  C09 JFJ B 201     -13.824 -13.954  90.882  1.00 31.81           C  
HETATM 2623 CL10 JFJ B 201     -16.317 -12.897  91.038  0.90 31.26          CL  
ATOM      1  O5'   G A   1      59.712  40.180-111.625  1.00  0.00           O  
ATOM      2  C5'   G A   1      61.014  39.722-111.985  1.00  0.00           C  
ATOM     33  P     G A   2      58.321  36.610-112.262  1.00  0.00           P  
ATOM     34  OP1   G A   2      58.167  35.355-113.032  1.00  0.00           O  
ATOM     35  OP2   G A   2      58.072  37.917-112.908  1.00  0.00           O  
TER
HETATM 2609  P   PO4 A 301      -4.401   0.932 108.199  1.00 45.91           P  
        """
        fp = tempfile.TemporaryFile('w+')
        fp.write(pdb_text)
        fp.seek(0)
        
        with pdbRead(fp) as f:
            pdb = f.read()[0]
            
        c0 = np.arange(3*pdb.natoms).reshape(-1, 3)
        pdb.coords = c0
        assert pdb[2][0].coords[2] == 3*pdb.natoms - 1