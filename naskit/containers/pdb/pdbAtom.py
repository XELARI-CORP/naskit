import numpy as np
from typing import Union
from .atom_name_mapping import ATOM_RENAME_MAP
from ...exceptions import InvalidPDB


    
class PdbAtom:
    __slots__ = ("is_hetatm", "anum", 
                 "aname", "altloc", 
                 "mname", "chain", "mnum", 
                 "occupancy", "temp", 
                 "segment", 
                 "element", "charge", 
                 "coords")
    
    def __init__(
                self,
                is_hetatm: bool, anum: int, 
                aname: str, altloc: str,
                mname: str, chain: str, mnum: int,
                x: float, y: float, z: float,
                occupancy: float, temp: float,
                segment: str,
                element: str, charge: int
                ):
        
        self.is_hetatm = is_hetatm
        self.anum = anum
        self.aname = aname
        self.altloc = altloc
        self.mname = mname
        self.chain = chain
        self.mnum = mnum
        self.occupancy = occupancy
        self.temp = temp
        self.segment = segment
        self.element = element
        self.charge = charge
        self.coords = np.array([x,y,z], dtype=np.float32)
    
    
    def __str__(self):
        
        atom_type = "HETATM" if self.is_hetatm else "ATOM"

        aname = self.aname
        if len(aname)==4:
            aname = aname
        elif len(aname)==1:
            aname = f" {aname}  "
        elif len(self.element)==1: # H C
            aname = f" {aname:<3}"    
        else: # two characters element (Cl, Fe ...)
            aname = aname.ljust(4, ' ')
            
        mname = f"{self.mname:>3}".ljust(4)
        mnum = f"{self.mnum:>4}".ljust(5)
        
        occupancy = f"{self.occupancy:>6.2f}" if isinstance(self.occupancy, float) else " "*6
        temp = f"{self.temp:>6.2f}" if isinstance(self.temp, float) else " "*6
        
        charge = ''
        if self.charge!=0:
            charge = ('-', '+')[int(self.charge>0)] + str(abs(self.charge))
            charge = charge.rstrip('1')
        
        return (f"{atom_type:<6}{self.anum:>5} "
                f"{aname}{self.altloc:>1}"
                f"{mname}{self.chain}{mnum}   "
                f"{self.coords[0]:>8.3f}{self.coords[1]:>8.3f}{self.coords[2]:>8.3f}"
                f"{occupancy}{temp}      "
                f"{self.segment:<4}{self.element:>2}{charge:<2}")

    
    @property
    def x(self):
        return self.coords[0]
    
    @property
    def y(self):
        return self.coords[1]
    
    @property
    def z(self):
        return self.coords[2]
    
    def copy(self):
        return self.__class__(is_hetatm=self.is_hetatm, anum=self.anum, 
                              aname=self.aname, altloc=self.altloc,
                              mname=self.mname, chain=self.chain, mnum=self.mnum,
                              x=self.x, y=self.y, z=self.z,
                              occupancy=self.occupancy, temp=self.temp,
                              segment=self.segment, element=self.element, charge=self.charge)
    
    def as_dict(self):
        return dict(is_hetatm=self.is_hetatm, anum=self.anum, 
                    aname=self.aname, altloc=self.altloc,
                    mname=self.mname, chain=self.chain, mnum=self.mnum,
                    coords=self.coords,
                    occupancy=self.occupancy, temp=self.temp,
                    segment=self.segment, element=self.element, charge=self.charge)

    def dist(self, a: Union["PdbAtom", "PdbMolecule"]):
        if len(a.coords.shape)==1:
            return np.linalg.norm(a.coords - self.coords)
        return np.linalg.norm((a.coords - self.coords), axis=1)
    
    @staticmethod
    def _default_element_derive_func(is_hetatm: bool, aname: str, mname: str, chain: str):
        return aname[0]
    
    @classmethod
    def from_pdb_line(cls,
                      line,
                      derive_element: bool = False,
                      element_derive_func = None
                     ):
        """
        http://www.wwpdb.org/documentation/file-format
        or
        https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf
        """
        if derive_element and element_derive_func is None:
            element_derive_func = PdbAtom._default_element_derive_func
            
        line = line.ljust(80, ' ')
        
        is_hetatm = line.startswith("HETATM")  # ATOM or HETATM
        anum = int(line[6:11].strip())         # Atom serial number
        aname = line[12:16].strip()            # Atom name
        altloc = line[16]                      # Alternate location indicator
        mname = line[17:21].strip()            # Residue/mol name. Must be [17:20], used extended range [17:21]
        chain = line[21]                       # Chain identifier
        mnum = int(line[22:27].strip())        # Residue sequence number
        # Ignore insertion code at 26
        x = float(line[30:38].strip())         # X
        y = float(line[38:46].strip())         # Y
        z = float(line[46:54].strip())         # Z
        
        occupancy = line[54:60].strip()        # Occupancy
        occupancy = float(occupancy) if occupancy else 1.
            
        temp = line[60:66].strip()             # Temperature factor
        temp = float(temp) if temp else 0.
        
        segment = line[72:76]                  # Segment identifier
        element = line[76:78].strip()          # Element symbol
        if element=="":
            if not derive_element:
                raise InvalidPDB(f"Atom {anum} {aname} has no element field.")
            element = element_derive_func(is_hetatm, aname, mname, chain)     
        
        charge = line[78:80].strip()           # Charge
        if charge=='': charge = "+0"
        if len(charge)==1: charge = "1"+charge
        sign, charge = sorted(charge)
        try:
            charge = int(charge)
            if sign=='+':
                pass
            elif sign=='-':
                charge *= -1
            else:
                raise InvalidPDB(f"Invalid atom charge sign '{sign}' in atom {anum} {aname}.")
        except:
            raise InvalidPDB(f"Invalid atom charge '{charge}' in atom {anum} {aname}.")

        return PdbAtom(is_hetatm, anum, 
                       aname, altloc, 
                       mname, chain, mnum, 
                       x, y, z, occupancy, temp, segment, element, charge)







