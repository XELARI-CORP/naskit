import numpy as np
from ...exceptions import InvalidPDB


    
class PdbAtom:
    __slots__ = ("is_hetatm", "atomn", "name", "altloc", 
                 "mol_name", "chain", "moln", 
                 "occupancy", "temp", 
                 "segment", 
                 "element", "charge", 
                 "coords")
    
    def __init__(
                self,
                is_hetatm: bool, atomn: int, name: str, altloc: str,
                mol_name: str, chain: str, moln: int,
                x: float, y: float, z: float,
                occupancy: float, temp: float,
                segment: str,
                element: str, charge: int
                ):
        
        self.is_hetatm = is_hetatm
        self.atomn = atomn
        self.name = name
        self.altloc = altloc
        self.mol_name = mol_name
        self.chain = chain
        self.moln = moln
        self.occupancy = occupancy
        self.temp = temp
        self.segment = segment
        self.element = element
        self.charge = charge
        self.coords = np.array([x,y,z], dtype=np.float32)
    
    
    def __str__(self):
        atom_type = "HETATM" if self.is_hetatm else "ATOM"
        
        if len(self.name)==4:
            name = self.name
        elif len(self.name)==1:
            name = f" {self.name}  "
        elif len(self.element)==1: # H C
            name = f" {self.name:<3}"    
        else: # two characters element (Cl, Fe ...)
            name = self.name.ljust(4, ' ')
            
        mol_name = f"{self.mol_name:>3}".ljust(4)
        moln = f"{self.moln:>4}".ljust(5)
        
        occupancy = f"{self.occupancy:>6.2f}" if isinstance(self.occupancy, float) else " "*6
        temp = f"{self.temp:>6.2f}" if isinstance(self.temp, float) else " "*6
        
        charge = ''
        if self.charge!=0:
            charge = ('-', '+')[int(self.charge>0)] + str(abs(self.charge))
            charge = charge.rstrip('1')
        
        return (f"{atom_type:<6}{self.atomn:>5} "
                f"{name}{self.altloc:>1}"
                f"{mol_name}{self.chain}{moln}   "
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
        return self.__class__(is_hetatm=self.is_hetatm, atomn=self.atomn, name=self.name, altloc=self.altloc,
                              mol_name=self.mol_name, chain=self.chain, moln=self.moln,
                              x=self.x, y=self.y, z=self.z,
                              occupancy=self.occupancy, temp=self.temp,
                              segment=self.segment, element=self.element, charge=self.charge)
    
    def as_dict(self):
        return dict(is_hetatm=self.is_hetatm, atomn=self.atomn, name=self.name, altloc=self.altloc,
                      mol_name=self.mol_name, chain=self.chain, moln=self.moln,
                      coords=self.coords,
                      occupancy=self.occupancy, temp=self.temp,
                      segment=self.segment, element=self.element, charge=self.charge)
    
    @staticmethod
    def _default_element_derive_func(is_hetatm: bool, name: str, mol_name: str, chain: str):
        return name[0]
    
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
        atomn = int(line[6:11].strip())        # Atom serial number
        name = line[12:16].strip()             # Atom name
        altloc = line[16]                      # Alternate location indicator
        mol_name = line[17:21].strip()         # Residue/mol name. Must be [17:20], used extended range [17:21]
        chain = line[21]                       # Chain identifier
        moln = int(line[22:27].strip())        # Residue sequence number
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
                raise InvalidPDB(f"Atom {atomn} {name} has no element field.")
            element = element_derive_func(is_hetatm, name, mol_name, chain)     
        
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
                raise InvalidPDB(f"Invalid atom charge sign '{sign}' in atom {atomn} {name}.")
        except:
            raise InvalidPDB(f"Invalid atom charge '{charge}' in atom {atomn} {name}.")
        
        return PdbAtom(is_hetatm, atomn, name, altloc, mol_name, chain, moln,
                        x, y, z, occupancy, temp, segment, element, charge)