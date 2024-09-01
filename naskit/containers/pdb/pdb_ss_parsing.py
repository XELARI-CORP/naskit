import numpy as np

from ...parse_na import NA
from ..nucleic_acid import NucleicAcid
from .pdbResidue import NucleicAcidResidue



DONOR_ACCEPTOR_GROUPS = {
    'A':{
        'donors':(('H61', 'N6'), ('H62', 'N6')),
        'acceptors':('N1', 'N3', 'N7')
    },
    'U':{
        'donors':(('H3', 'N3'),),
        'acceptors':('O2', 'O4')
    },
    'G':{
        'donors':(('H1', 'N1'), ('H21', 'N2'), ('H22', 'N2')),
        'acceptors':('O6', 'N3', 'N7')
    },
    'C':{
        'donors':(('H41', 'N4'), ('H42', 'N4')),
        'acceptors':('O2', 'N3')
    },
    'T':{
        'donors':(('H3', 'N3'),),
        'acceptors':('O2', 'O4')
    },
}

APPROXIMATE_H_BOND_DIST = 0.96
EPSILON = {'N':0.71128, 'O':0.87864, 'H':0.06569}
SIGMA = {'N':0.325, 'O':0.295992, 'H':0.106908}

COULOMB_CONST = 138.935458 / 78
RNA_CHARGE = {
    "A":{'N1':-0.76150, 'N3':-0.69970, 'N7':-0.60730, 
         'N6':-0.90190, 'H61': 0.4115, 'H62': 0.4115},
    "G":{'O6':-0.55970, 'N3':-0.63230, 'N7':-0.57090, 
         'N1':-0.47870, 'H1': 0.3424, 'N2':-0.96720, 'H21': 0.4364, 'H22': 0.4364},
    "C":{'O2':-0.62520, 'N3':-0.75840, 
         'N4':-0.95300, 'H41': 0.4234, 'H42': 0.4234},
    "U":{'O2':-0.54770, 'O4':-0.57610, 
         'N3':-0.35490, 'H3': 0.3154}
}

DNA_CHARGE = {
    "A":{'N1':-0.76240, 'N3':-0.74170, 'N7':-0.61750, 
         'N6':-0.91230, 'H61': 0.4167, 'H62': 0.4167},
    "G":{'O6':-0.56990, 'N3':-0.66360, 'N7':-0.57250, 
         'N1':-0.50530, 'H1': 0.352, 'N2':-0.92300, 'H21': 0.4235, 'H22': 0.4235},
    "C":{'O2':-0.65480, 'N3':-0.77480, 
         'N4':-0.97730, 'H41': 0.4314, 'H42': 0.4314},
    "T":{'O2':-0.58810, 'O4':-0.55630, 
         'N3':-0.43400, 'H3': 0.342}
}

MIN_H_ENERGY_THRESHOLD = -1.0 #-1.0493852218717947

class SSParsing:
    
    def to_na(self, approximate_hs: bool = False):
        """
        10.1093/nar/gkv716
        https://academic.oup.com/nar/article/43/21/e142/2468098
        """
        
        energy_matrix = self.get_ss_energy_matrix(approximate_hs)
        adj = self.parse_ss_adjacency(energy_matrix, threshold=MIN_H_ENERGY_THRESHOLD)
        na = NucleicAcid.from_adjacency(adj, seq=self.seq)
        return na

    
    def parse_ss_adjacency(self, M, threshold):
        N = M.shape[-1]
        out = np.zeros(M.shape, dtype=np.int32)
        
        while True:
            mi = np.argmin(M)
            r, c = mi//N, mi%N
            
            v = M[r, c]
            if v>threshold:
                break    
            
            out[r,c] = 1
            M[r] = np.inf
            M[c] = np.inf
            M[:, r] = np.inf
            M[:, c] = np.inf
    
        out = (out + out.T)
        return out

    def can_form_pair(self, 
                      nt1: NucleicAcidResidue, 
                      nt2: NucleicAcidResidue,
                      verbose: bool = False
                     ) -> bool:
        # max distance between the two origins
        origin1 = nt1["N9"] if nt1.is_purine() else nt1["N1"]
        origin2 = nt2["N9"] if nt2.is_purine() else nt2["N1"]
        origin_dist = origin1.dist(origin2)
        if verbose: print(f"{origin_dist:.4f} - distance between the two origins")
        if origin_dist > 15:
            if verbose: print("FAILED")
            return False

        # min distance between direction atoms
        a1dir, a2dir = (nt1["N9"], nt1["N1"]) if nt1.is_purine() else (nt1["C6"], nt1["N3"])
        b1dir, b2dir = (nt2["N9"], nt2["N1"]) if nt2.is_purine() else (nt2["C6"], nt2["N3"])
        dir_dist = a1dir.dist(b1dir)
        min_dir_dist = a1dir.dist(a2dir) + b1dir.dist(b2dir)
        if verbose: print(f"{dir_dist:.4f} - direction distance (molecule lengths - {min_dir_dist:.4f})")
        if dir_dist < min_dir_dist:
            if verbose: print("FAILED")
            return False
        
        # angle between the base normal vectors
        norm1 = nt1.base_normal_vec()
        norm2 = nt2.base_normal_vec()
        normals_angle = np.arccos(abs(np.dot(norm1, norm2))) * 180 / np.pi
        if verbose: print(f"{normals_angle:.4f} - angle between the base normal vectors")
        if normals_angle > 65:
            if verbose: print("FAILED")
            return False
        
        # vertical separation between the base planes
        # shifted_origin_coords = origin2.coords - origin1.coords
        # sep_dist = abs(np.dot(shifted_origin_coords, norm1))
        # if verbose: print(f"{sep_dist:.4f} - vertical separation between the base planes")
        # if sep_dist > 2.5:
        #     if verbose: print("FAILED")
        #     return False
        
        # absence of stacking between the two bases
        
        # presence of at least one hydrogen bond


        return True

    
    def _add_h_bonds(self, 
                     bonds: dict, 
                     donor: NucleicAcidResidue, 
                     acceptor: NucleicAcidResidue, 
                     approximate_hs: bool
                    ): 
        
        datoms = DONOR_ACCEPTOR_GROUPS[donor.name.lstrip('D')]["donors"]
        aatoms = DONOR_ACCEPTOR_GROUPS[acceptor.name.lstrip('D')]["acceptors"]
        
        for hdname, dname in datoms:
            need_approximate = False
            donor_atom_name = hdname
            
            # if hdname not in donor:
            #     if not approximate_hs:
            #         raise ValueError(f"Residue {donor.name} {donor.moln} does not contain hydrogen in '{dname}' donor group.")

            #     need_approximate = True
            #     donor_atom_name = dname

            donor_atom = donor[donor_atom_name]
            
            for acceptor_atom_name in aatoms:
                acceptor_atom = acceptor[acceptor_atom_name]
                
                dist = donor_atom.dist(acceptor_atom)
                # if need_approximate:
                #     dist -= APPROXIMATE_H_BOND_DIST / 2

                eps = np.sqrt(EPSILON[donor_atom.element] * EPSILON[acceptor_atom.element])
                sigma = (SIGMA[donor_atom.element] + SIGMA[acceptor_atom.element]) / 2
                sigma6 = sigma**6
                sigma12 = sigma6**2
                dist6 = dist**6
                dist12 = dist6**2
                e_lj = 4*eps*(sigma12/dist12 + sigma6/dist6)

                if donor.natype=="rna":
                    qd = RNA_CHARGE[donor.name][donor_atom_name]
                else:
                    qd = DNA_CHARGE[donor.name][donor_atom_name]

                if acceptor.natype=="rna":
                    qa = RNA_CHARGE[acceptor.name][acceptor_atom_name]
                else:
                    qa = DNA_CHARGE[acceptor.name][acceptor_atom_name]
                
                e_c = COULOMB_CONST * qa * qd / dist
                
                bonds[(donor_atom_name, acceptor_atom_name)] = (dist, e_lj+e_c, e_lj, e_c)

    
    def calculate_h_bonds(self, 
                          nt1: NucleicAcidResidue, 
                          nt2: NucleicAcidResidue, 
                          approximate_hs: bool = False
                         ) -> dict:
        
        bonds = {} # {('donor atom', 'acceptor atom'): (dist, energy), ...}
        self._add_h_bonds(bonds, nt1, nt2, approximate_hs)
        self._add_h_bonds(bonds, nt2, nt1, approximate_hs)
        
        return bonds
        
        
    def get_ss_energy_matrix(self, approximate_hs: bool = False):
        E = np.full((len(self), len(self)), np.inf, dtype=np.float32)

        for i in range(len(self)):
            nt1 = self[i]
            for j in range(i+3, len(self)):
                nt2 = self[j]

                if not self.can_form_pair(nt1, nt2):
                    continue

                h_bonds = self.calculate_h_bonds(nt1, nt2, approximate_hs) # {('donor atom', 'acceptor atom'): (dist, energy), ...}
                bond_energy = sum([b[1] for b in h_bonds.values()])
                
                E[i, j] = bond_energy
                E[j, i] = bond_energy
                
        return E












