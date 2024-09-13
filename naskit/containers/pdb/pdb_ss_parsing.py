from dataclasses import dataclass
from typing import List
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

APPROXIMATE_H_BOND_DIST = 1.0
H_BOND_DISTANCE_CUTOFF = 4.5
MIN_H_ENERGY_THRESHOLD = -0.125

ORIGIN_DISTANCE_THRESHOLD = 15
NORMALS_ANGLE_THRESHOLD = 75
COPLANAR_ANGLE_THRESHOLD = 25


@dataclass
class HBond:
    donor_res_name: str
    acceptor_res_name: str
    donor_name: str
    acceptor_name: str
    dist: float
    approximated: bool
    bond_e: float
    LJ_e: float
    C_e: float
    

class SSParsing:
    
    def to_na(self, 
              approximate_hs: bool = False,
              verbose: bool = False
             ):
        """
        Geometry analysis inspired by DSSR - 10.1093/nar/gkv716
        Atomic charges and Lennard-Jones parameters are from amber99bsc1 force field
        """
        
        energy_matrix = self.get_ss_energy_matrix(approximate_hs, verbose)
        adj = self.parse_ss_adjacency(energy_matrix, threshold= 2*MIN_H_ENERGY_THRESHOLD)
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
        if origin_dist > ORIGIN_DISTANCE_THRESHOLD:
            if verbose: print(f"FAILED - must be < {ORIGIN_DISTANCE_THRESHOLD}")
            return False

        # min distance between direction atoms
        a1dir, a2dir = (nt1["N9"], nt1["N1"]) if nt1.is_purine() else (nt1["C6"], nt1["N3"])
        b1dir, b2dir = (nt2["N9"], nt2["N1"]) if nt2.is_purine() else (nt2["C6"], nt2["N3"])
        dir_dist = a1dir.dist(b1dir)
        min_dir_dist = a1dir.dist(a2dir) + b1dir.dist(b2dir)
        
        if verbose: print(f"{dir_dist:.4f} - direction distance")
        if dir_dist < min_dir_dist:
            if verbose: print(f"FAILED - must be > molecule lengths = {min_dir_dist:.4f}")
            return False

        # angle between direction atoms
        # dirv1 = a2dir.coords - a1dir.coords
        # dirv2 = b2dir.coords - b1dir.coords
        # dirv1 /= np.linalg.norm(dirv1)
        # dirv2 /= np.linalg.norm(dirv2)
        # dir_angle = np.arccos(np.dot(dirv2, dirv1)) * 180 / np.pi
        # if verbose: print(f"{dir_angle:.4f} - direction angle")
        # if dir_angle < 90:
        #     if verbose: print("FAILED")
        #     return False
        
        # angle between base normal vectors
        norm1 = nt1.base_normal_vec()
        norm2 = nt2.base_normal_vec()
        normals_angle = np.arccos(abs(np.dot(norm1, norm2))) * 180 / np.pi
        
        if verbose: print(f"{normals_angle:.4f} - angle between base normal vectors")
        if normals_angle > NORMALS_ANGLE_THRESHOLD: # 67 min required
            if verbose: print(f"FAILED - must be < {NORMALS_ANGLE_THRESHOLD}")
            return False

        # vertical plane separation distance
        if normals_angle<COPLANAR_ANGLE_THRESHOLD: # check if bases are coplanar
            if verbose: print(f"Use vertical plane separation distance filter (norm angle < 25)")
            orig_dirv = origin1.coords - origin2.coords
            vert_dist1 = np.abs(np.dot(norm1, orig_dirv))
            vert_dist2 = np.abs(np.dot(norm2, orig_dirv))
            vert_dist = (vert_dist1 + vert_dist2) / 2
            
            if verbose: print(f"{vert_dist:.4f} - vertical plane separation distance")
            if vert_dist > 2.5:
                if verbose: print("FAILED")
                return False

        
        return True

    
    def _add_h_bonds(self, 
                     bonds: list, 
                     donor: NucleicAcidResidue, 
                     acceptor: NucleicAcidResidue, 
                     approximate_hs: bool
                    ): 
        
        datoms = DONOR_ACCEPTOR_GROUPS[donor.mname.lstrip('D')]["donors"]
        aatoms = DONOR_ACCEPTOR_GROUPS[acceptor.mname.lstrip('D')]["acceptors"]
        
        for hdname, dname in datoms:
            need_approximate = False
            if hdname not in donor:
                if not approximate_hs:
                    raise ValueError(f"Residue {donor.mname} {donor.mnum} does not contain hydrogen in '{dname}' donor group.")

                need_approximate = True
            
            for acceptor_atom_name in aatoms:
                acceptor_atom = acceptor[acceptor_atom_name]
                
                if need_approximate:
                    dist = acceptor_atom.dist(donor[dname]) - APPROXIMATE_H_BOND_DIST / 2
                else:
                    dist = acceptor_atom.dist(donor[hdname])

                eps = np.sqrt(EPSILON['H'] * EPSILON[acceptor_atom.element])
                sigma = (SIGMA['H'] + SIGMA[acceptor_atom.element]) / 2
                sigma6 = sigma**6
                sigma12 = sigma6**2
                dist6 = dist**6
                dist12 = dist6**2
                e_lj = 4*eps*(sigma12/dist12 - sigma6/dist6)

                if donor.is_rna():
                    qd = RNA_CHARGE[donor.mname][hdname]
                else:
                    qd = DNA_CHARGE[donor.mname.lstrip('D')][hdname]

                if acceptor.is_rna():
                    qa = RNA_CHARGE[acceptor.mname][acceptor_atom_name]
                else:
                    qa = DNA_CHARGE[acceptor.mname.lstrip('D')][acceptor_atom_name]
                
                e_c = COULOMB_CONST * qa * qd / dist
                
                bond = HBond(donor_res_name = donor.mname,
                               acceptor_res_name = acceptor.mname,
                               donor_name = hdname,
                               acceptor_name = acceptor_atom_name,
                               dist = dist,
                               approximated = need_approximate,
                               bond_e = e_lj+e_c,
                               LJ_e = e_lj,
                               C_e = e_c)
                bonds.append(bond)

    
    def filter_h_bonds(self, bonds: List[HBond]) -> List[HBond]:
        return list(filter(lambda x: x.dist <= H_BOND_DISTANCE_CUTOFF, bonds))
        
    def calculate_h_bonds(self, 
                          nt1: NucleicAcidResidue, 
                          nt2: NucleicAcidResidue, 
                          approximate_hs: bool = False
                         ) -> List[HBond]:
        
        bonds = []
        self._add_h_bonds(bonds, nt1, nt2, approximate_hs)
        self._add_h_bonds(bonds, nt2, nt1, approximate_hs)
        bonds = sorted(bonds, key=lambda x: x.bond_e)
        bonds = self.filter_h_bonds(bonds)
        
        return bonds
        
        
    def get_ss_energy_matrix(self, 
                             approximate_hs: bool = False,
                             verbose: bool = False
                            ) -> np.ndarray:
        
        E = np.full((len(self), len(self)), np.inf, dtype=np.float32)

        for i in range(len(self)):
            nt1 = self[i]
            for j in range(i+3, len(self)):
                nt2 = self[j]

                if not self.can_form_pair(nt1, nt2, verbose):
                    continue

                h_bonds = self.calculate_h_bonds(nt1, nt2, approximate_hs)
                bond_energy = sum([b.bond_e for b in h_bonds])

                if len(h_bonds)==0:
                    continue
                    
                if h_bonds[0].bond_e > MIN_H_ENERGY_THRESHOLD:
                    if verbose:
                        print(f"Strongest H bond - {h_bonds[0].bond_e:.4f} is too week")
                    continue
                
                E[i, j] = bond_energy
                E[j, i] = bond_energy
                
        return E












