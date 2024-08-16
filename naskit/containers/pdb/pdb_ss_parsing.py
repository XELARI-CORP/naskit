import numpy as np

from ...parse_na import NA
from ..nucleic_acid import NucleicAcid



class SSParsing:
    
    def to_na(self):
        """
        10.1093/nar/gkv716
        https://academic.oup.com/nar/article/43/21/e142/2468098
        """
        
        energy_matrix = self.get_ss_energy_matrix()
        adj = self.parse_ss_adjacency(energy_matrix, threshold=4.5)
        na = NucleicAcid.from_adjacency(adj, seq=self.seq)
        return na

    
    def parse_ss_adjacency(self, M, threshold=0.5):
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


    def _origin_dist(self, nt1, nt2):
        o1 = nt1["N9"] if nt1.is_purine() else nt1["N1"]
        o2 = nt2["N9"] if nt2.is_purine() else nt2["N1"]
        return o1.dist(o2)
        
        
    def get_ss_energy_matrix(self):
        E = np.full((len(self), len(self)), np.inf, dtype=np.float32)

        for i in range(len(self)):
            nt1 = self[i]
            for j in range(i+3, len(self)):
                nt2 = self[j]

                # distance between the two origins
                if self._origin_dist(nt1, nt2) > 15:
                    continue

                # angle between the base normal vectors
                norm1 = nt1.base_normal_vec()
                norm2 = nt2.base_normal_vec()
                normals_angle = np.arccos(abs(np.dot(norm1, norm2)))
                if normals_angle > (65*np.pi/180):
                    continue
                
                # vertical separation between the base planes
                
                # absence of stacking between the two bases
                
                # presence of at least one hydrogen bond
                
                d = nt1.dist(nt2).min()
                E[i, j] = d
                E[j, i] = d
                
        return E












