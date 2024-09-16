from typing import Union, List, Iterable, Tuple
import numpy as np



def align(
            target: np.ndarray,
            aligned: np.ndarray,
            target_indices: Iterable[int],
            aligned_indices: Iterable[int],
            target_origin_idx: int = None,
            aligned_origin_idx: int = None
         ) -> np.ndarray:
    """
    https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    
    Aligns 'aligned' onto 'target' based on correspondence indices. 
    If origin index is None, a mean of aligning point coordinates will be used. 
    :param target: target coords
    :param aligned: aligned coords
    :param target_indices: target coords indices for aligning
    :param aligned_indices: aligned coords indices for aligning
    :param target_origin_idx: index of target point which coordinates will be used for origin shift
    :param aligned_origin_idx: index of aligned point which coordinates will be used for origin shift
    """
    
    if len(target_indices)!=len(aligned_indices):
        raise ValueError(f"Number of target and aligned indices must be the same.")
    
    a = target[target_indices]
    a_shift = a.mean(0) if (target_origin_idx is None) else target[target_origin_idx]
    a -= a_shift
    
    b = aligned[aligned_indices]
    b_shift = b.mean(0) if (aligned_origin_idx is None) else aligned[aligned_origin_idx]
    b -= b_shift
    
    C = a.T@b # dimension correlation matrix
    U, S, Vt = np.linalg.svd(C)
    R = np.dot(U, Vt).T # columns oriented matrix

    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = np.dot(U, Vt).T
    
    aligned = np.dot((aligned - b_shift), R) + a_shift
    return aligned