from typing import Union, List, Iterable, Tuple
import numpy as np



def align(
            target: np.ndarray,
            aligned: np.ndarray,
            target_indices: Iterable[int],
            aligned_indices: Iterable[int]
         ) -> np.ndarray:
    """
    https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    
    Aligns 'aligned' onto 'target' based on correspondence indices.
    :param target: target coords
    :param aligned: aligned coords
    :param target_indices: target coords indices for aligning
    :param aligned_indices: aligned coords indices for aligning
    """
    
    if len(target_indices)!=len(aligned_indices):
        raise ValueError(f"Number of target and aligned indices must be the same.")
    
    a = target[target_indices]
    a_mean = a.mean(0)
    a -= a_mean
    
    b = aligned[aligned_indices]
    b_mean = b.mean(0)
    b -= b_mean
    
    C = a.T@b # dimension correlation matrix
    U, S, Vt = np.linalg.svd(C)
    R = np.dot(U, Vt).T # columns oriented matrix

    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = np.dot(U, Vt).T
    
    aligned = np.dot(aligned, R) + a_mean
    return aligned