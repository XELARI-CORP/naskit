from typing import Union, List, Iterable, Tuple
import numpy as np



def align(
            a: np.ndarray, 
            b: np.ndarray, 
            correspondence: Iterable[Tuple[int, int]]
         ):
    """
    Aligns b onto a based on correspondence indices.
    :param a: target coords
    :param b: aligned coords
    :param correspondence: pairs of indices (ai, bi) for aligning
    """
    
    b = b - b.mean(0) + a.mean(0)
    return b