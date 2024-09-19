from typing import Optional, Union
import random
from ..containers.nucleic_acid import NucleicAcid
from ..parse_na import NA



def find_free_motifs(na: NucleicAcid):
    last_nb_has_bond = True
    motifs = []
    for i in range(len(na)):
        if na.complnb(i) is None:
            if last_nb_has_bond:
                motifs.append([])
            motifs[-1].append(i)
            last_nb_has_bond = False
        else:
            last_nb_has_bond = True
    
    return motifs


def add_helix(na: NucleicAcid, n: int):
    motifs = [m for m in find_free_motifs(na) if len(m)>=n]
    if len(motifs)==0:
        return False
    elif len(motifs)>1:
        first_motif = motifs[0]
        second_motif = random.choice(motifs[1:])
    else:
        m = motifs[0]
        if len(m)//2 < n:
            return False
        first_motif = m[:len(m)//2]
        second_motif = m[len(m)//2:]

    l = len(first_motif)
    a = first_motif[(l-n)//2 : (l-n)//2 + n]
    l = len(second_motif)
    b = second_motif[(l-n)//2 : (l-n)//2 + n]

    for o, e in zip(a, b[::-1]):
        na.join(o,e)
    
    return True


def compl_ratio(na: NucleicAcid) -> float:
    return 2*len(na.pairs)/len(na)


def generate_ss(na: Union[NucleicAcid, str],
                min_helix_size: int = 2,
                max_helix_size: int = 6,
                max_compl_ratio: float = 0.5,
                patience: int = 1,
                seed: Optional[int] = None
               ) -> NucleicAcid:
    """
    Generates random secondary structure.

    :param na: NucleicAcid object or sequence string.
    :param min_helix_size: minimal size of generated helix. Default = 2.
    :param max_helix_size: maximum size of generated helix. Default = 6.
    :param max_compl_ratio: maximum fraction of nucleotides to form bonds [0, 1]. Default = 0.5.
    :param patience: number of tries to fit helix of random size before function return. Default = 1.
    :param seed: random seed. Default = None - every function call gives different result.
    
    :return: new NucleicAcid object with generated structure.
    """

    na = NA(na) if isinstance(na, str) else na.copy()
    if seed is not None:
        random.seed(seed)

    attempts = patience
    while (compl_ratio(na)<max_compl_ratio and attempts>0):
        n = random.randint(min_helix_size, max_helix_size)
        success = add_helix(na, n)
        if not success:
            attempts-=1
            continue
        attempts = patience

    return na







