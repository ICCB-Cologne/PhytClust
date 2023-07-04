import numpy as np
from itertools import combinations
import pandas as pd

def calculate_pdm_from_tree(tree, sample_labels, normal_name='diploid'):
    '''This function takes a tree and creates a pairwise distance matrix using distances along the tree branches'''
    depths = {k.name: v for k, v in tree.depths().items()}
    dm = np.zeros((len(sample_labels), len(sample_labels)))

    paths = {t: np.array([x.name for x in tree.get_path(t)][:-1]) for t in sample_labels}

    for i, j in combinations(np.arange(len(sample_labels)), r=2):
        s1 = sample_labels[i]
        s2 = sample_labels[j]

        if s1 == 'diploid':
            cur_dist = depths[s2]
        elif s2 == 'diploid':
            cur_dist = depths[s1]
        else:
            cma = paths[s1][np.in1d(paths[s1], paths[s2])][-1]
            cur_dist = depths[s1] + depths[s2] - 2* depths[cma]
        dm[i, j] = cur_dist
        dm[j, i] = cur_dist    

    dm = pd.DataFrame(dm, index = sample_labels, columns = sample_labels)
    dm = dm.astype(float)

    return dm