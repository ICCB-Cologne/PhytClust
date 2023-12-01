import numpy as np
import pandas as pd


def pairwise_distances(tree, mode="terminals", as_dataframe=False, mrca=None):
    if mrca:
        clade = next((cl for cl in tree.find_clades() if cl.name == mrca), None)
        if clade is None:
            raise ValueError(f"No node found with name {mrca}")
        terms = list(clade.get_terminals())
    elif mode == "terminals":
        terms = list(tree.get_terminals())
    elif mode == "nonterminals":
        terms = list(tree.get_nonterminals())
    elif mode == "all":
        terms = list(tree.get_terminals()) + list(tree.get_nonterminals())
    else:
        print("Invalid mode")
        return None

    n = len(terms)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            bl = tree.distance(terms[i], terms[j])
            dist[i, j] = bl
            dist[j, i] = bl

    if as_dataframe:
        dist_df = pd.DataFrame(
            dist,
            index=[term.name for term in terms],
            columns=[term.name for term in terms],
        )
        return dist_df
    else:
        return dist
