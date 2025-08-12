from __future__ import annotations

from typing import Dict, List, Optional, Tuple
import heapq

from Bio.Phylo.BaseTree import Tree


#  terminal -> path of internal nodes
def map_terminal_to_internal(tree: Tree) -> Dict[str, List]:
    terminal_to_internal: Dict[str, List] = {}

    def traverse(clade, path):
        if not clade.clades:
            terminal_to_internal[getattr(clade, "name", None)] = path
        else:
            new_path = path + [clade]
            for child in clade.clades:
                traverse(child, new_path)

    traverse(tree.root, [])
    # Filter out None-named leaves if any
    return {k: v for k, v in terminal_to_internal.items() if k is not None}


# PD pick
def _get_output_maximizing_pd(ranked_nodes: List[Tuple[str, float]]) -> Tuple[str, str]:
    sum_distances = sum(distance for _, distance in ranked_nodes)
    node_names = [name for name, _ in ranked_nodes]
    maximizing_pd_output = f"Maximizing PD to get {sum_distances}"
    chosen_leaves_output = f"Chosen leaves: {', '.join(node_names)}"
    return maximizing_pd_output, chosen_leaves_output


def maximize_pd(tree: Tree, num_species: Optional[int] = None) -> List[Tuple[str, str]]:
    """
    Convenience wrapper: greedy selection by 'all' distances, maximizing sum.
    Returns human-readable summaries for each incremental selection set.
    """
    ranked_nodes = rank_terminal_nodes(
        tree, num_species=num_species, mode="maximize", distance_ref="all"
    )
    outputs: List[Tuple[str, str]] = []
    maximizing_pd_output, chosen_leaves_output = _get_output_maximizing_pd(ranked_nodes)
    outputs.append((maximizing_pd_output, chosen_leaves_output))
    return outputs


# Distances used for ranking
def compute_species_distance(
    tree: Tree,
    terminals: List[str],
    *,
    distance_ref: str = "all",
) -> Dict[str, float]:
    """
    Compute a distance measure for each terminal.

    distance_ref:
        - 'all': sum of distances from this terminal to all other terminals  (O(n^2))
        - 'mrca': distance from this terminal to the MRCA of all 'terminals'

    Returns: dict {terminal_name: distance_value}

    NOTE: For large trees, 'all' is O(n^2) over terminals since Bio.Phylo.distance
    is invoked pairwise. Consider caching paths or using a distance matrix if this
    becomes slow at scale.
    """
    distances: Dict[str, float] = {}

    if distance_ref == "all":
        for i, t in enumerate(terminals):
            s = 0.0
            for j, other in enumerate(terminals):
                if i == j:
                    continue
                s += float(tree.distance(t, other))
            distances[t] = s

    elif distance_ref == "mrca":
        if not terminals:
            return {}
        mrca = tree.common_ancestor(terminals)
        for t in terminals:
            distances[t] = float(tree.distance(t, mrca))

    else:
        raise ValueError("distance_ref must be 'all' or 'mrca'")

    return distances


# Greedy ranking
def rank_terminal_nodes(
    tree: Tree,
    num_species: Optional[int] = None,
    *,
    mode: str = "maximize",
    distance_ref: str = "all",
) -> List[Tuple[str, float]]:
    """
    Rank/Select terminals by a simple scalar "distance" measure.

    mode:
        - 'maximize': pick highest distance first (greedy)
        - 'minimize': pick lowest distance first

    Returns a list of (terminal_name, distance_value) in chosen order.
    """
    terminals = [t.name for t in tree.get_terminals() if t.name]
    base = compute_species_distance(tree, terminals, distance_ref=distance_ref)

    if mode not in {"maximize", "minimize"}:
        raise ValueError("mode must be 'maximize' or 'minimize'")

    if mode == "maximize":
        queue = [(-base[t], t) for t in terminals]
    else:
        queue = [(base[t], t) for t in terminals]
    heapq.heapify(queue)

    selected: List[Tuple[str, float]] = []
    seen = set()
    while queue and (num_species is None or len(selected) < num_species):
        priority, terminal = heapq.heappop(queue)
        if terminal in seen:
            continue
        seen.add(terminal)
        actual = -priority if mode == "maximize" else priority
        selected.append((terminal, float(actual)))

    return selected


# One representative per cluster
def select_representative_species(
    tree: Tree,
    clusters: Dict[str, int],
    *,
    mode: str = "maximize",
    distance_ref: str = "mrca",
) -> List[str]:
    """
    Pick a single representative species per cluster using the chosen criterion.

    clusters: mapping {species_name -> cluster_id}
    """
    # group by cluster
    by_cluster: Dict[int, List[str]] = {}
    for sp, cid in clusters.items():
        by_cluster.setdefault(cid, []).append(sp)

    reps: List[str] = []

    for cid, species_list in by_cluster.items():
        if len(species_list) == 1:
            reps.append(species_list[0])
            continue

        # rank within each cluster by the chosen criterion
        # Use the full tree; MRCA is computed from the species_list subset
        mrca = tree.common_ancestor(species_list)

        # Build a temporary subtree rooted at MRCA for distance calc context
        sub_tree = Tree(root=mrca, rooted=True)

        ranked = rank_terminal_nodes(
            sub_tree,
            num_species=1,
            mode=mode,
            distance_ref=distance_ref,
        )
        reps.append(ranked[0][0] if ranked else species_list[0])

    return reps
