from phytclust import PhytClust
from phytclust import plot_tree
from phytclust import pairwise_distances
import io
import os
from Bio import Phylo
#to be added: finidng minimal pd

# without clusters
def map_terminal_to_internal(tree):
    terminal_to_internal = {}

    def traverse(clade, path):
        if not clade.clades:
            terminal_to_internal[clade.name] = path
        else:
            new_path = path + [clade]
            for child in clade.clades:
                traverse(child, new_path)

    traverse(tree.root, [])

    return terminal_to_internal

import heapq


def rank_terminal_nodes(tree, num_species=None):
    terminal_to_internal = map_terminal_to_internal(tree)

    terminal_distances = {
        terminal: -tree.distance(terminal) for terminal in terminal_to_internal
    }

    subtracted_nodes = {terminal: set() for terminal in terminal_to_internal}
    choice_ranking = []
    output_lists = []

    # Use a priority queue for terminal_distances
    terminal_queue = [(distance, terminal) for terminal, distance in terminal_distances.items()]
    heapq.heapify(terminal_queue)

    while terminal_queue:
        chosen_distance, chosen_terminal = heapq.heappop(terminal_queue)
        chosen_distance = -chosen_distance  # Convert back to positive distance

        print(f"Chosen terminal: {chosen_terminal}, Distance: {chosen_distance}")

        shared_nodes = set(terminal_to_internal[chosen_terminal])
        for terminal in terminal_distances:
            if terminal != chosen_terminal:
                terminal_nodes = set(terminal_to_internal[terminal])
                for node in shared_nodes:
                    if (
                        node in terminal_nodes
                        and node.branch_length is not None
                        and node not in subtracted_nodes[terminal]
                    ):
                        print(f"Subtracting {node.branch_length} from {terminal}")
                        terminal_distances[terminal] += node.branch_length
                        subtracted_nodes[terminal].add(node)

        choice_ranking.append((chosen_terminal, chosen_distance))
        output_lists.append(choice_ranking.copy())

        if num_species is not None and len(choice_ranking) >= num_species:
            break

    # Check for ties and mark them
    previous_distance = None
    for i, (terminal, distance) in enumerate(choice_ranking):
        if distance == previous_distance:
            choice_ranking[i] = (terminal, distance, "tie")
        else:
            choice_ranking[i] = (terminal, distance)
        previous_distance = distance

    return output_lists

def _get_output_maximizing_pd(ranked_nodes):
    sum_distances = sum(distance for _, distance in ranked_nodes)

    node_names = [name for name, _ in ranked_nodes]

    maximizing_pd_output = f"Maximizing PD to get {sum_distances}"
    chosen_leaves_output = f"Chosen leaves: {', '.join(node_names)}"

    return maximizing_pd_output, chosen_leaves_output


def maximize_pd(tree, num_species=None, outgroup=None, clusters=None):
    ranked_nodes_lists = rank_terminal_nodes(
        tree, num_species, outgroup, clusters=clusters
    )
    outputs = []
    for ranked_nodes in ranked_nodes_lists:
        maximizing_pd_output, chosen_leaves_output = _get_output_maximizing_pd(
            ranked_nodes
        )
        outputs.append((maximizing_pd_output, chosen_leaves_output))
    return outputs


def select_representative_species(tree, clusters):
    cluster_to_species = {}
    for species, cluster in clusters.items():
        if cluster not in cluster_to_species:
            cluster_to_species[cluster] = []
        cluster_to_species[cluster].append(species)

    representatives = []

    for cluster, species_in_cluster in cluster_to_species.items():
        if len(species_in_cluster) == 1:
            representatives.append(species_in_cluster[0])
        else:
            mrca = tree.common_ancestor(species_in_cluster)

            sub_tree = Phylo.BaseTree.Tree(root=mrca, rooted=True)
            ranked_species = rank_terminal_nodes(sub_tree, num_species=1)
            representatives.append(ranked_species[0][0][0])

    return representatives
