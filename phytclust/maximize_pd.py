from phytclust import PhytClust
from phytclust import plot_tree
from phytclust import pairwise_distances
import io
import os
from Bio import Phylo


# # without clusters
# def map_terminal_to_internal(tree):
#     terminal_to_internal = {}

#     def traverse(clade, path):
#         if not clade.clades:
#             terminal_to_internal[clade.name] = path
#         else:
#             new_path = path + [clade]
#             for child in clade.clades:
#                 traverse(child, new_path)

#     traverse(tree.root, [])

#     return terminal_to_internal


# def rank_terminal_nodes(tree, num_species=None, outgroup=None, clusters=None):
#     terminal_to_internal = map_terminal_to_internal(tree)

#     terminal_distances = {
#         terminal: tree.distance(terminal)
#         for terminal in terminal_to_internal
#         if terminal != outgroup
#     }

#     subtracted_nodes = {terminal: set() for terminal in terminal_to_internal}
#     choice_ranking = []
#     output_lists = []

#     while terminal_distances:
#         sorted_terminals = sorted(
#             terminal_distances.items(), key=lambda x: x[1], reverse=True
#         )

#         chosen_terminal, _ = sorted_terminals[0]

#         same_distance_terminals = [
#             terminal
#             for terminal, distance in sorted_terminals
#             if distance == terminal_distances[chosen_terminal]
#         ]
#         if len(same_distance_terminals) > 1:
#             print(
#                 f"Terminals with the same largest distance: {same_distance_terminals}"
#             )

#         print(
#             f"Chosen terminal: {chosen_terminal}, Distance: {terminal_distances[chosen_terminal]}"
#         )

#         shared_nodes = set(terminal_to_internal[chosen_terminal])
#         for terminal in list(
#             terminal_distances
#         ):  # Use list to avoid runtime modification issues
#             if terminal != chosen_terminal:
#                 for node in shared_nodes:
#                     if (
#                         node in terminal_to_internal[terminal]
#                         and node.branch_length is not None
#                         and node not in subtracted_nodes[terminal]
#                     ):
#                         print(f"Subtracting {node.branch_length} from {terminal}")
#                         terminal_distances[terminal] -= node.branch_length
#                         subtracted_nodes[terminal].add(node)

#         choice_ranking.append(
#             (chosen_terminal, terminal_distances.pop(chosen_terminal))
#         )
#         output_lists.append(choice_ranking.copy())

#         if num_species is not None and len(choice_ranking) >= num_species:
#             break

#     # Check for ties and mark them
#     previous_distance = None
#     for i, (terminal, distance) in enumerate(choice_ranking):
#         if distance == previous_distance:
#             choice_ranking[i] = (terminal, distance, "tie")
#         else:
#             choice_ranking[i] = (terminal, distance)
#         previous_distance = distance

#     return output_lists


# def _get_output_maximizing_pd(ranked_nodes):
#     sum_distances = sum(distance for _, distance in ranked_nodes)

#     node_names = [name for name, _ in ranked_nodes]

#     maximizing_pd_output = f"Maximizing PD to get {sum_distances}"
#     chosen_leaves_output = f"Chosen leaves: {', '.join(node_names)}"

#     return maximizing_pd_output, chosen_leaves_output


# def maximize_pd(tree, num_species=None, outgroup=None):
#     ranked_nodes_lists = rank_terminal_nodes(tree, num_species, outgroup)
#     outputs = []
#     for ranked_nodes in ranked_nodes_lists:
#         maximizing_pd_output, chosen_leaves_output = _get_output_maximizing_pd(
#             ranked_nodes
#         )
#         outputs.append((maximizing_pd_output, chosen_leaves_output))
#     return outputs

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


def rank_terminal_nodes(tree, num_species=None):
    terminal_to_internal = map_terminal_to_internal(tree)

    terminal_distances = {
        terminal: tree.distance(terminal) for terminal in terminal_to_internal
    }

    subtracted_nodes = {terminal: set() for terminal in terminal_to_internal}
    choice_ranking = []
    output_lists = []

    while terminal_distances:
        sorted_terminals = sorted(
            terminal_distances.items(), key=lambda x: x[1], reverse=True
        )

        chosen_terminal, chosen_distance = sorted_terminals[0]

        same_distance_terminals = [
            terminal
            for terminal, distance in sorted_terminals
            if distance == chosen_distance
        ]
        if len(same_distance_terminals) > 1:
            print(
                f"Terminals with the same largest distance: {same_distance_terminals}"
            )

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
                        terminal_distances[terminal] -= node.branch_length
                        subtracted_nodes[terminal].add(node)

        choice_ranking.append(
            (chosen_terminal, terminal_distances.pop(chosen_terminal))
        )
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


# def rank_terminal_nodes(tree, num_species=None, outgroup=None, clusters=None):
#     terminal_to_internal = map_terminal_to_internal(tree)

#     terminal_distances = {
#         terminal: tree.distance(terminal)
#         for terminal in terminal_to_internal
#         if terminal != outgroup
#     }

#     subtracted_nodes = {terminal: set() for terminal in terminal_to_internal}
#     choice_ranking = []
#     output_lists = []

#     best_per_cluster = {}
#     if clusters:
#         for terminal, cluster_id in clusters.items():
#             if terminal.name in terminal_distances:  # Match by name
#                 distance = terminal_distances[terminal.name]
#                 if (
#                     cluster_id not in best_per_cluster
#                     or distance > best_per_cluster[cluster_id][1]
#                 ):
#                     best_per_cluster[cluster_id] = (terminal.name, distance)

#         choice_ranking = list(best_per_cluster.values())

#     else:
#         while terminal_distances:
#             sorted_terminals = sorted(
#                 terminal_distances.items(), key=lambda x: x[1], reverse=True
#             )

#             chosen_terminal, _ = sorted_terminals[0]

#             same_distance_terminals = [
#                 terminal
#                 for terminal, distance in sorted_terminals
#                 if distance == terminal_distances[chosen_terminal]
#             ]
#             if len(same_distance_terminals) > 1:
#                 print(
#                     f"Terminals with the same largest distance: {same_distance_terminals}"
#                 )

#             print(
#                 f"Chosen terminal: {chosen_terminal}, Distance: {terminal_distances[chosen_terminal]}"
#             )

#             shared_nodes = set(terminal_to_internal[chosen_terminal])
#             for terminal in list(terminal_distances):
#                 if terminal != chosen_terminal:
#                     for node in shared_nodes:
#                         if (
#                             node in terminal_to_internal[terminal]
#                             and node.branch_length is not None
#                             and node not in subtracted_nodes[terminal]
#                         ):
#                             print(f"Subtracting {node.branch_length} from {terminal}")
#                             terminal_distances[terminal] -= node.branch_length
#                             subtracted_nodes[terminal].add(node)

#             choice_ranking.append(
#                 (chosen_terminal, terminal_distances.pop(chosen_terminal))
#             )
#             output_lists.append(choice_ranking.copy())

#             if num_species is not None and len(choice_ranking) >= num_species:
#                 break

#     return choice_ranking if clusters else output_lists


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
    # Map each terminal node to its path from the root

    # Initialize the list of representative species
    representatives = []

    # For each cluster...
    for cluster in set(clusters.values()):
        # Get the species in the cluster
        species_in_cluster = [
            species for species, clust in clusters.items() if clust == cluster
        ]

        # If the cluster has only one species, that species is the representative
        if len(species_in_cluster) == 1:
            representatives.append(species_in_cluster[0])
        else:
            # Otherwise, find the MRCA of the species in the cluster
            mrca = tree.common_ancestor(species_in_cluster)

            # Create a new tree with the MRCA as the root
            sub_tree = Phylo.BaseTree.Tree(root=mrca, rooted=True)

            # Select the representative species from the sub-tree
            ranked_species = rank_terminal_nodes(sub_tree, num_species=1)
            representatives.append(ranked_species[0][0][0])

    return representatives
