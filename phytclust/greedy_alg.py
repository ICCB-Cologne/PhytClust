import io
import os
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def calculate_total_potential_score(node, memo={}):
    if node in memo:
        return memo[node]
    total_score = node.self_cost + sum(
        calculate_total_potential_score(child, memo) for child in node.children
    )
    memo[node] = total_score
    return total_score


def split_tree(node, max_splits=3):
    nodes = [node]
    remaining_potential_scores = []
    splits_performed = 0
    states = []
    total_potential_score = calculate_total_potential_score(node)
    total = total_potential_score
    remaining_potential_scores.append(total_potential_score)

    while nodes and splits_performed < max_splits:
        # Sort nodes based on the absolute value of self_cost, descending
        nodes.sort(key=lambda n: abs(n.self_cost), reverse=True)

        # Select the node to split
        node_to_split = nodes.pop(0)
        total_potential_score -= node_to_split.self_cost
        remaining_potential_scores.append(total_potential_score)
        splits_performed += 1

        nodes.extend(node_to_split.children)

        # Update the states list with the current snapshot of nodes
        states.append([n.name for n in nodes])  # Storing node names as state

        # is_terminal = "Yes" if len(node_to_split.children) == 0 else "No"
        # print(
        #     f"Splitting at node: {node_to_split.name}, "
        #     f"Remaining Potential Score: {total_potential_score}, "
        #     f"Is Terminal (Diploid): {is_terminal}"
        # )
        # print(
        #     f"Reason for choice: Node selected based on being non-terminal and having the highest self_cost/total_recursive_cost ratio; "
        #     f"Ratio: {float(node_to_split.self_cost) / node_to_split.total_recursive_cost if node_to_split.total_recursive_cost else 0}"
        # )

    return remaining_potential_scores, states, total


def find_node_by_name(root, name):
    for clade in root.find_clades():
        if clade.name == name:
            return clade
    return None


def map_nodes(root, remaining_nodes_lists):
    all_mappings = []

    # Function to find a node by its name
    def find_node_by_name(node, name):
        if node.name == name:
            return node
        for (
            child
        ) in node.clades:  # clades is the correct attribute in Bio.Phylo for children
            result = find_node_by_name(child, name)
            if result is not None:
                return result
        return None

    for node_list in remaining_nodes_lists:
        mapping = {}
        unique_number = 0  # Reset the unique number for each list
        for node_name in node_list:
            node = find_node_by_name(root, node_name)
            if node:
                terminals = node.get_terminals()  # Get terminal nodes directly
                for terminal in terminals:
                    mapping[terminal.name] = unique_number
            unique_number += 1  # Increment for the next node in the list
        all_mappings.append(mapping)

    return all_mappings


def calculate_beta_scores(scores, total_terminal_nodes):
    beta_scores = []
    beta_list = []
    den_list = []
    N = total_terminal_nodes
    beta_1 = scores[0]
    beta_list.append(beta_1)

    for k, beta_k in enumerate(scores):
        if beta_k < 0.00005 or beta_k == 0:
            beta_scores.append(float("nan"))
            beta_list.append(beta_k)
        else:
            K = k + 1
            if K == 1:
                beta_scores.append(float("nan"))
                # beta_list.append(0)
            else:
                num = (beta_1 - beta_k) / beta_k
                den = (N - K) / (K - 1)
                ranking = beta_k / beta_1
                beta_score = num * den
                beta_scores.append(beta_score)
                beta_list.append(beta_k)
                den_list.append(ranking)

    return beta_scores, beta_list, den_list


def traverse_and_map(node, name_to_number_mapping):
    mapping = {}
    if node.name in name_to_number_mapping:
        mapping[node] = name_to_number_mapping[node.name]
    for child in node.clades:
        mapping.update(traverse_and_map(child, name_to_number_mapping))
    return mapping


def convert_to_clade_dicts(tree, output):
    clade_dicts = []

    # Create an initial dictionary with all nodes having a cluster number of 0
    initial_dict = traverse_and_map(
        tree.root, {node.name: 0 for node in tree.find_clades()}
    )
    clade_dicts.append(initial_dict)

    for cluster in output:
        clade_dict = traverse_and_map(tree.root, cluster)
        clade_dicts.append(clade_dict)
    return clade_dicts
