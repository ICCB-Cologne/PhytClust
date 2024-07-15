import io
import os
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Node:
    def __init__(self, clade):
        self.clade = clade
        self.children = [Node(child) for child in clade.clades]
        self.branch_cost = clade.branch_length if clade.branch_length is not None else 0
        self.self_cost = sum(child.terminal_child_count * child.branch_cost for child in self.children)
        self.terminal_child_count = self.calculate_terminal_child_count()
        self.total_recursive_cost = self.self_cost + sum(
            child.total_recursive_cost for child in self.children
        )
        # print(
        #     f"Node: {self.clade.name}, Self Cost: {self.self_cost}, Total Recursive Cost: {self.total_recursive_cost}"

    def calculate_terminal_child_count(self):
        if not self.children:
            return 1  # Leaf node
        return sum(child.calculate_terminal_child_count() for child in self.children)


# def calculate_self_cost(clade):
#     if not clade.clades:
#         return 0
#     return sum(
#         calculate_terminal_child_count(child)
#         * (child.branch_length if child.branch_length else 0)
#         for child in clade.clades
#     )


# def calculate_terminal_child_count(clade):
#     if not clade.clades:
#         return 1
#     return sum(calculate_terminal_child_count(child) for child in clade.clades)


# def calculate_total_recursive_cost(clade):
#     if not clade.clades:
#         return calculate_self_cost(clade)
#     return calculate_self_cost(clade) + sum(
#         calculate_total_recursive_cost(child) for child in clade.clades
#     )


def calculate_total_potential_score(node):  # totalBranchCost
    return node.self_cost + sum(
        calculate_total_potential_score(child) for child in node.children
    )


# def split_tree(root, max_splits=3):
#     nodes = [root]
#     remaining_potential_scores = []
#     splits_performed = 0
#     states = []
#     total_potential_score = calculate_total_recursive_cost(root)
#     remaining_potential_scores.append(total_potential_score)

#     while nodes and splits_performed < max_splits:
#         nodes.sort(key=lambda clade: abs(calculate_self_cost(clade)), reverse=True)
#         node_to_split = nodes.pop(0)
#         total_potential_score -= calculate_self_cost(node_to_split)
#         remaining_potential_scores.append(total_potential_score)
#         splits_performed += 1
#         nodes.extend(node_to_split.clades)
#         states.append([child.name for child in node_to_split.clades])

#     return remaining_potential_scores, states, total_potential_score


def split_tree(node, max_splits=3):
    nodes = [node]
    remaining_potential_scores = []
    splits_performed = 0
    states = []
    total_potential_score = calculate_total_potential_score(node)
    total = total_potential_score
    remaining_potential_scores.append(total_potential_score)

    # while nodes and splits_performed < max_splits:
    # Sort nodes based on the ratio of self_cost to total_recursive_cost, descending
    # Nodes with higher ratios are prioritized
    # nodes.sort(
    #     key=lambda n: (
    #         len(n.children) > 0,
    #         (
    #             float(n.self_cost) / n.total_recursive_cost
    #             if n.total_recursive_cost
    #             else 0
    #         ),
    #     ),
    #     reverse=True,
    # )
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

        is_terminal = "Yes" if len(node_to_split.children) == 0 else "No"
        print(
            f"Splitting at node: {node_to_split.name}, "
            f"Remaining Potential Score: {total_potential_score}, "
            f"Is Terminal (Diploid): {is_terminal}"
        )
        print(
            f"Reason for choice: Node selected based on being non-terminal and having the highest self_cost/total_recursive_cost ratio; "
            f"Ratio: {float(node_to_split.self_cost) / node_to_split.total_recursive_cost if node_to_split.total_recursive_cost else 0}"
        )

    return remaining_potential_scores, states, total


def find_node_by_name(root, name):
    for clade in root.find_clades():
        if clade.name == name:
            return clade
    return None


# def map_nodes(root, remaining_nodes_lists):
#     all_mappings = {}
#     for node_list in remaining_nodes_lists:
#         for node_name in node_list:
#             node = find_node_by_name(root, node_name)
#             if node:
#                 all_mappings[node.name] = node.terminal_child_count
#     return all_mappings


def map_nodes(root, remaining_nodes_lists):
    all_mappings = []

    # Function to recursively collect terminal nodes and assign them the same unique number
    def assign_terminal_numbers(node, mapping, unique_number):
        if not node.children:  # It's a terminal node
            mapping[node.name] = unique_number
        else:  # It's an internal node
            for child in node.children:
                assign_terminal_numbers(child, mapping, unique_number)

    # Function to find a node by its name
    def find_node_by_name(node, name):
        if node.name == name:
            return node
        for child in node.children:
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
                assign_terminal_numbers(node, mapping, unique_number)
            unique_number += 1  # Increment for the next node in the list
        all_mappings.append(mapping)

    return all_mappings


def calculate_beta_scores(scores, total_terminal_nodes):
    beta_scores = []
    num_list = []
    den_list = []
    N = total_terminal_nodes

    beta_1 = scores[0]

    for k in range(len(scores)):
        beta_k = scores[k]

        # Check if beta_k is very small (close to zero) or exactly zero
        if beta_k < 0.005 or beta_k == 0:
            beta_scores.append(float("nan"))
        else:
            K = k + 1
            # Ensure that we're not dividing by zero when calculating beta_score
            if K == 1:
                beta_scores.append(float("nan"))
            else:
                num = ((beta_1 - beta_k) / beta_k) 
                den = ((N - K) / (K - 1))
                beta_score = num*den
                beta_scores.append(beta_score)
                num_list.append(num)
                den_list.append(den)

    return beta_scores, num_list, den_list


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
