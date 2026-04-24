# How It Works

This page explains what PhytClust does under the hood. No equations required — just intuition about trees and clustering.

---

## The problem

You have a rooted phylogenetic tree and you want to split its leaves into groups. Sounds simple, but there's a catch: you want each group to be a **clade** — a complete subtree where every descendant is included. This property is called **monophyly**.

Why does monophyly matter? Because a clade has a clear evolutionary interpretation: it's "everything descended from this ancestor." If your clusters aren't monophyletic, they're mixing lineages in ways that are hard to justify biologically.

The question, then, is: *of all possible ways to partition a tree into k monophyletic groups, which one is best?*

## Where cuts can happen

Think of the tree as a branching structure where each internal edge is a potential "cut point." If you cut an edge, you separate the subtree below it from the rest. A clustering is just a set of cuts that divides all leaves into *k* groups.

Not every combination of cuts produces monophyletic clusters — you need the cuts to be *consistent*. PhytClust's dynamic programming guarantees this by working bottom-up through the tree, so it only ever considers valid partitions.

## The objective: minimize within-cluster branch length

Given a valid partition, PhytClust measures quality by summing up the branch lengths *within* each cluster. Shorter within-cluster branches mean tighter, more cohesive groups. The DP finds the partition that minimizes this total for a given *k*.

Intuitively: long branches between clusters (where the cuts happen) are good — they separate well-differentiated lineages. Short branches within clusters are good — they keep similar taxa together.

## Dynamic programming on the tree

The algorithm works bottom-up:

1. **Start at the leaves.** Each leaf is trivially a cluster of size 1.
2. **At each internal node**, decide how to allocate its *children* across clusters. The DP table records the best achievable cost for "partition the subtree rooted here into *j* groups" for every valid *j*.
3. **At the root**, read off the optimal partition for your target *k*.

Because the DP visits each node once and considers all valid splits at each node, the result is **exact** — not an approximation or heuristic. For a binary tree with *n* leaves, the time complexity is roughly O(*n* · *k*²), which is fast enough for trees with thousands of leaves.

## The score curve: finding good *k* values

If you don't know the right *k*, PhytClust computes the optimal cost for every *k* from 2 up to some maximum. Plotting cost against *k* gives a **score curve**.

The score curve typically drops sharply at *k* values where the tree has natural breakpoints — places where a long branch separates two subtrees. These sharp drops appear as **peaks** in the rate of improvement.

PhytClust uses scipy's peak detection to find these peaks, then ranks them by **prominence** — how much a peak stands out from its neighbors. The most prominent peaks correspond to the *k* values where the tree structure most strongly suggests a cluster boundary.

### Raw vs. adjusted ranking

By default, PhytClust uses "adjusted" ranking, which accounts for how the number of outlier clusters changes with *k*. This prevents the ranker from favoring *k* values that only look good because they produce many tiny noise clusters. You can control the balance between raw and adjusted ranking with the `lambda_weight` parameter.

## Polytomies

A polytomy (multifurcation) is an internal node with more than two children. They're common in real trees — either because the true branching order is unresolved, or because multiple speciation events happened in rapid succession.

PhytClust handles polytomies in two ways:

- **Hard mode** (default): each child of a polytomy goes entirely into one cluster. This is fast and deterministic — the DP just has more children to allocate at that node.

- **Soft mode**: children can be partially merged across clusters. This is more flexible but computationally expensive — the number of possible allocations grows exponentially with the node's degree. A safety guardrail (`soft_polytomy_max_degree`, default 18) falls back to hard mode for very high-degree nodes.

For most trees, hard mode works well. Soft mode is worth trying when you have moderate-degree polytomies (say, 5-15 children) and suspect that the hard constraint is forcing unnatural groupings.

## Outlier handling

Small clusters are a fact of life in real phylogenies — isolated taxa on long branches, contaminant sequences, or just biological noise. PhytClust gives you two ways to deal with them:

- **Hard constraint** (`min_cluster_size`): the DP enforces a minimum cluster size during optimization. If a *k* can't be achieved without violating this, it's simply not returned.

- **Soft handling** (`outlier_size_threshold` + `prefer_fewer_outliers`): clusters below the threshold are *marked* as outliers (cluster ID = -1) in the output, but they still exist. With `prefer_fewer_outliers`, the DP slightly favors solutions that concentrate noise into fewer groups rather than scattering it.

## Zero-length edges

Some trees have internal edges of length zero — typically from collapsing poorly-supported nodes or from identical sequences. By default, PhytClust treats these like any other edge and may split at them if the DP says so.

With `--no-split-zero-length`, the DP is forbidden from cutting at zero-length edges. This keeps groups of near-identical taxa together, which makes sense when the zero-length branch has no biological meaning.

## Summary

| Concept | What it means in PhytClust |
|---------|---------------------------|
| **Monophyly** | Every cluster is a complete clade — guaranteed, not post-hoc filtered |
| **Objective** | Minimize total within-cluster branch length |
| **DP** | Exact bottom-up algorithm on the tree; visits each node once |
| **Score curve** | Cost vs. *k*; peaks indicate natural cluster boundaries |
| **Peak ranking** | Prominence-based, with optional outlier adjustment |
| **Polytomies** | Hard mode (fast, per-child) or soft mode (flexible, exponential) |
| **Outliers** | Hard constraint (min size) or soft marking (threshold + preference) |
