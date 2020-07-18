# Implement Felsenstein algorithm with the exemplary JC69 model.

Nucleotides: A, C, T, G represented as 'A', 'C' etc. Stationary distribution stored in the dictionary, eg q['A'] = .... 
Assume uniform distribution. Set the alpha parameter to 0.25.
The alignment is a list of strings x[0, ..., n-1], where x [i] is the i-th sequence of m length. Tree:
- nodes numbered from 0 to 2n2, nodes from 0 to n1 are leaves
- tree structure represented by the list s [0, ..., 2n-2] list of tuples, empty for leaves, pair of numbers for internal nodes
- edge lengths in the table t [0, ..., 2n-2]

## Functions:
pr(x, y, time) - returns the probability of transforming the nucleotide x into y on the edge of the length T (formula according to the sequence evolution model)

prL(k, a, u, x, t, s) - returns the position in the u-th sequence for the subtree of the root of the number k, where the nucleotide a 
stands in k (according to the Felsenstein algorithm P(Lk|a)).

prLBis(u, x, t, s) - returns the probability of a tree for u-th position (column); uses the prL.

treeLogLikelihood (x, t, s) - main procedure calculating the tree likelihood log for: x - alignment, t - edge lengths in the tree, 
s - father-sons relationship; uses prLBis.

