#!/usr/bin/env python3
# Simulate trees with msprime, convert to a matrix of binary-vector encoded
# nodes, and calculate RNNI distance.
import msprime
import numpy as np


def node_encodings(ts):
    """
    Associate a binary-vector identifier to each node according to
    the identities of the leaves descended from the node.
    """
    # XXX: assumes samples == leaves
    n = ts.num_samples
    W = np.eye(n)
    f = lambda x: x
    T = ts.general_stat(
        W,
        f,
        output_dim=n,
        mode="node",
        windows="trees",
        strict=False,
        polarised=True,
        span_normalise=True,
    )
    T = T.astype(np.uint8)
    T = T[:, n:, :]  # exclude the leaves, we don't really need them
    return T


def bitwise_subset(a, b):
    """
    Return true if all bits in ``a`` are set in ``b``.
    """
    return all(a & b == a)


def rank(T, x):
    """
    Return the rank of node x in T (the index of the first node in T that
    contains all leaves under x). This is the same as:
    """
    for j, t in enumerate(T):
        if all(t & x == x):
            # T[j] is a common ancestor of the members of x
            return j
    raise ValueError(f"{x} not contained in {T}")


def rnni_findpath(T, R):
    """
    RNNI distance between trees T and R.
    Collienne & Gavryushkin (2020) https://arxiv.org/abs/2007.12307
    """
    d = 0
    T1 = T.copy()
    # path = [T1.copy()]
    for k in range(len(R) - 1):
        Ck = R[k]
        r = rank(T1, Ck)
        while r > k:
            v = T1[r]
            u = T1[r - 1]
            if bitwise_subset(u, v):
                # u is a child of v, so do a NNI to reduce rank(T1, Ck).
                w = v - u
                # Find children of u.
                # XXX: assumes binary trees
                leaves = np.nonzero(u)[0]
                if len(leaves) == 2:
                    # Both children are leaves.
                    x = np.zeros_like(u)
                    x[leaves[0]] = 1
                else:
                    for x in reversed(T1[: r - 1]):
                        if bitwise_subset(x, u):
                            # x is a child of u
                            break
                    else:
                        raise ValueError(f"{u} has no children in {T1}")
                y = u - x
                # Currently we have u = x + y, and the two NNI options are:
                #    u = x + w, or,
                #    u = y + w.
                # Only one of these will reduce rank(T1, Ck).
                if bitwise_subset(Ck, x + w):
                    T1[r - 1] = x + w
                else:
                    T1[r - 1] = y + w
            else:
                # Swap nodes v and u.
                T1[[r, r - 1]] = T1[[r - 1, r]]
            r -= 1  # Both operations reduce the rank by 1.
            d += 1
            # path.append(T1.copy())
    return d  # , path


if __name__ == "__main__":
    n = 50  # number of leaves
    # Simulate two tree sequences. We have recombination_rate=0 so each
    # "sequence" of trees will only be a single tree.
    ts1, ts2 = msprime.simulate(
        sample_size=n,
        recombination_rate=0,
        random_seed=1,
        num_replicates=2,
        Ne=1000,
    )
    tree1 = node_encodings(ts1)[0]
    tree2 = node_encodings(ts2)[0]
    d = rnni_findpath(tree1, tree2)
    print(d)
