"""Independent reference implementations used only for testing.

These build AMGs as NetworkX multigraphs directly from the DSB counts, without
using the partner maps of :class:`InitialConfiguration`, so they can check the
library rather than repeat it. Edge keys record the edge type, so parallel
edges of different types are kept apart.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence

import networkx as nx


def amg_multigraph(breaks: Sequence[int], rejoin_pairs) -> nx.MultiGraph:
    """The full AMG: telomeres, free ends, chromatin, DSB and rejoin edges.

    Chromosome ``i`` is laid out as ``T_i^L - e - e - ... - e - T_i^R`` where
    the free ends ``e`` are numbered consecutively over all chromosomes and
    consecutive ends alternate chromatin and DSB edges.
    """
    g = nx.MultiGraph()
    end = 0
    for i, b in enumerate(breaks):
        path = [("tel", i, "L")] + list(range(end, end + 2 * b)) + [("tel", i, "R")]
        end += 2 * b
        g.add_nodes_from(path)
        for pos in range(len(path) - 1):
            kind = "chromatin" if pos % 2 == 0 else "dsb"
            g.add_edge(path[pos], path[pos + 1], key=kind)
    for u, v in rejoin_pairs:
        g.add_edge(u, v, key="rejoin")
    return g


def exchange_multigraph(num_dsbs: int, rejoin_pairs) -> nx.MultiGraph:
    """Ξ: free ends with DSB edges (2j, 2j+1) and rejoin edges."""
    g = nx.MultiGraph()
    g.add_nodes_from(range(2 * num_dsbs))
    for j in range(num_dsbs):
        g.add_edge(2 * j, 2 * j + 1, key="dsb")
    for u, v in rejoin_pairs:
        g.add_edge(u, v, key="rejoin")
    return g


def reference_cycle_parts(num_dsbs: int, rejoin_pairs) -> tuple[int, ...]:
    """Cycle lengths ``l`` (edges / 2) of every component of Ξ, non-increasing."""
    g = exchange_multigraph(num_dsbs, rejoin_pairs)
    assert all(d == 2 for _, d in g.degree()), "Ξ must be 2-regular"
    parts = []
    for comp in nx.connected_components(g):
        edges = g.subgraph(comp).number_of_edges()
        assert edges == len(comp), "each component of a 2-regular graph is a cycle"
        parts.append(edges // 2)
    return tuple(sorted(parts, reverse=True))


def all_perfect_matchings(num_ends: int) -> Iterator[tuple[tuple[int, int], ...]]:
    """Every perfect matching of ``range(num_ends)`` as a tuple of pairs
    (brute force).
    """

    def rec(remaining: tuple[int, ...]):
        if not remaining:
            yield ()
            return
        v, rest = remaining[0], remaining[1:]
        for i, w in enumerate(rest):
            for tail in rec(rest[:i] + rest[i + 1 :]):
                yield ((v, w),) + tail

    yield from rec(tuple(range(num_ends)))


def compositions(n: int) -> Iterator[tuple[int, ...]]:
    """All ordered DSB distributions (b_1, ..., b_k) with every b_i >= 1 and sum n."""
    if n == 0:
        yield ()
        return
    for first in range(1, n + 1):
        for rest in compositions(n - first):
            yield (first,) + rest


def reference_is_deranged(num_dsbs: int, rejoin_pairs) -> bool:
    """No parallel edges in Ξ (a rejoin edge parallel to a DSB edge)."""
    xi = exchange_multigraph(num_dsbs, rejoin_pairs)
    return all(xi.number_of_edges(u, v) == 1 for u, v in xi.edges())


def reference_is_proper(breaks: Sequence[int], rejoin_pairs) -> bool:
    """Paper definition, on explicit multigraphs: connected AMG, no parallel
    edges in Ξ.
    """
    return reference_is_deranged(sum(breaks), rejoin_pairs) and nx.is_connected(
        amg_multigraph(breaks, rejoin_pairs)
    )


def reference_states(
    breaks: Sequence[int], space: str
) -> set[tuple[tuple[int, int], ...]]:
    """Brute-force state set as sorted pair tuples; ``space`` is 'all',
    'deranged' or 'proper'.
    """
    n = sum(breaks)
    keep = {
        "all": lambda pairs: True,
        "deranged": lambda pairs: reference_is_deranged(n, pairs),
        "proper": lambda pairs: reference_is_proper(breaks, pairs),
    }[space]
    return {
        tuple(sorted(pairs)) for pairs in all_perfect_matchings(2 * n) if keep(pairs)
    }
