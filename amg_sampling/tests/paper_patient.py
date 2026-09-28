"""Structure of patient P05-1657 (paper section 7.1, figure 13, table 2).

Derived by ``amg_sampling.data.sheth.convert`` from the reference copy of the
dataset (SHA-256 13a70c63…) and recorded here as labels only (DSB counts per
chromosome and observed rejoin end pairs), with no positions or other source
fields, so the regression tests run without the data file.
The DSB counts per chromosome (4: 2, 7: 4, 8: 9, 12: 1, 21: 2) match the paper.
"""

P05_1657 = {
    (8, 12): {
        "breaks": (9, 1),
        "observed": ((1, 12), (2, 10), (4, 18), (7, 9), (14, 17)),
    },
    (7,): {"breaks": (4,), "observed": ((0, 6), (3, 4))},
    (4,): {"breaks": (2,), "observed": ((1, 2),)},
    (21,): {"breaks": (2,), "observed": ((1, 2),)},
}

# Paper, table 2: cycle structures of the 945 completions for chromosomes 8 and 12.
TABLE_2 = {
    "C10": 384,
    "C8+C2": 240,
    "C6+C4": 160,
    "C6+2C2": 80,
    "2C4+C2": 60,
    "C4+3C2": 20,
    "5C2": 1,
}
