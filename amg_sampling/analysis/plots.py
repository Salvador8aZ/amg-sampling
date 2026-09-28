"""Optional figures (requires the ``plots`` extra: matplotlib).

``make_plots`` returns the files written, or an empty list with a message when
matplotlib is not installed. Figures read only the results dictionary.
"""

from __future__ import annotations

from pathlib import Path


def matplotlib_available() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def make_plots(results: dict, directory: str | Path) -> list[Path]:
    if not matplotlib_available():
        return []
    import matplotlib

    matplotlib.use("Agg")
    directory = Path(directory)
    written = []
    analysis = results["analysis"]
    if analysis == "distribution":
        written.append(_cycle_count_plot(results, directory / "cycle_count_distribution.png"))
    if analysis == "observed":
        written.append(_observed_structures_plot(results, directory / "observed_cycle_structures.png"))
    if analysis in ("observed", "rejoin_probability") and any("null_probability" in r for r in results["rejoins"]):
        written.append(_rejoin_plot(results, directory / "rejoin_probabilities.png"))
    return [p for p in written if p is not None]


def _title(results: dict) -> str:
    problem = results["problem"]
    who = problem.get("patient_id")
    theta = problem["theta"]["notation"]
    return f"{who}  {theta}" if who else theta


def _cycle_count_plot(results, path):
    import matplotlib.pyplot as plt

    rows = [r for r in results["cycle_count_distribution"] if r.get("sampled_count", 0) or _exact(r) >= 1e-4]
    if not rows:
        return None
    xs = [r["num_cycles"] for r in rows]
    fig, ax = plt.subplots(figsize=(7, 4))
    if any("sampled_probability" in r for r in rows):
        ax.bar(xs, [r.get("sampled_probability", 0) for r in rows], color="#4C72B0", alpha=0.8,
               yerr=[1.96 * r.get("standard_error", 0) for r in rows], capsize=3, label="IID sample (±1.96 s.e.)")
    if any(r.get("exact_probability") for r in rows):
        ax.plot(xs, [_exact(r) for r in rows], "o", color="#C44E52", label="exact")
    ax.set_xlabel("number of exchange cycles")
    ax.set_ylabel("probability under the null")
    ax.set_title(f"{_title(results)}: {results['null_model']}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _observed_structures_plot(results, path):
    import matplotlib.pyplot as plt

    rows = results["observed_cycle_structures_under_null"]
    if not rows:
        return None
    labels = [r["cycle_structure"] for r in rows]
    xs = range(len(rows))
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(rows)), 4))
    width = 0.4
    share = [r["share_among_completions"]["value"] if r.get("share_among_completions") else 0 for r in rows]
    null = [
        r["null_exact"]["value"] if r.get("null_exact") else (r.get("null_sampled") or {}).get("estimate", 0)
        for r in rows
    ]
    if any(share):
        ax.bar([x - width / 2 for x in xs], share, width, color="#55A868", label="share of completions of the observed rejoins")
    ax.bar([x + width / 2 for x in xs], null, width, color="#4C72B0", label=f"null: {results['null_model']}")
    ax.set_xticks(list(xs), labels)
    ax.set_ylabel("probability")
    observed = results.get("observed_cycle_structure")
    ax.set_title(f"{_title(results)}" + (f"  (observed: {observed})" if observed else ""))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _rejoin_plot(results, path):
    import matplotlib.pyplot as plt

    rows = [r for r in results["rejoins"] if "null_probability" in r]
    labels = [f"{r['edge']}: {r['site_a']} — {r['site_b']}" for r in rows]
    values = [r["null_probability"]["estimate"] for r in rows]
    errors = [1.96 * r["null_probability"]["standard_error"] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 0.5 * len(rows) + 1.5))
    ax.barh(range(len(rows)), values, xerr=errors, color="#4C72B0", capsize=3)
    n = results["problem"]["theta"]["num_dsbs"]
    ax.axvline(1 / (2 * n - 1), color="#C44E52", linestyle="--", label="1/(2n−1): any given pair under uniform ALL")
    ax.set_yticks(range(len(rows)), labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"P(edge present) under {results['null_model']}  (±1.96 s.e.)")
    ax.set_title(_title(results))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _exact(row) -> float:
    p = row.get("exact_probability")
    return p["value"] if p else 0.0
