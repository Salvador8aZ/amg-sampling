"""Human-readable terminal report built from an analysis results dictionary."""

from __future__ import annotations

from amg_sampling.analysis.estimates import format_probability as fp

SPACE_DESCRIPTIONS = {
    "all": "ALL (every perfect matching of the free ends)",
    "deranged": "DERANGED (no DSB repaired to itself, i.e. no C1)",
    "proper": "PROPER (deranged and connected: the paper's proper AMGs)",
}
RULE = "─" * 72


def format_report(results: dict, top: int = 12) -> str:
    lines: list[str] = []
    add = lines.append
    problem = results["problem"]
    kind = "Patient" if problem["kind"] == "patient" else "Configuration"
    add(f"AMG Sampling — {kind} analysis: {results['analysis']}")
    add(RULE)
    _problem_block(problem, add)
    add("")
    add(f"State space:      {SPACE_DESCRIPTIONS[results['state_space']]}")
    add(f"Null model:       {results['null_model']}  (uniform combinatorial null)")
    sampling = results.get("sampling") or (results.get("null") or {}).get("sampling")
    _sampling_block(results, sampling, add)
    exact = results.get("exact") or (results.get("null") or {}).get("exact")
    _exact_block(exact, results["state_space"], add)

    if results["analysis"] == "distribution":
        _distribution_sections(results, top, add)
    if results["analysis"] == "observed":
        _observed_sections(results, add)
    if results["analysis"] in ("observed", "rejoin_probability"):
        _rejoin_sections(results, add)
    add("")
    add("Interpretation")
    for chunk in _wrap(results["interpretation"], 70):
        add(f"  {chunk}")
    return "\n".join(lines)


def _problem_block(problem, add):
    theta = problem["theta"]
    if problem["kind"] == "patient":
        add(f"Patient:          {problem['patient_id']}")
        chroms = problem["chromosomes"]
        add(
            f"Chromosomes:      {', '.join(map(str, chroms))}  "
            f"({theta['num_chromosomes']} chromosomes)"
        )
        alloc = ", ".join(f"chr{c}: {b}" for c, b in zip(chroms, theta["breaks"]))
        source = problem.get("data_source")
        if source:
            tag = (
                "reference copy"
                if source["is_reference_copy"]
                else "NOT the reference copy"
            )
            add(
                f"Data:             {source['file_name']} "
                f"({tag}, sha256 {source['sha256'][:12]}…)"
            )
    else:
        add(f"Configuration:    {problem['label']}")
        add(f"Chromosomes:      {theta['num_chromosomes']}")
        alloc = ", ".join(map(str, theta["breaks"]))
    add(f"Initial config.:  {theta['notation']}")
    add(f"DSB allocation:   {alloc}")
    add(f"Total DSBs:       {theta['num_dsbs']}")
    add(f"Free DSB ends:    {theta['num_free_ends']}")


def _sampling_block(results, sampling, add):
    sampler = results["sampler"]
    if sampler["name"] not in ("iid", "mcmc") or sampling is None:
        add(f"Sampling method:  exact only (no sampling); seed {sampler['seed']}")
        return
    if sampler["name"] == "mcmc":
        add(
            "Sampling method:  MCMC, Metropolis–Hastings with 2-switch (reversal) moves"
        )
        add(f"Samples:          {sampling['samples']:,}")
        add(
            f"Burn-in, thin:    {sampler['burn_in']:,} steps, "
            f"then every {sampler['thin']:,} steps"
        )
        add(f"Moves proposed:   {sampling['proposals']:,}")
        add(f"Move acceptance:  {fp(sampling['acceptance_rate'])}")
        add(
            f"Effective size:   {sampling['effective_samples_num_cycles']:,.0f} "
            "(number of cycles; autocorrelation time "
            f"{sampling['integrated_autocorrelation_time']:.3g})"
        )
        add(f"Seed:             {sampler['seed']}")
        add("                  standard errors use batch means (correlated samples)")
        return
    theory = sampling.get("theoretical_acceptance")
    theory_text = f" (exact {fp(theory['value'])})" if theory else ""
    add(
        "Sampling method:  IID rejection sampling "
        "(uniform proposals, exact acceptance test)"
    )
    add(f"Samples:          {sampling['samples']:,}")
    add(f"Proposals:        {sampling['proposals']:,}")
    add(f"Acceptance rate:  {fp(sampling['acceptance_rate'])}{theory_text}")
    add(f"Seed:             {sampler['seed']}")


def _exact_block(exact, space, add):
    if exact is None:
        return
    total = exact.get("total_states")
    if total is not None:
        add(
            f"|{space.upper()}(Θ)|:{' ' * max(1, 12 - len(space))}{_big(total)} (exact)"
        )
    if exact.get("cycle_distribution_available"):
        add("Exact reference:  full cycle-structure distribution available")
    else:
        reasons = "; ".join(exact.get("unavailable", {}).values())
        add(f"Exact reference:  cycle-structure distribution unavailable ({reasons})")
    if exact.get("total_variation_bound"):
        add(
            "                  scalar summaries use exact Uniform(DERANGED); "
            "it differs from "
            f"Uniform(PROPER) by at most "
            f"{fp(exact['total_variation_bound']['value'])} in any probability"
        )


def _distribution_sections(results, top, add):
    rows = results["cycle_distribution"]
    add("")
    add(f"Cycle structures (top {min(top, len(rows))} of {len(rows)} observed/exact)")
    add(_table(rows[:top], "cycle_structure", "structure"))
    means = results.get("means", {})
    rows = [
        r
        for r in results["cycle_count_distribution"]
        if r.get("sampled_count", 0) > 0 or _exact_value(r) >= 1e-4
    ]
    add("")
    add("Number of cycles |C| (values with a sample or exact probability ≥ 1e-4)")
    add(_table(rows, "num_cycles", "|C|"))
    add(_mean_line(means.get("num_cycles")))
    rows = sorted(
        results["largest_cycle_distribution"],
        key=lambda r: -(_exact_value(r) or r.get("sampled_probability", 0.0)),
    )[:5]
    add("")
    add(
        "Largest cycle C_l (five most likely values; "
        "full table in largest_cycle_distribution.csv)"
    )
    add(_table(rows, "largest_cycle", "l"))
    add(_mean_line(means.get("largest_cycle")))


def _mean_line(mean) -> str:
    if not mean:
        return ""
    parts = []
    if mean.get("sampled_mean") is not None:
        parts.append(
            f"sampled {mean['sampled_mean']:.4g} ± {mean['sampled_standard_error']:.2g}"
        )
    if mean.get("exact_mean") is not None:
        parts.append(f"exact {mean['exact_mean']:.4g}")
    return f"  mean: {'; '.join(parts)}"


def _observed_sections(results, add):
    obs = results["observed"]
    add("")
    add("Observed data")
    add(f"  Observed rejoins (junctions):  {obs['num_observed_rejoins']}")
    add(
        f"  Unmatched DSB ends:            {obs['num_unmatched_ends']}"
        f"  →  {obs['num_completions_all']:,} possible completions"
    )
    comp = obs["completions"]
    if comp.get("exact_count") is not None:
        add(
            f"  Completions in {comp['space'].upper()}:         "
            f"{comp['exact_count']:,} (exact enumeration)"
        )
    elif comp.get("unavailable"):
        add(
            f"  Completions in {comp['space'].upper()}:         "
            f"not enumerated ({comp['unavailable']})"
        )
    rec = obs.get("reconstructed")
    if rec:
        add("")
        add(
            "Reconstructed configuration "
            "(the observed rejoins leave a unique completion)"
        )
        for edge in rec["rejoins"]:
            add(f"  {edge['site_a']} — {edge['site_b']}   [{edge['status']}]")
        member = ", ".join(
            f"{k.upper()}: {'yes' if v else 'no'}" for k, v in rec["membership"].items()
        )
        add(f"  Cycle structure: {rec['cycle_structure']}    ({member})")
    rows = results["observed_cycle_structures_under_null"]
    add("")
    if rec:
        add("Observed cycle structure under the null")
    else:
        add(
            "Cycle structures of AMGs consistent with the observed rejoins, "
            "and their null probability"
        )
    add(
        "  Share of completions: among AMGs that contain every observed rejoin "
        "(exact; sampled by"
    )
    add("  conditional completion). Null: probability under the uniform null for Θ.")
    header = (
        f"  {'structure':<10} {'share (exact)':>19} {'share (sampled)':>20} "
        f"{'null (exact)':>13} {'null (sampled)':>24}"
    )
    add(header)
    for r in rows:
        share = r.get("share_among_completions")
        share_text = (
            f"{fp(share['value'])} ({share['fraction']})"
            if share and share.get("fraction")
            else "—"
        )
        exact = r.get("null_exact")
        add(
            f"  {r['cycle_structure']:<10} {share_text:>19} "
            f"{r.get('share_among_completions_sampled_text', '—'):>20} "
            f"{fp(exact['value']) if exact else '—':>13} "
            f"{r.get('null_sampled_text', '—'):>24}"
        )


def _sequence_features(row) -> str:
    """``; homology 11 bp; insertion 0 bp`` for junctions whose features are known."""
    parts = []
    for key, name in (
        ("homology_length", "homology"),
        ("foreign_sequence_length", "insertion"),
    ):
        if row.get(key) is not None:
            parts.append(f"{name} {row[key]} bp")
    return "".join(f"; {p}" for p in parts)


def _rejoin_sections(results, add):
    add("")
    add("Observed rejoin edges: probability of each edge under the null")
    for r in results["rejoins"]:
        text = r.get("null_probability_text", "not sampled")
        add(f"  {r['edge']:<4} {r['site_a']:>26} — {r['site_b']:<26} {text}")
        add(
            f"       [{r['status']}; source line {r['source_line']}; "
            f"class {r['sv_class']}{_sequence_features(r)}]"
        )
    joint = results.get("joint_observed_rejoins")
    if joint:
        add("")
        add(
            f"All {joint['num_observed_rejoins']} observed rejoins "
            "present simultaneously"
        )
        if joint.get("exact"):
            add(
                f"  exact:   {fp(joint['exact']['value'])}  "
                f"({_big(joint['exact_numerator'])} / "
                f"{_big(joint['exact_denominator'])})"
            )
        if joint.get("sampled_text"):
            add(f"  sampled: {joint['sampled_text']}")


def _table(rows, key, label):
    has_exact = any("exact_probability" in r for r in rows)
    has_sampled = any("sampled_probability" in r for r in rows)
    head = f"  {label:<14}"
    if has_sampled:
        head += f" {'sampled':>10} {'± s.e.':>10}"
    if has_exact:
        head += f" {'exact':>12}"
    out = [head]
    for r in rows:
        line = f"  {str(r[key]):<14}"
        if has_sampled:
            line += (
                f" {fp(r.get('sampled_probability')):>10} "
                f"{fp(r.get('standard_error')):>10}"
            )
        if has_exact:
            line += f" {fp(_exact_value(r)) if r.get('exact_probability') else '—':>12}"
        out.append(line)
    return "\n".join(out)


def _exact_value(row) -> float:
    p = row.get("exact_probability")
    return p["value"] if p else 0.0


def _big(value) -> str:
    try:
        number = int(value)
    except (TypeError, ValueError):
        return str(value)
    if number < 10**15:
        return f"{number:,}"
    digits = str(number)
    return f"{digits[0]}.{digits[1:5]}e+{len(digits) - 1}"


def _wrap(text, width):
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        out.append(line)
    return out
