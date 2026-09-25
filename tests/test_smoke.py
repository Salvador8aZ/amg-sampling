import amg_sampling


def test_package_imports():
    assert amg_sampling.__version__ == "0.1.0.dev0"


def test_core_imports_no_third_party_packages():
    import subprocess
    import sys

    code = (
        "import sys, amg_sampling.core; "
        "bad = {'networkx', 'numpy', 'hydra', 'mlflow', 'pennylane'} & set(sys.modules); "
        "assert not bad, bad"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def test_only_the_application_imports_hydra():
    import subprocess
    import sys

    code = (
        "import sys, amg_sampling.core, amg_sampling.exact, amg_sampling.samplers, "
        "amg_sampling.data.sheth, amg_sampling.analysis.runs; "
        "assert 'hydra' not in sys.modules and 'omegaconf' not in sys.modules"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
