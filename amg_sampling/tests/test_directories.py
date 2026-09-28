from pathlib import Path

import pytest

from amg_sampling import directories


def _is_source_checkout():
    return directories.base("pyproject.toml").exists()


class TestDirectories:
    def test_directories_exist(self):
        assert directories.base().is_dir()
        assert directories.code().is_dir()
        assert directories.tests().is_dir()
        assert directories.test_data().is_dir()

    def test_filenames(self):
        assert directories.code("__init__.py").exists()
        assert directories.tests("__init__.py").exists()
        assert directories.test_data("README.md").exists()

    def test_returns_path_objects(self):
        assert isinstance(directories.code(), Path)
        assert isinstance(directories.code("__init__.py"), Path)

    def test_package_data_resolves_inside_the_package(self):
        assert directories.package_data() == directories.code("data")
        assert directories.package_data("x.csv") == directories.code("data") / "x.csv"

    def test_data_resolves_at_the_repository_root(self):
        assert directories.data() == directories.base("data")

    @pytest.mark.skipif(
        not _is_source_checkout(),
        reason="base() and data() address the repository, not an installed package",
    )
    def test_source_checkout_layout(self):
        assert directories.base("README.md").exists()
        assert directories.data().is_dir()
