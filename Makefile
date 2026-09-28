# Package name, read from pyproject.toml so this Makefile is reusable across
# projects. This works because the distribution name matches the importable
# package directory -- keep them identical (underscores, not hyphens) when
# renaming. Any mismatch shows up immediately as 0% coverage.
PACKAGE := $(shell awk -F'"' '/^name = / {print $$2; exit}' pyproject.toml)

# Minimum coverage percentage required for tests to pass
COVERAGE_FAIL = 50

# Run the test suite
test:
	poetry run pytest

# Format the code using Ruff. Applies safe lint fixes (import sorting,
# pyupgrade rewrites) before formatting, so `make check` passes afterwards.
format:
	poetry run ruff check --fix .
	poetry run ruff format .

# Verify formatting without rewriting anything (this is what CI runs)
format-check:
	poetry run ruff format --check .

# Lint the code using Ruff (configured in pyproject.toml [tool.ruff])
lint:
	poetry run ruff check .

# Run all quality checks. Does not modify files; run `make format` to fix.
check: format-check lint test

# Run tests with coverage enforcement (terminal output only)
# Omit patterns are configured in pyproject.toml [tool.coverage.run].
coverage:
	poetry run coverage run --source=$(PACKAGE) -m pytest
	poetry run coverage report --fail-under=$(COVERAGE_FAIL)

# Run tests with coverage and produce an HTML report
coverage-html:
	poetry run coverage run --source=$(PACKAGE) -m pytest
	poetry run coverage report --fail-under=$(COVERAGE_FAIL)
	poetry run coverage html
	@echo "HTML coverage report generated at htmlcov/index.html"

# Verify every optional dependency is isolated behind a single module.
import-boundaries:
	poetry run python scripts/import_boundaries.py .

# Verify imported packages are declared, and declared in the right group.
deps-check:
	poetry run deptry .

# Report unused code. ADVISORY: a library's public API is uncalled by
# construction, so read the output rather than trusting it.
# deadcode-whitelist.py is picked up automatically once it exists, so a
# baselined project reports only newly dead code.
deadcode:
	poetry run vulture $(PACKAGE) scripts $(wildcard deadcode-whitelist.py)

# Baseline the existing dead code in a project that already carries some, so it
# stops blocking work while new dead code still surfaces. Commit the result.
# Two traps are handled here, both of which turn `make check` red if this is
# done by hand: vulture exits 3 whenever it finds anything, and it emits a
# trailing blank line that `ruff format --check` rejects.
deadcode-baseline:
	-poetry run vulture --make-whitelist $(PACKAGE) scripts > deadcode-whitelist.py
	poetry run ruff format deadcode-whitelist.py
	@if [ -s deadcode-whitelist.py ]; then \
	    echo "baseline written to deadcode-whitelist.py - review it, then commit it"; \
	else \
	    rm -f deadcode-whitelist.py; \
	    echo "no dead code found; no baseline needed"; \
	fi

# Build a wheel, install it into a throwaway virtualenv and run the shipped
# tests against the INSTALLED package. The recipe changes directory first
# because from the repository root `import $(PACKAGE)` resolves to the source
# tree, which would silently test the wrong code.
test-wheel:
	rm -rf dist .wheeltest
	poetry build -f wheel
	python -m venv .wheeltest
	@# The shipped tests need the `plots` extra for the figure tests, and
	@# hypothesis and networkx, which are dev dependencies.
	./.wheeltest/bin/pip install --quiet "$$(ls dist/*.whl)[plots]" pytest hypothesis networkx
	cd "$$(mktemp -d)" && $(CURDIR)/.wheeltest/bin/python -m pytest --pyargs $(PACKAGE).tests -q
	rm -rf .wheeltest
	@echo "shipped tests pass against the installed wheel"

.PHONY: test format format-check lint check coverage coverage-html \
        import-boundaries deps-check deadcode deadcode-baseline test-wheel
