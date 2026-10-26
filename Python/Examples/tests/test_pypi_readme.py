from importlib.metadata import PackageNotFoundError, metadata
from pathlib import Path

import pytest

REFERENCE = Path(__file__).resolve().parents[3] / "CI" / "pypi_readme.md"


@pytest.mark.pypi
def test_pypi_readme():
    """The readme in the wheel metadata, i.e. the PyPI landing page, matches the
    checked-in CI/pypi_readme.md, so changes to the page show up in review."""
    try:
        readme = metadata("pyacts").json["description"]
    except PackageNotFoundError:
        pytest.skip("pyacts wheel not installed")

    expected = REFERENCE.read_text().rstrip("\n")
    assert (
        readme.rstrip("\n") == expected
    ), "PyPI readme changed, run `uvx hatch-fancy-pypi-readme -o CI/pypi_readme.md`"
