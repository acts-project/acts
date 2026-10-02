"""Exercise uv installs with view-only packages and an inherited native package."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import pytest

SETUP = Path(__file__).with_name("setup.sh")


def run(*args, env=None):
    return subprocess.check_output(args, env=env, text=True).strip()


def site_packages(python):
    return Path(
        run(str(python), "-c", "import sysconfig; print(sysconfig.get_path('purelib'))")
    )


@pytest.mark.parametrize("full_install", ["false", "true"])
def test_packages_visible_in_view_but_not_venv(tmp_path, full_install):
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("uv is required")
    wheels = tmp_path / "wheels"
    wheels.mkdir()
    packages = [
        ("pyyaml", "yaml", "1.0"),
        ("jinja2", "jinja2", "1.0"),
        ("histcmp", "histcmp", "0.10.0"),
        ("matplotlib", "matplotlib", "1.0"),
        ("pytest_md_report", "pytest_md_report", "1.0"),
    ]
    for name, module, version in packages:
        info = f"{name}-{version}.dist-info"
        with zipfile.ZipFile(
            wheels / f"{name}-{version}-py3-none-any.whl", "w"
        ) as wheel:
            wheel.writestr(f"{module}.py", "fixture = True\n")
            wheel.writestr(
                f"{info}/METADATA",
                f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n"
                "Requires-Dist: acts-native-stack\n",
            )
            wheel.writestr(
                f"{info}/WHEEL",
                "Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n",
            )
            wheel.writestr(f"{info}/RECORD", "")
    env = {
        **os.environ,
        "UV_CACHE_DIR": str(tmp_path / "cache"),
        "UV_OFFLINE": "1",
        "UV_NO_INDEX": "1",
        "UV_FIND_LINKS": str(wheels),
    }
    view = tmp_path / "view"
    venv = tmp_path / "venv"
    run(sys.executable, "-m", "venv", "--without-pip", str(view))
    view_python = view / "bin/python3"
    run(
        uv,
        "pip",
        "install",
        "--python",
        str(view_python),
        "--no-deps",
        "pyyaml",
        "jinja2",
        env=env,
    )
    run(str(view_python), "-m", "venv", "--system-site-packages", str(venv))
    python = venv / "bin/python3"
    # A nested venv loses the view's site-packages, reproducing the CI failure.
    if (
        run(
            str(python),
            "-c",
            "import importlib.util; print(importlib.util.find_spec('yaml'))",
        )
        != "None"
    ):
        pytest.skip("base Python already supplies PyYAML")
    inherited = tmp_path / "native"
    inherited.mkdir()
    (inherited / "acts_native_stack.py").write_text("native = True\n")
    metadata = inherited / "acts_native_stack-1.0.dist-info"
    metadata.mkdir()
    (metadata / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: acts-native-stack\nVersion: 1.0\n"
    )
    for interpreter in (view_python, python):
        (site_packages(interpreter) / "native.pth").write_text(str(inherited) + "\n")
    requirements = tmp_path / "Python/Examples/tests/requirements.txt"
    requirements.parent.mkdir(parents=True)
    requirements.write_text("acts-native-stack==1.0\n")
    script_dir = tmp_path / "CI/dependencies"
    script_dir.mkdir(parents=True)
    section = SETUP.read_text().split(
        'start_section "Prepare python environment"\n', 1
    )[1]
    section = section.split('checkpoint "Python environment prepared"', 1)[0]
    run(
        "bash",
        "-euc",
        'retry_transient() { "$@"; }\n' + section,
        env={
            **env,
            "view_dir": str(view),
            "venv_dir": str(venv),
            "SCRIPT_DIR": str(script_dir),
            "full_install": full_install,
        },
    )
    run(
        str(python), "-c", "import yaml, jinja2; assert yaml.fixture and jinja2.fixture"
    )
    assert run(
        str(python), "-c", "import acts_native_stack; print(acts_native_stack.__file__)"
    ) == str(inherited / "acts_native_stack.py")
    assert not list(site_packages(python).glob("acts_native_stack*"))
    if full_install == "true":
        run(
            str(python),
            "-c",
            "import histcmp, matplotlib, pytest_md_report; assert histcmp.fixture and matplotlib.fixture and pytest_md_report.fixture",
        )
