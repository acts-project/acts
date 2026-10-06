"""Check the repair hook selects uv's patchelf ahead of the image's version."""

import os
from pathlib import Path
import subprocess


def test_repair_uses_uv_tool_bin(tmp_path):
    bin_dir = tmp_path / "bin"
    tool_bin = tmp_path / "uv tools"
    bin_dir.mkdir()
    tool_bin.mkdir()
    scripts = {
        bin_dir / "uv": """#!/bin/sh
case "$2" in
  run) printf '%s' "$CIBW_REPAIR_WHEEL_COMMAND_LINUX" ;;
  install) exit 0 ;;
  dir) printf '%s' "$TEST_TOOL_BIN" ;;
  *) exit 1 ;;
esac
""",
        bin_dir / "patchelf": "#!/bin/sh\necho old\n",
        tool_bin / "patchelf": "#!/bin/sh\necho pinned\n",
        bin_dir / "auditwheel": '#!/bin/sh\n[ "$(patchelf)" = pinned ]\n',
    }
    for path, content in scripts.items():
        path.write_text(content)
        path.chmod(0o755)
    env = {
        **os.environ,
        "PATH": f"{bin_dir}:{os.environ['PATH']}",
        "TEST_TOOL_BIN": str(tool_bin),
        "CIBW_BUILD": "cp313-manylinux_x86_64",
    }
    command = subprocess.check_output(
        ["bash", str(Path(__file__).with_name("cibuildwheel.sh"))], env=env, text=True
    )
    subprocess.run(["sh", "-ec", command], env=env, check=True)
