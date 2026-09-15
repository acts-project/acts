#!/usr/bin/env python3
"""Check the wheel wrapper's container environment without building wheels."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class WheelEnvironmentTests(unittest.TestCase):
    def test_compiler_cache_budget_is_forwarded(self):
        with tempfile.TemporaryDirectory() as directory:
            stub = Path(directory) / "uv"
            stub.write_text(
                f"#!{sys.executable}\n"
                "import json, os\n"
                "print(json.dumps({k: os.environ.get(k) for k in "
                "['CCACHE_MAXSIZE', 'CIBW_ENVIRONMENT_PASS', 'CIBW_ENVIRONMENT_LINUX']}))\n"
            )
            stub.chmod(0o755)
            env = os.environ | {
                "PATH": directory + os.pathsep + os.environ["PATH"],
                "CIBW_BUILD": "cp314-manylinux_x86_64",
                "CCACHE_MAXSIZE": "750M",
                "CCACHE_DIR": "/tmp/wheel-cache",
            }
            result = subprocess.check_output(
                ["bash", str(Path(__file__).with_name("cibuildwheel.sh"))],
                env=env,
                text=True,
            )
            configured = json.loads(result)
            self.assertEqual(configured["CCACHE_MAXSIZE"], "750M")
            self.assertIn("CCACHE_MAXSIZE", configured["CIBW_ENVIRONMENT_PASS"].split())
            self.assertIn(
                "CCACHE_DIR=/host/tmp/wheel-cache", configured["CIBW_ENVIRONMENT_LINUX"]
            )


if __name__ == "__main__":
    unittest.main()
