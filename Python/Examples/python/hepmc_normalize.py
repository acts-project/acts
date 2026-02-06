#!/usr/bin/env python3
"""Standalone script for HepMC3 file normalization.

This is a convenience wrapper that calls the hepmc3 module's CLI.
It's equivalent to running: python -m acts.examples.hepmc3
"""

import sys
import os

# Import and run the CLI from the hepmc3 package
from acts.examples.hepmc3.__main__ import main

if __name__ == "__main__":
    sys.exit(main(prog=os.path.basename(sys.argv[0])))
