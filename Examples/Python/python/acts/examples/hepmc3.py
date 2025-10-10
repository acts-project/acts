"""
Backward compatibility shim for acts.examples.hepmc3.

This module maintains backward compatibility by re-exporting everything
from the hepmc3 package. Users can still do:
    from acts.examples import hepmc3
    hepmc3.normalize(...)
    hepmc3.Compression.zstd
    etc.
"""

from .hepmc3 import *  # noqa: F401, F403
