import uproot
import awkward as ak
import numpy as np
from pathlib import Path

import math
import argparse

# Arguments
p = argparse.ArgumentParser()
p.add_argument("--input", type=Path, default=Path("pythia8_output.root/particles.root"))
p.add_argument("--tree", type=str, default="particles")
args = p.parse_args()

ufile = uproot.open(args.input)
utree = ufile["particles"]

# events = utree.arrays()

events = utree.arrays(library="ak")
print(f"Number of events: {len(events)}")

for px, py, pz, particle_type in zip(
    events["px"], events["py"], events["pz"], events["particle_type"]
):
    # filter
    mask = particle_type == 22  # PDG ID 22 corresponds to photons
    photons_px = px[mask]
    photons_py = py[mask]
    photons_pz = pz[mask]
    photons = ak.zip({"px": photons_px, "py": photons_py, "pz": photons_pz})
    # Count how many are bith E > 20
    photons_selected = photons[
        np.sqrt(photons["px"] ** 2 + photons["py"] ** 2 + photons["pz"] ** 2) > 60
    ]
    print(f"Number of photons in event: {len(photons)}")
    print(f"Number of selected photons in event: {len(photons_selected)}")

# print(f"Number of photons: {len(photons)}")
