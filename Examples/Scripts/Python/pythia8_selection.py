#!/usr/bin/env python3
from pathlib import Path

import argparse
from typing import Optional

import acts
import acts.examples
from acts.examples.simulation import addPythia8
from acts.examples.pythia8 import Pythia8Event, Pythia8Particle


u = acts.UnitConstants

# Arguments
p = argparse.ArgumentParser()
p.add_argument("--events", "-n", type=int, default=100)
p.add_argument("--threads", "-j", type=int, default=-1)
p.add_argument("--processes", type=str, nargs="+", help="Pythia8 hard processes")
p.add_argument("--nhard", type=int, default=1, help="Number of hard scatterings")
p.add_argument("--seed", type=int, default=42, help="Random seed")
p.add_argument("--output-hepmc3", type=Path, default=Path("pythia8_output.hepmc3"))
p.add_argument("--output-root", type=Path, default=Path("pythia8_output.root"))
p.add_argument("--output-csv", type=Path, default=Path("pythia8_output.csv"))
p.add_argument("--npileup", type=int, default=0, help="Number of pileup interactions")
args = p.parse_args()

# Preliminaries
rnd = acts.examples.RandomNumbers(seed=args.seed)

# Sequencer
s = acts.examples.Sequencer(
    events=args.events, numThreads=args.threads, logLevel=acts.logging.INFO
)

# Higgs physics example:
# Generate Higgs bosons via gluon fusion (process 102) with Pythia
# and force them to decay to Z0 Z0 (PDG ID 23) which then decay to
# a pair of top quarks (PDG ID 6).
# This can be achieved by setting the following Pythia8 commands:
# 102:gg2H = on # to switch on gluon fusion Higgs production
# 25:m0 = 125.0 # set Higgs mass to 125 GeV
#
# H -> Z0 Z0
# hardProcess = ["HiggsSM:gg2H = on", "25:m0 = 125.0", "25:onMode = off", "25:onIfMatch = 23 23"  ],


def select_higgs_to_zz(
    event: Pythia8Event,
    final_state: str = "4l",
    eta_sel: Optional[tuple[float, float]] = (-1, 1),
) -> bool:
    # loop over particles in the event and look for a Higgs boson decaying to two Z bosons
    for p in event.particles():
        if p.id() == 25:  # Higgs boson PDG ID
            daughters = [event.particles()[dau_idx] for dau_idx in p.daughterList()]
            if len(daughters) != 2:
                continue
            if all(d.id() == 23 for d in daughters):  # Both daughters are Z0 bosons
                # Now check the decay products of the Z bosons
                lepton_count = 0
                for z in daughters:
                    z_daughters = [
                        event.particles()[dau_idx] for dau_idx in z.daughterList()
                    ]
                    for zd in z_daughters:
                        if abs(zd.id()) in [11, 13]:  # Electron or muon PDG IDs
                            if final_state == "4l":
                                lepton_count += 1
                            if eta_sel and (
                                zd.eta() < eta_sel[0] or zd.eta() > eta_sel[1]
                            ):
                                return False
                if final_state == "4l" and lepton_count == 4:
                    return True

    return False


#
# H -> gamma gamma
# hardProcess = ["HiggsSM:gg2H = on", "25:m0 = 125.0", "25:onMode = off", "25:onIfMatch = 22 22"  ],


addPythia8(
    s,
    rnd=rnd,
    nhard=args.nhard,
    npileup=args.npileup,
    hardProcess=[
        "HiggsSM:gg2H = on",
        "25:m0 = 125.0",
        "25:onMode = off",
        "25:onIfMatch = 23 23",
        "23:onMode = off",
        "23:onIfMatch = 11 -11",
        "23:onIfMatch = 13 -13",
    ],
    # hardProcess = ["HiggsSM:gg2H = on", "25:m0 = 125.0", "25:onMode = off", "25:onIfMatch = 22 22"  ],
    pileupProcess=["SoftQCD:all = on"],
    outputDirRoot=args.output_root,
    outputDirCsv=args.output_csv,
    writeHepMC3=args.output_hepmc3,
    hardProcessSelectors=[select_higgs_to_zz],
)

s.run()
