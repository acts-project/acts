#!/usr/bin/env python3
import acts
import acts.examples
from acts.examples.simulation import addPythia8, addTruthJetAlg, TruthJetConfig


outputDir = "./truth_jet_test_output"
u = acts.UnitConstants
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=10, outputDir=outputDir)

addPythia8(
    s,
    hardProcess=["Top:qqbar2ttbar=on"],
    npileup=50,
    vtxGen=acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
    ),
    rnd=rnd,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir,
)

addTruthJetAlg(
    s,
    TruthJetConfig(
        inputTruthParticles="particles_generated",
        outputJets="output_jets",
        jetPtMin=20 * u.GeV,
    ),
    loglevel=acts.logging.INFO,
)

s.run()
