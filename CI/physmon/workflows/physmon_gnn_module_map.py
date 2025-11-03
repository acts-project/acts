#!/usr/bin/env python3

"""
Physmon workflow for GNN module map on ODD.

Tests GNN track finding with module map (geometry-based) graph construction
on the OpenDataDetector with simulated data.
"""

from pathlib import Path
import sys

# Add physmon dir to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from physmon_common import makeSetup

import acts
from acts import UnitConstants as u
from acts.examples import Sequencer
from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import addGnn


def run(events=1000):
    """Run GNN module map physmon workflow"""
    setup = makeSetup()

    s = Sequencer(events=events, numThreads=1, outputDir=str(setup.outdir))

    # Random number generator
    rnd = acts.examples.RandomNumbers(seed=42)

    # Particle gun: muons, 1-10 GeV, eta -3 to 3
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0, uniform=True),
        ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
        logLevel=acts.logging.INFO,
    )

    # FATRAS simulation
    addFatras(
        s,
        setup.trackingGeometry,
        setup.field,
        rnd=rnd,
        logLevel=acts.logging.INFO,
    )

    # Digitization
    addDigitization(
        s,
        setup.trackingGeometry,
        setup.field,
        digiConfigFile=setup.digiConfig,
        rnd=rnd,
        logLevel=acts.logging.INFO,
    )

    # Configure GNN stages for module map workflow
    # All parameters hardcoded

    # Stage 1: Graph construction via module map
    graphConstructor = acts.examples.ModuleMapCuda(
        level=acts.logging.INFO,
        moduleMapPath="module_map_odd_2k_events",
        rScale=1000.0,
        phiScale=3.141592654,
        zScale=1000.0,
        etaScale=1.0,
        gpuDevice=0,
        gpuBlocks=512,
        moreParallel=True,
    )

    # Stage 2: Single-stage edge classification (Torch)
    edgeClassifiers = [
        acts.examples.TorchEdgeClassifier(
            level=acts.logging.INFO,
            modelPath="gnn_odd.pt",
            cut=0.5,
            useEdgeFeatures=True,
        )
    ]

    # Stage 3: GPU track building
    trackBuilder = acts.examples.CudaTrackBuilding(
        level=acts.logging.INFO,
        useOneBlockImplementation=False,
        doJunctionRemoval=True,
    )

    # Node features: ITk 12-feature configuration
    e = acts.examples.NodeFeature
    nodeFeatures = [
        e.R, e.Phi, e.Z, e.Eta,
        e.Cluster1R, e.Cluster1Phi, e.Cluster1Z, e.Cluster1Eta,
        e.Cluster2R, e.Cluster2Phi, e.Cluster2Z, e.Cluster2Eta,
    ]
    featureScales = [1000.0, 3.141592654, 1000.0, 1.0] * 3

    # Add GNN tracking
    addGnn(
        s,
        graphConstructor=graphConstructor,
        edgeClassifiers=edgeClassifiers,
        trackBuilder=trackBuilder,
        nodeFeatures=nodeFeatures,
        featureScales=featureScales,
        trackingGeometry=setup.trackingGeometry,
        geometrySelection=setup.geoSel,
        inputClusters="clusters",
        outputDirRoot=str(setup.outdir),
        logLevel=acts.logging.INFO,
    )

    s.run()


if __name__ == "__main__":
    run(events=1000)
