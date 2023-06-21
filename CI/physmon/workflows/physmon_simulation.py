#!/usr/bin/env python3

import tempfile
from pathlib import Path
import shutil

import acts
from acts.examples.simulation import (
    addFatras,
    addGeant4,
)

from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()


with tempfile.TemporaryDirectory() as temp:
    tp = Path(temp)

    rnd = acts.examples.RandomNumbers(seed=42)

    s = acts.examples.Sequencer(
        events=10000,
        numThreads=-1,
        logLevel=acts.logging.INFO,
        fpeMasks=acts.examples.Sequencer.FpeMask.fromFile(
            Path(__file__).parent.parent / "fpe_masks.yml"
        ),
    )

    for d in setup.decorators:
        s.addContextDecorator(d)

    s.addReader(
        acts.exampels.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.exampels.EventGenerator.Generator(
                    multiplicity=acts.exampels.FixedMultiplicityGenerator(n=1),
                    vertex=acts.examples.GaussianVertexGenerator(
                        mean=acts.Vector4(0, 0, 0, 0),
                        stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 1 * u.ns),
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 100 * u.GeV),
                        pTransverse=True,
                        eta=(-3.0, 3.0),
                        phi=(0.0 * u.degree, 360.0 * u.degree),
                        pdg=pdg,
                        randomizeCharge=True,
                    ),
                )
                for pdg in [
                    acts.PdgParticle.eMuon,
                    acts.PdgParticle.ePionPlus,
                    acts.PdgParticle.eElectron,
                ]
            ],
            outputParticles="particles_input",
            randomNumbers=rnd,
        )
    )

    addFatras(
        s,
        setup.trackingGeometry,
        setup.field,
        rnd,
        detector=setup.detector,
        enableInteractions=True,
    )

    addGeant4(
        s,
        setup.trackingGeometry,
        setup.field,
        rnd,
        detector=setup.detector,
    )

    s.run()
    del s

    perf_file = tp / "performance_gsf.root"
    assert perf_file.exists(), "Performance file not found"
    shutil.copy(perf_file, setup.outdir / "performance_gsf.root")
