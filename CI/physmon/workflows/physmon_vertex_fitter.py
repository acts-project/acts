#!/usr/bin/env python3
import tempfile
from pathlib import Path
import shutil
import datetime

import acts
from acts.examples.simulation import (
    addParticleGun,
    addParticleSelection,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
)

from acts.examples.reconstruction import (
    addParticleSmearing,
    addVertexFitting,
    VertexFinder,
)

import sys

sys.path.append("/home/frusso/hep/acts/CI/physmon/")
sys.path.append("/home/frusso/hep/acts/Examples/Scripts/Python/")
from physmon_common import makeSetup

u = acts.UnitConstants

setup = makeSetup()


def run_vertex_fitting(mu, events):
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=events, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in setup.decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=42)

        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0),
            # PhiConfig(-u.pi, u.pi-1e-8),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=mu,
            rnd=rnd,
        )
        # Dummy particle selection to put "particles_selected" on the whiteboard
        addParticleSelection(
            s,
            inputParticles="particles",
            outputParticles="particles_selected",
        )
        addParticleSmearing(
            s,
            rnd=rnd,
            inputParticles="particles_selected",
            outputTrackParameters="smeared_params",
        )
        addVertexFitting(
            s,
            setup.field,
            trajectories=None,
            trackParameters="smeared_params",
            associatedParticles="particles_selected",
            outputProtoVertices="proto_vertices",
            outputVertices="fitted_vertices",
            vertexFinder=VertexFinder.Truth,
            outputDirRoot=tp,
        )

        s.run()

        del s

        perf_file = tp / f"performance_vertexing.root"
        assert perf_file.exists(), "Performance file not found"
        shutil.copy(
            perf_file,
            setup.outdir / f"performance_vertex_fitting_mu{mu}.root",
        )


with acts.FpeMonitor():
    mu = 200
    start = datetime.datetime.now()

    events = 5
    run_vertex_fitting(mu, events)

    delta = datetime.datetime.now() - start

    duration = delta.total_seconds() / events

    (setup.outdir / f"performance_vertex_fitting_mu{mu}_time.txt").write_text(
        str(duration)
    )
