#!/usr/bin/env python3
# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os
from pathlib import Path

os.environ["ACTS_SEQUENCER_DISABLE_FPEMON"] = "1"

import acts
import acts.examples
from acts import UnitConstants as u


def runTrackFindingPythonOnly(
    trackingGeometry,
    field,
    digiConfigFile,
    geoSelectionConfigFile,
    outputDir,
    decorators=[],
    s=None,
):
    from acts.examples.simulation import (
        addParticleGun,
        MomentumConfig,
        EtaConfig,
        PhiConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
    )

    s = s or acts.examples.Sequencer(events=1, numThreads=1, logLevel=acts.logging.INFO)
    outputDir = Path(outputDir)
    rnd = acts.examples.RandomNumbers(seed=42)

    for d in decorators:
        s.addContextDecorator(d)

    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-2.0, 2.0, uniform=True),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    s.addAlgorithm(
        acts.examples.SpacePointMaker(
            level=acts.logging.INFO,
            trackingGeometry=trackingGeometry,
            inputMeasurements="measurement_subset",
            outputSpacePoints="spacepoints",
            geometrySelection=acts.examples.json.readJsonGeometryList(
                str(geoSelectionConfigFile)
            ),
        )
    )

    class PythonTrackFinder(acts.examples.IAlgorithm):
        def __init__(self, name, level):
            acts.examples.IAlgorithm.__init__(self, name, level)

            self.spacepoints = acts.examples.ReadDataHandle(
                self, acts.SpacePointContainer2, "Spacepoints"
            )
            self.spacepoints.initialize("spacepoints")

            self.prototracks = acts.examples.WriteDataHandle(
                self, acts.examples.ProtoTrackContainer, "Prototracks"
            )
            self.prototracks.initialize("prototracks")

        def execute(self, context):
            spacepoints = self.spacepoints(context.eventStore)

            track = acts.examples.ProtoTrack()
            for sp in sorted(spacepoints, key=lambda sp: sp.r):
                for sl in sp.sourceLinks:
                    isl = acts.examples.IndexSourceLink.FromSourceLink(sl)
                    track.append(isl.index())

            prototracks = acts.examples.ProtoTrackContainer()
            prototracks.append(track)

            self.prototracks(context, prototracks)
            return acts.examples.ProcessCode.SUCCESS

    s.addAlgorithm(PythonTrackFinder("PythonTrackFinder", acts.logging.INFO))

    class PythonTrackFitter(acts.examples.IAlgorithm):
        def __init__(self, name, level):
            acts.examples.IAlgorithm.__init__(self, name, level)

            self.prototracks = acts.examples.ReadDataHandle(
                self, acts.examples.ProtoTrackContainer, "Prototracks"
            )
            self.prototracks.initialize("prototracks")

            self.tracks = acts.examples.WriteDataHandle(
                self, acts.examples.ConstTrackContainer, "Tracks"
            )
            self.tracks.initialize("fitted_tracks")

        def execute(self, context):
            prototracks = self.prototracks(context.eventStore)

            container = acts.examples.TrackContainer()
            for prototrack in prototracks:
                track = container.makeTrack()
                track.parameters = acts.BoundVector(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
                track.nMeasurements = len(prototrack)

            self.tracks(context, container.makeConst())
            return acts.examples.ProcessCode.SUCCESS

    s.addAlgorithm(PythonTrackFitter("PythonTrackFitter", acts.logging.INFO))

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="fitted_tracks",
            inputParticles="particles",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="track_particle_matching",
            outputParticleTrackMatching="particle_track_matching",
            doubleMatching=True,
        )
    )

    cfg = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg.inputTracks = "fitted_tracks"
    cfg.inputParticles = "particles"
    cfg.inputTrackParticleMatching = "track_particle_matching"
    cfg.inputParticleTrackMatching = "particle_track_matching"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    perfWriter = acts.examples.PythonTrackFinderPerformanceWriter(
        cfg, acts.logging.INFO
    )
    s.addWriter(perfWriter)

    return s, perfWriter


if __name__ == "__main__":
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector = acts.examples.GenericDetector(acts.examples.GenericDetector.Config())
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))

    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"
    geoSelectionConfigFile = srcdir / "Examples/Configs/generic-seeding-config.json"

    outputDir = Path.cwd() / "output_track_finding_python_only"
    outputDir.mkdir(exist_ok=True)

    s, perfWriter = runTrackFindingPythonOnly(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        geoSelectionConfigFile=geoSelectionConfigFile,
        outputDir=outputDir,
        decorators=decorators,
    )
    s.run()

    histograms = perfWriter.histograms()
    print(
        f"Retrieved {len(histograms)} performance histograms: {list(histograms.keys())}"
    )
