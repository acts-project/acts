#!/usr/bin/env python3
# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from pathlib import Path

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

    # Option 1: Implement a track finding algorithm...
    class PythonTrackFinder(acts.examples.IAlgorithm):
        def __init__(self, name, level):
            acts.examples.IAlgorithm.__init__(self, name, level)

            self.spacepoints = acts.examples.ReadDataHandle(
                self, acts.SpacePointContainer, "Spacepoints"
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

    # ... or option 2: use truth values
    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=acts.logging.INFO,
        inputParticles="particles_generated_selected",
        inputMeasurements="measurements",
        inputParticleMeasurementsMap="particle_measurements_map",
        inputSimHits="simhits",
        inputMeasurementSimHitsMap="measurement_simhits_map",
        outputProtoTracks="prototracks",
    )
    # s.addAlgorithm(truthTrkFndAlg)

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

            self.spacepoints = acts.examples.ReadDataHandle(
                self, acts.SpacePointContainer, "Spacepoints"
            )
            self.spacepoints.initialize("spacepoints")

        def execute(self, context):
            prototracks = self.prototracks(context.eventStore)
            spacepoints = self.spacepoints(context.eventStore)

            measurement_to_spacepoint = {}
            measurement_to_sourcelink = {}
            for sp in spacepoints:
                for sl in sp.sourceLinks:
                    isl = acts.examples.IndexSourceLink.FromSourceLink(sl)
                    meas_id = isl.index()
                    measurement_to_spacepoint[meas_id] = sp
                    measurement_to_sourcelink[meas_id] = sl
            surface_map = trackingGeometry.geoIdSurfaceMap()

            container = acts.examples.TrackContainer()
            for prototrack in prototracks:
                track = container.makeTrack()
                track.parameters = acts.BoundVector(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
                track.nMeasurements = len(prototrack)

                for meas_id in prototrack:
                    sp = measurement_to_spacepoint[meas_id]
                    sl = measurement_to_sourcelink[meas_id]
                    isl = acts.examples.IndexSourceLink.FromSourceLink(sl)
                    sf = surface_map[isl.geometryId()]

                    trackState = track.appendTrackState()
                    trackState.typeFlags.isMeasurement = True
                    trackState.uncalibratedSourceLink = sl
                    trackState.referenceSurface = sf

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

    # Add track finder performance writer
    cfg_finder = acts.examples.PythonTrackFinderPerformanceWriter.Config()
    cfg_finder.inputTracks = "fitted_tracks"
    cfg_finder.inputParticles = "particles"
    cfg_finder.inputTrackParticleMatching = "track_particle_matching"
    cfg_finder.inputParticleTrackMatching = "particle_track_matching"
    cfg_finder.inputParticleMeasurementsMap = "particle_measurements_map"
    perfWriterFinder = acts.examples.PythonTrackFinderPerformanceWriter(
        cfg_finder, acts.logging.INFO
    )
    s.addWriter(perfWriterFinder)

    # Add track fitter performance writer
    cfg_fitter = acts.examples.PythonTrackFitterPerformanceWriter.Config()
    cfg_fitter.inputTracks = "fitted_tracks"
    cfg_fitter.inputParticles = "particles"
    cfg_fitter.inputTrackParticleMatching = "track_particle_matching"
    perfWriterFitter = acts.examples.PythonTrackFitterPerformanceWriter(
        cfg_fitter, acts.logging.INFO
    )
    s.addWriter(perfWriterFitter)

    return s, perfWriterFinder, perfWriterFitter


if __name__ == "__main__":
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector = acts.examples.GenericDetector(acts.examples.GenericDetector.Config())
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))

    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"
    geoSelectionConfigFile = (
        srcdir / "Examples/Configs/generic-pixel-sstrips-lstrips-spacepoints.json"
    )

    outputDir = Path.cwd() / "output_track_finding_python_only"
    outputDir.mkdir(exist_ok=True)

    s, perfWriterFinder, perfWriterFitter = runTrackFindingPythonOnly(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=digiConfigFile,
        geoSelectionConfigFile=geoSelectionConfigFile,
        outputDir=outputDir,
        decorators=decorators,
    )
    s.run()

    histograms_finder = perfWriterFinder.histograms()
    print(
        f"Retrieved {len(histograms_finder)} finder performance histograms: {list(histograms_finder.keys())}"
    )

    histograms_fitter = perfWriterFitter.histograms()
    print(
        f"Retrieved {len(histograms_fitter)} fitter performance histograms: {list(histograms_fitter.keys())}"
    )
