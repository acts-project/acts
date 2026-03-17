#!/usr/bin/env python3

import argparse

import acts
from acts import (
    MaterialValidater,
    IntersectionMaterialAssigner,
)

from acts.examples.simulation import (
    addParticleGun,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
)

from acts.examples import (
    MaterialValidation,
    ParticleTrackParamExtractor,
)

from acts.examples.odd import getOpenDataDetector
from acts.examples.root import RootMaterialTrackWriter

u = acts.UnitConstants


def runMaterialValidation(
    surfaces,
    s,
    tracksPerEvent=100,
    etaRange=(-4.0, 4.0),
    phiRange=(0.0, 360.0 * u.degree),
    materialTrackCollectionName="material_tracks",
    outputFile="material_validation",
    trackingGeometry=None,
):

    rnd = acts.examples.RandomNumbers(seed=42)

    addParticleGun(
        s,
        ParticleConfig(
            num=tracksPerEvent, pdg=acts.PdgParticle.eMuon, randomizeCharge=True
        ),
        EtaConfig(*etaRange),
        PhiConfig(*phiRange),
        rnd=rnd,
    )

    trkParamExtractor = ParticleTrackParamExtractor(
        level=acts.logging.INFO,
        inputParticles="particles_generated",
        outputTrackParameters="params_particles_generated",
    )
    s.addAlgorithm(trkParamExtractor)

    # Assignment setup : Intersection assigner
    materialAssingerConfig = IntersectionMaterialAssigner.Config()
    materialAssingerConfig.surfaces = surfaces
    materialAssinger = IntersectionMaterialAssigner(
        materialAssingerConfig, acts.logging.INFO
    )

    # Validater setup
    materialValidaterConfig = MaterialValidater.Config()
    materialValidaterConfig.materialAssigner = materialAssinger
    materialValidater = MaterialValidater(materialValidaterConfig, acts.logging.INFO)

    # Validation Algorithm
    materialValidationConfig = MaterialValidation.Config()
    materialValidationConfig.inputTrackParameters = "params_particles_generated"
    materialValidationConfig.materialValidater = materialValidater
    materialValidationConfig.outputMaterialTracks = materialTrackCollectionName
    materialValidation = MaterialValidation(materialValidationConfig, acts.logging.INFO)
    s.addAlgorithm(materialValidation)

    # Add the mapped material tracks writer
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=materialValidationConfig.outputMaterialTracks,
            filePath=str(outputFile) + ".root",
            storeSurface=True,
            storeVolume=True,
        )
    )

    # This is configured to
    if trackingGeometry is not None:
        nav = acts.Navigator(trackingGeometry=trackingGeometry)
        stepper = acts.StraightLineStepper()
        propagator = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

        propagationAlgorithm = acts.examples.PropagationAlgorithm(
            propagatorImpl=propagator,
            level=acts.logging.INFO,
            sterileLogger=False,
            inputTrackParameters="params_particles_generated",
            outputSummaryCollection="propagation_summary",
            outputMaterialCollection=materialTrackCollectionName + "_propagated",
        )
        s.addAlgorithm(propagationAlgorithm)

        # Add the mapped material tracks writer
        s.addWriter(
            RootMaterialTrackWriter(
                level=acts.logging.INFO,
                inputMaterialTracks=materialTrackCollectionName + "_propagated",
                filePath=str(outputFile) + "_propagated.root",
                storeSurface=True,
                storeVolume=True,
            )
        )

    return s


def main():
    p = argparse.ArgumentParser()

    p.add_argument(
        "-n", "--events", type=int, default=1000, help="Number of events to process"
    )
    p.add_argument(
        "-t", "--tracks", type=int, default=100, help="Number of tracks per event"
    )
    p.add_argument("-j", "--threads", type=int, default=-1, help="Number of threads")
    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )
    p.add_argument(
        "--eta-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(-4.0, 4.0),
        help="Eta range for generated particles",
    )
    p.add_argument(
        "--phi-range",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(0.0, 360.0),
        help="Phi range in degree for generated particles",
    )
    p.add_argument(
        "--material-track-collection",
        type=str,
        default="material_tracks",
        help="Output material track collection name",
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="propagation-material",
        help="Output file stem (without extension)",
    )
    p.add_argument(
        "-p",
        "--propagate",
        action="store_true",
        help="Enable propagation validation",
    )

    args = p.parse_args()
    materialDecorator = None
    if args.map != "":
        materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    detector = getOpenDataDetector(materialDecorator)
    trackingGeometry = detector.trackingGeometry()

    materialSurfaces = trackingGeometry.extractMaterialSurfaces()

    s = acts.examples.Sequencer(events=args.events, numThreads=args.threads)

    runMaterialValidation(
        surfaces=materialSurfaces,
        s=s,
        tracksPerEvent=args.tracks,
        etaRange=tuple(args.eta_range),
        phiRange=(args.phi_range[0] * u.degree, args.phi_range[1] * u.degree),
        materialTrackCollectionName=args.material_track_collection,
        outputFile=args.output,
        trackingGeometry=trackingGeometry if args.propagate else None,
    ).run()


if "__main__" == __name__:
    main()
