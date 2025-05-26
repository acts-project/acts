from pathlib import Path
import argparse

import acts
import acts.examples
from acts.examples.reconstruction import addTrackSelection, TrackSelectorConfig

u = acts.UnitConstants


def runGNN4ITk(
    inputRootDump: Path,
    moduleMapPath: str,
    gnnModel: Path,
    events: int = 1,
    logLevel=acts.logging.INFO,
):
    assert inputRootDump.exists()
    assert Path(moduleMapPath + ".doublets.root").exists()
    assert Path(moduleMapPath + ".triplets.root").exists()
    assert gnnModel.exists()

    moduleMapConfig = {
        "level": logLevel,
        "moduleMapPath": moduleMapPath,
        "rScale": 1000.0,
        "phiScale": 3.141592654,
        "zScale": 1000.0,
        "gpuDevice": 0,
        "gpuBlocks": 512,
        "moreParallel": True,
    }

    gnnConfig = {
        "level": logLevel,
        "cut": 0.5,
        "modelPath": str(gnnModel),
        "useEdgeFeatures": True,
    }

    builderCfg = {
        "level": logLevel,
        "useOneBlockImplementation": False,
        "doJunctionRemoval": True,
    }

    graphConstructor = acts.examples.ModuleMapCuda(**moduleMapConfig)
    if gnnModel.suffix == ".pt":
        edgeClassifier = acts.examples.TorchEdgeClassifier(**gnnConfig)
    elif gnnModel.suffix == ".onnx":
        del gnnConfig["useEdgeFeatures"]
        edgeClassifier = acts.examples.OnnxEdgeClassifier(**gnnConfig)
    elif gnnModel.suffix == ".engine":
        edgeClassifier = acts.examples.TensorRTEdgeClassifier(**gnnConfig)
    trackBuilder = acts.examples.CudaTrackBuilding(**builderCfg)

    s = acts.examples.Sequencer(
        events=events,
        numThreads=1,
    )

    s.addReader(
        acts.examples.RootAthenaDumpReader(
            level=logLevel,
            treename="GNN4ITk",
            inputfiles=[str(inputRootDump)],
            outputSpacePoints="spacepoints",
            outputClusters="clusters",
            outputMeasurements="measurements",
            outputMeasurementParticlesMap="measurement_particles_map",
            outputParticleMeasurementsMap="particle_measurements_map",
            outputParticles="particles",
            skipOverlapSPsPhi=True,
            skipOverlapSPsEta=False,
            absBoundaryTolerance=0.01 * u.mm,
        )
    )

    e = acts.examples.NodeFeature
    s.addAlgorithm(
        acts.examples.TrackFindingAlgorithmExaTrkX(
            level=logLevel,
            graphConstructor=graphConstructor,
            edgeClassifiers=[edgeClassifier],
            trackBuilder=trackBuilder,
            nodeFeatures=[
                e.R,
                e.Phi,
                e.Z,
                e.Eta,
                e.Cluster1R,
                e.Cluster1Phi,
                e.Cluster1Z,
                e.Cluster1Eta,
                e.Cluster2R,
                e.Cluster2Phi,
                e.Cluster2Z,
                e.Cluster2Eta,
            ],
            featureScales=[1000.0, 3.14159265359, 1000.0, 1.0] * 3,
            inputSpacePoints="spacepoints",
            inputClusters="clusters",
            outputProtoTracks="prototracks",
        )
    )

    s.addAlgorithm(
        acts.examples.PrototracksToTracks(
            level=logLevel,
            inputProtoTracks="prototracks",
            inputMeasurements="measurements",
            outputTracks="gnn_only_tracks",
        )
    )

    s.addAlgorithm(
        acts.examples.ParticleSelector(
            level=logLevel,
            ptMin=1 * u.GeV,
            rhoMax=26 * u.cm,
            measurementsMin=7,
            removeSecondaries=True,
            removeNeutral=True,
            excludeAbsPdgs=[
                11,
            ],
            inputParticles="particles",
            outputParticles="particles_selected",
            inputParticleMeasurementsMap="particle_measurements_map",
        )
    )

    addTrackSelection(
        s,
        TrackSelectorConfig(nMeasurementsMin=7, requireReferenceSurface=False),
        inputTracks="gnn_only_tracks",
        outputTracks="gnn_only_tracks_selected",
        logLevel=logLevel,
    )

    # NOTE: This is not standard ATLAS matching
    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=logLevel,
            inputTracks="gnn_only_tracks_selected",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="tpm",
            outputParticleTrackMatching="ptm",
            doubleMatching=True,
        )
    )

    s.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=logLevel,
            inputParticles="particles_selected",
            inputParticleMeasurementsMap="particle_measurements_map",
            inputTrackParticleMatching="tpm",
            inputParticleTrackMatching="ptm",
            inputTracks="gnn_only_tracks_selected",
            filePath="performance_gnn4itk.root",
        )
    )

    s.run()


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Run the GNN4ITk example")

    argparser.add_argument(
        "--inputRootDump",
        type=Path,
        required=True,
        help="Path to the input ROOT dump file",
    )
    argparser.add_argument(
        "--moduleMapPath",
        type=str,
        required=True,
        help="Path to the module map file (without .doublets.root or .triplets.root suffixes)",
    )
    argparser.add_argument(
        "--gnnModel",
        type=Path,
        required=True,
        help="Path to the GNN model file (ONNX or Torch)",
    )
    argparser.add_argument(
        "--events",
        type=int,
        default=1,
        help="Number of events to process",
    )

    args = argparser.parse_args()

    runGNN4ITk(
        inputRootDump=args.inputRootDump,
        moduleMapPath=args.moduleMapPath,
        gnnModel=args.gnnModel,
        events=args.events,
    )
