"""ddsim steering file for a material recording with geantinos.

Shoots geantinos through the Geant4 geometry that DDG4 builds from a DD4hep
detector and records the material they cross with the `ActsMaterialRecorder`
DDG4 stepping action. The recorded material tracks are written to
`--outputFile` by the `ActsMaterialTrackWriter` event action, in the format the
ACTS material mapping reads, instead of the usual simulation output.

The `ActsPluginDD4hepGeant4` library has to be on `LD_LIBRARY_PATH`, e.g.
by sourcing `this_acts.sh`. Usage:

    ddsim --steeringFile material_recording_ddsim.py \\
          --compactFile detector.xml \\
          --numberOfEvents 1000 \\
          --outputFile geant4_material_tracks.root

Gun settings can be overridden on the command line, e.g. `--gun.etaMin`,
`--gun.etaMax`, `--gun.multiplicity` or `--random.seed`. Multi-threaded
running with `--numberOfThreads`, available from DD4hep 1.38, writes a single
output file.
"""

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import GeV

SIM = DD4hepSimulation()
SIM.runType = "batch"
SIM.outputFile = "geant4_material_tracks.root"

SIM.enableGun = True
SIM.gun.particle = "geantino"
SIM.gun.multiplicity = 1000
SIM.gun.distribution = "eta"
SIM.gun.etaMin = -4.0
SIM.gun.etaMax = 4.0
SIM.gun.momentumMin = 1 * GeV
SIM.gun.momentumMax = 1 * GeV
SIM.gun.position = (0.0, 0.0, 0.0)


def materialRecordingOutput(dd4hepSimulation):
    """Record material tracks to the output file, instead of the simulation output"""
    from DDG4 import EventAction, Kernel, SteppingAction

    kernel = Kernel().worker()

    # One recorder per worker thread, collecting the material tracks per event
    recorder = SteppingAction(kernel, "ActsMaterialRecorder/MaterialRecorder")
    recorder.ExcludeMaterials = ["Air", "Vacuum"]
    recorder.RecordElementFractions = False
    kernel.steppingAction().add(recorder)

    # One writer, shared by all worker threads in multi-threaded running. The
    # kernel knows the number of threads independently of the ddsim version.
    shared = Kernel().NumberOfThreads > 1
    writer = EventAction(kernel, "ActsMaterialTrackWriter/MaterialTrackWriter", shared)
    writer.OutputFile = dd4hepSimulation.outputFile
    writer.TreeName = "material_tracks"
    kernel.eventAction().add(writer)


SIM.outputConfig.userOutputPlugin = materialRecordingOutput
