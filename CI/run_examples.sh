#!/bin/bash
# Run a reasonably large subset of the Acts examples

# Force bash to catch more errors
set -euo pipefail

# Common setup
NUM_EVENTS=100
SRC_DIR=`pwd`
BUILD_DIR=`pwd`/build
DD4HEP_INPUT="--dd4hep-input file:${SRC_DIR}/Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml"
run_example() {
    echo ""
    echo "=== Running $* -n ${NUM_EVENTS} ==="
    echo ""
    time ${BUILD_DIR}/bin/$* -n ${NUM_EVENTS}
}

# Run hello world example
run_example ActsExampleHelloWorld

# Run geometry examples for all geometries
#
# We must try to avoid running examples for all geometries because
# that results in a combinatorial explosion of CI running time. But
# these examples are fast enough.
#
run_example ActsExampleGeometryAligned
run_example ActsExampleGeometryDD4hep ${DD4HEP_INPUT}
run_example ActsExampleGeometryEmpty
run_example ActsExampleGeometryGeneric
run_example ActsExampleGeometryPayload
run_example ActsExampleGeometryTelescope
# TODO: Add TGeo geometry example (needs an input file + knowhow)

# Run propagation examples for all geometries
#
# In addition to the geometry examples, we also need one slightly more
# complex example to run with all geometries because geometries with
# conditions have code paths that are only exercised when the
# Sequencer actually runs.
#
run_example ActsExamplePropagationAligned
run_example ActsExamplePropagationDD4hep ${DD4HEP_INPUT}
run_example ActsExampleGeometryEmpty
run_example ActsExampleGeometryGeneric
run_example ActsExampleGeometryPayload
# TODO: Add TGeo propagation example (needs an input file + knowhow)

# Run event generation examples
run_example ActsExampleParticleGun
run_example ActsExamplePythia8

# Run Geantino recording example
#
# FIXME: Currently only works in single-threaded mode, even though it
#        should theoretically be thread-safe as Geant4 usage is
#        protected by a mutex. What most likely happens is that some
#        thread-unsafe Geant4 code is accidentally run outside of the
#        mutex-protected region of the code. See issue #207 .
#
run_example ActsExampleGeantinoRecordingDD4hep ${DD4HEP_INPUT} -j1
# TODO: Add GDML version (needs an input file + knowhow)

# Run material validation example (generic-only, see above)
run_example ActsExampleMaterialValidationGeneric

# TODO: Run CKF examples (needs some setup, documented in Acts howto?)
# TODO: Run EventRecording example (requires CSV input, source to be found. Fatras?)
# TODO: Run Fatras examples (needs some setup, documented in Acts howto?)
# TODO: Run HepMC3 examples (needs input files + knowhow)
# TODO: Run MagneticField examples (needs input files + knowhow)
# TODO: Run MaterialMapping examples (needs input files + knowhow)
# TODO: Run Seeding examples (needs input files + knowhow)
# TODO: Run Show(FatrasGeneric|Particles) examples (needs input files + knowhow)
# TODO: Run TruthTracks examples (needs input files + knowhow)
# TODO: Run (AdaptiveMulti|Iterative)?Vertex(Finder|Fitter) examples (currently broken, not ready for CI)

# TODO: Bring back multi-threaded output reproducibility tests

# Run vertex finder tutorial
run_example ActsTutorialVertexFinder
