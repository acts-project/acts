from pathlib import Path
import subprocess
import sys

import pytest

from helpers import dd4hepEnabled, geant4Enabled, hepmc3Enabled, rootEnabled


@pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled and dd4hepEnabled),
    reason="Geant4 material recording dependencies not available",
)
def test_persistent_material_recording_matches_beam_on(tmp_path):
    import awkward as ak
    import uproot

    scripts = Path(__file__).resolve().parents[3] / "Examples/Scripts/Python"
    # Each fresh process owns one Geant4 manager. Use a nonzero event offset to
    # check that ACTS-derived seeds still produce identical material records.
    script = """
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import acts
import acts.examples
from material_recording import runMaterialRecording

output = Path(sys.argv[2])
output.mkdir()
detector = acts.examples.TelescopeDetector(
    acts.examples.TelescopeDetector.Config(
        positions=[30, 60, 90], stereos=[0, 0, 0]
    )
)
sequencer = acts.examples.Sequencer(
    events=3, skip=5, numThreads=1, outputDir=str(output)
)
runMaterialRecording(
    detector=detector,
    s=sequencer,
    tracksPerEvent=20,
    outputFileBase=str(output / "material"),
    persistentRun=sys.argv[3] == "true",
).run()
"""
    outputs = []
    for persistent in (False, True):
        output = tmp_path / str(persistent)
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                script,
                str(scripts),
                str(output),
                str(persistent).lower(),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        with uproot.open(output / "material.root") as file:
            tree = file["material_tracks"]
            assert tree.num_entries == 60
            outputs.append(tree.arrays(library="ak"))
    assert ak.almost_equal(outputs[0], outputs[1], rtol=0, atol=0)
