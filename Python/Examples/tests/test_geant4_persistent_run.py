from pathlib import Path
import subprocess
import sys

import pytest

from helpers import dd4hepEnabled, geant4Enabled, hepmc3Enabled, rootEnabled


@pytest.mark.skipif(
    not (geant4Enabled and hepmc3Enabled and rootEnabled and dd4hepEnabled),
    reason="Geant4 material recording dependencies not available",
)
def test_material_recording_event_seeds_independent_of_run_length(tmp_path):
    import awkward as ak
    import uproot

    scripts = Path(__file__).resolve().parents[3] / "Examples/Scripts/Python"
    # Compare a multi-event run with the same events in separate processes.
    # A nonzero offset checks that seeds depend on the ACTS event number rather
    # than the event's position in the Geant4 run.
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
    events=int(sys.argv[3]), skip=int(sys.argv[4]), numThreads=1,
    outputDir=str(output)
)
runMaterialRecording(
    detector=detector,
    s=sequencer,
    tracksPerEvent=20,
    outputFileBase=str(output / "material"),
).run()
"""

    def record(name, events, skip):
        output = tmp_path / name
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                script,
                str(scripts),
                str(output),
                str(events),
                str(skip),
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        with uproot.open(output / "material.root") as file:
            tree = file["material_tracks"]
            assert tree.num_entries == events * 20
            return tree.arrays(library="ak")

    combined = record("combined", 3, 5)
    separate = ak.concatenate([record(str(event), 1, event) for event in range(5, 8)])
    assert ak.almost_equal(combined, separate, rtol=0, atol=0)
