import os
import tempfile

import detray.core
import detray.examples


def test_generate_and_read_detector():
    # Use the toy generator to produce detector files to read back.
    out_dir = tempfile.mkdtemp() + os.sep
    detray.examples.generate_toy_detector(output_dir=out_dir)

    geometry_file = os.path.join(out_dir, "toy_detector_geometry.json")
    assert os.path.exists(geometry_file)

    detector, names = detray.core.read_detector(geometry_file)
    assert isinstance(detector, detray.core.Detector)
    assert isinstance(names, detray.core.NameMap)
