import os
import tempfile

import detray.core
import detray.io
import detray.tests


def test_generate_and_read_detector():
    # Use the toy generator to produce detector files to read back.
    out_dir = tempfile.mkdtemp() + os.sep

    det, det_names = detray.tests.buildToyDetector(
        detray.core.HostMemoryResource(), detray.tests.ToyDetectorConfig()
    )

    writer_config = detray.io.DetectorWriterConfig()
    writer_config.path = out_dir
    detray.io.writeDetector(det, det_names, writer_config)

    geometry_file = os.path.join(out_dir, "toy_detector_geometry.json")
    assert os.path.exists(geometry_file)

    reader_config = detray.io.DetectorReaderConfig().addFile(geometry_file)

    detector, names = detray.io.readDetector(
        detray.core.HostMemoryResource(), reader_config
    )
    assert isinstance(names, detray.core.NameMap)
    assert len(detector.volumes) == 22
    assert len(detector.surfaces) == 3230
