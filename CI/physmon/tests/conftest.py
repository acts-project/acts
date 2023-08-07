from pathlib import Path
import collections
from typing import List, IO, Tuple
import threading
import csv
from datetime import datetime
import time

import pytest
import psutil

import acts
import acts.examples
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


@pytest.fixture()
def output_path(request):
    path: Path = request.config.getoption("--physmon-output-path").resolve()
    path.mkdir(parents=True, exist_ok=True)

    return path


PhysmonSetup = collections.namedtuple(
    "Setup",
    [
        "detector",
        "trackingGeometry",
        "decorators",
        "field",
        "digiConfig",
        "geoSel",
    ],
)


@pytest.fixture(scope="session")
def setup():
    u = acts.UnitConstants
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    matDeco = acts.IMaterialDecorator.fromFile(
        srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
        level=acts.logging.INFO,
    )

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(), matDeco
    )
    setup = PhysmonSetup(
        detector=detector,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        digiConfig=srcdir
        / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        geoSel=srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json",
        field=acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T)),
    )
    return setup


# @TODO: try using spyral directly


class Monitor:
    interval: float

    max_rss: float = 0
    max_vms: float = 0

    terminate: bool = False

    exception = None

    def __init__(
        self,
        output: IO[str],
        interval: float = 0.5,
    ):
        self.interval = interval
        self.writer = csv.writer(output)
        self.writer.writerow(("time", "rss", "vms"))


        self.time: List[float] = [0]
        self.rss: List[float] = [0]
        self.vms: List[float] = [0]

    @staticmethod
    def _get_memory(p: psutil.Process) -> Tuple[float, float]:
        rss = p.memory_info().rss
        vms = p.memory_info().vms
        for subp in p.children(recursive=True):
            try:
                rss += subp.memory_info().rss
                vms += subp.memory_info().vms
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        return rss, vms

    def run(self, p: psutil.Process):

        rss_offset, vms_offset = self._get_memory(p)
        self.writer.writerow((0, 0, 0))

        try:
            start = datetime.now()
            while p.is_running() and p.status() in (
                psutil.STATUS_RUNNING,
                psutil.STATUS_SLEEPING,
            ):
                if self.terminate:
                    return

                delta = (datetime.now() - start).total_seconds()

                rss, vms = self._get_memory(p)
                rss -= rss_offset
                vms -= vms_offset

                self.rss.append(rss / 1e6)
                self.vms.append(vms / 1e6)
                self.time.append(delta)
                self.max_rss = max(rss, self.max_rss)
                self.max_vms = max(vms, self.max_vms)

                self.writer.writerow((delta, rss, vms))

                time.sleep(self.interval)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return
        except Exception as e:
            self.exception = e
            raise e


@pytest.fixture(autouse=True)
def monitor(output_path: Path, request, capsys):
    import psutil

    p = psutil.Process()

    memory = output_path / "memory"
    memory.mkdir(exist_ok=True)
    with (memory / f"mem_{request.node.name}.csv").open("w") as fh:
        mon = Monitor(output=fh, interval=0.1)

        t = threading.Thread(target=mon.run, args=(p,))
        t.start()

        yield

        mon.terminate = True
        t.join()

        # @TODO: Add plotting

        with capsys.disabled():
            print("MONITORING")


