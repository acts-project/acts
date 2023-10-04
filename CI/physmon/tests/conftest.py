from pathlib import Path
import re
from typing import List, IO, Tuple, Optional
import threading
import csv
from datetime import datetime
import time
import shutil

import pytest
import psutil
from pytest_check import check

import acts
import acts.examples
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


@pytest.fixture(scope="session")
def output_path(request):
    path: Path = request.config.getoption("--physmon-output-path").resolve()
    path.mkdir(parents=True, exist_ok=True)

    return path


@pytest.fixture(scope="session")
def reference_path(request):
    path: Path = request.config.getoption("--physmon-reference-path").resolve()

    return path


class Physmon:
    detector: "acts.examples.dd4hep.DD4hepDetector"
    trackingGeometry: acts.TrackingGeometry
    decorators: List[acts.IMaterialDecorator]
    field: acts.MagneticFieldProvider
    digiConfig: Path
    geoSel: Path
    output_path: Path
    reference_path: Path
    update_references: bool
    tmp_path: Path
    name: str

    def __init__(
        self,
        detector,
        trackingGeometry,
        decorators,
        field,
        digiConfig,
        geoSel,
        output_path,
        reference_path,
        update_references,
        tmp_path,
        name,
    ):
        self.detector = detector
        self.trackingGeometry = trackingGeometry
        self.decorators = decorators
        self.field = field
        self.digiConfig = digiConfig
        self.geoSel = geoSel
        self.output_path = output_path
        self.reference_path = reference_path
        self.update_references = update_references
        self.tmp_path = tmp_path
        self.name = name

        self.test_output_path.mkdir(exist_ok=True, parents=True)
        self.test_reference_path.mkdir(exist_ok=True, parents=True)

    @property
    def test_output_path(self) -> Path:
        return self.output_path / self.name

    @property
    def test_reference_path(self) -> Path:
        return self.reference_path / self.name

    def add_output_file(self, filename: str, rename: Optional[str] = None):
        __tracebackhide__ = True
        tmp = self.tmp_path / filename
        assert tmp.exists(), f"Output file {tmp} does not exist"
        outname = rename if rename else filename
        shutil.copy(tmp, self.test_output_path / outname)

    def histogram_comparison(
        self, filename: Path, title: str, config_path: Optional[Path] = None
    ):
        __tracebackhide__ = True
        monitored = self.test_output_path / filename
        reference = self.test_reference_path / filename

        assert monitored.exists(), f"Output file {monitored} does not exist"

        if self.update_references:
            shutil.copy(monitored, reference)
        assert reference.exists(), f"Reference file {reference} does not exist"

        from histcmp.console import Console
        from histcmp.report import make_report
        from histcmp.checks import Status
        from histcmp.config import Config
        from histcmp.github import is_github_actions, github_actions_marker
        from histcmp.cli import print_summary

        from histcmp.compare import compare, Comparison

        from rich.panel import Panel
        from rich.console import Group
        from rich.pretty import Pretty

        import yaml

        console = Console()

        console.print(
            Panel(
                Group(f"Monitored: {monitored}", f"Reference: {reference}"),
                title="Comparing files:",
            )
        )

        if config_path is None:
            config = Config.default()
        else:
            with config_path.open() as fh:
                config = Config(**yaml.safe_load(fh))

        console.print(Panel(Pretty(config), title="Configuration"))

        #  filter_path = Path(_filter)
        #  if filter_path.exists():
        #  with filter_path.open() as fh:
        #  filters = fh.read().strip().split("\n")
        #  else:
        #  filters = [_filter]
        filters = []
        comparison = compare(
            config, monitored, reference, filters=filters, console=console
        )

        comparison.label_monitored = "monitored"
        comparison.label_reference = "reference"
        comparison.title = title

        status = print_summary(comparison, console)

        plots = self.test_output_path / "plots"
        plots.mkdir(exist_ok=True, parents=True)
        report_file = self.test_output_path / f"{monitored.stem}.html"
        make_report(comparison, report_file, console, plots, format="pdf")

        for item in comparison.items:
            msg = f"{item.key} failures: " + ", ".join(
                [c.name for c in item.checks if c.status == Status.FAILURE]
            )
            check.equal(item.status, Status.SUCCESS, msg=msg)


@pytest.fixture(scope="session")
def _physmon_prereqs():
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    matDeco = acts.IMaterialDecorator.fromFile(
        srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
        level=acts.logging.INFO,
    )

    detector, trackingGeometry, decorators = getOpenDataDetector(
        getOpenDataDetectorDirectory(), matDeco
    )

    return srcdir, detector, trackingGeometry, decorators


def _sanitize_test_name(name: str):
    name = re.sub(r"^test_", "", name)
    name = re.sub(r"\]$", "", name)
    name = re.sub(r"[\[\]]", "_", name)
    return name


@pytest.fixture()
def physmon(
    output_path: Path, reference_path: Path, _physmon_prereqs, tmp_path, request
):
    u = acts.UnitConstants

    srcdir, detector, trackingGeometry, decorators = _physmon_prereqs

    setup = Physmon(
        detector=detector,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        digiConfig=srcdir
        / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        geoSel=srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json",
        field=acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T)),
        output_path=output_path,
        reference_path=reference_path,
        update_references=request.config.getoption("--physmon-update-references"),
        tmp_path=tmp_path,
        name=_sanitize_test_name(request.node.name),
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
def monitor(physmon: Physmon, request, capsys):
    import psutil

    p = psutil.Process()

    memory = physmon.output_path / "memory"
    memory.mkdir(exist_ok=True)
    name = _sanitize_test_name(request.node.name)
    with (memory / f"mem_{name}.csv").open("w") as fh:
        mon = Monitor(output=fh, interval=0.1)

        t = threading.Thread(target=mon.run, args=(p,))
        t.start()

        yield

        mon.terminate = True
        t.join()

        # @TODO: Add plotting

        #  with capsys.disabled():
        #  print("MONITORING")
