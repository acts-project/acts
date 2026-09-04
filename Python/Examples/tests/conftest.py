import contextlib
import multiprocessing
from pathlib import Path
import sys
import os
import tempfile
import shutil
from typing import Dict
import pytest_check as check
from collections import namedtuple
import filelock

sys.path = [
    str(Path(__file__).parent.parent.parent.parent / "Examples/Scripts/Python/"),
    str(Path(__file__).parent),
] + sys.path


import helpers

import pytest

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import addParticleGun, EtaConfig, ParticleConfig

try:
    import ROOT

    ROOT.gSystem.ResetSignals()
except ImportError:
    pass


u = acts.UnitConstants


class RootHashAssertionError(AssertionError):
    def __init__(
        self, file: Path, key: str, exp_hash: str, act_hash: str, *args, **kwargs
    ):
        super().__init__(f"{exp_hash} != {act_hash}", *args, **kwargs)
        self.file = file
        self.key = key
        self.exp_hash = exp_hash
        self.act_hash = act_hash


hash_assertion_failures = []


def _parse_hash_file(file: Path) -> Dict[str, str]:
    res = {}
    with file.open() as f:
        for line in f:
            if line.strip() == "" or line.strip().startswith("#"):
                continue
            key, h = line.strip().split(":", 1)
            res[key.strip()] = h.strip()
    return res


@pytest.fixture(scope="session")
def root_file_exp_hashes():
    path = Path(
        os.environ.get("ROOT_HASH_FILE", Path(__file__).parent / "root_file_hashes.txt")
    )
    return _parse_hash_file(path)


@pytest.fixture(name="assert_root_hash")
def assert_root_hash(request, root_file_exp_hashes):
    if not helpers.doHashChecks:

        def fn(*args, **kwargs):
            pass

        return fn

    def fn(key: str, file: Path):
        """
        Assertion helper function to check the hashes of root files.
        Do NOT use this function directly by importing, rather use it as a pytest fixture

        Arguments you need to provide:
        key: Explicit lookup key for the expected hash, should be unique per test function
        file: Root file to check the expected hash against
        """
        __tracebackhide__ = True
        # Lazy import: pulls uproot which pulls pyarrow; only relevant for
        # tests that actually use this fixture.
        import helpers.hash_root

        gkey = f"{request.node.name}__{key}"
        act_hash = helpers.hash_root.hash_root_file(file)
        if not gkey in root_file_exp_hashes:
            check.equal(
                act_hash,
                "[MISSING]",
                msg=f'Hash lookup key "{key}" not found for test "{request.node.name}"',
            )
            exc = RootHashAssertionError(file, gkey, "[MISSING]", act_hash)
            hash_assertion_failures.append(exc)

        else:
            refhash = root_file_exp_hashes[gkey]
            check.equal(act_hash, refhash)
            if act_hash != refhash:
                exc = RootHashAssertionError(file, gkey, refhash, act_hash)
                hash_assertion_failures.append(exc)

    return fn


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    docs_url = "https://acts.readthedocs.io/en/latest/examples/python_bindings.html#root-file-hash-regression-checks"
    if len(hash_assertion_failures) > 0:
        terminalreporter.ensure_newline()
        terminalreporter.section(
            "RootHashAssertionErrors", sep="-", red=True, bold=True
        )
        terminalreporter.line(
            "The ROOT files produced by tests have changed since the last recorded reference."
        )
        terminalreporter.line(
            "This can be be expected if e.g. the underlying algorithm changed, or it can be a test failure symptom."
        )
        terminalreporter.line(
            "Please manually check the output files listed below and make sure that their content is correct."
        )
        terminalreporter.line(
            "If it is, you can update the test reference file Python/Examples/tests/root_file_hashes.txt with the new hashes below."
        )
        terminalreporter.line(f"See {docs_url} for more details")
        terminalreporter.line("")

        for e in hash_assertion_failures:
            terminalreporter.line(f"{e.key}: {e.act_hash}")

    if not helpers.doHashChecks:
        terminalreporter.section("Root file hash checks", sep="-", blue=True, bold=True)
        terminalreporter.line(
            "NOTE: Root file hash checks were skipped, enable with ROOT_HASH_CHECKS=on"
        )
        terminalreporter.line(f"See {docs_url} for more details")


def kwargsConstructor(cls, *args, **kwargs):
    return cls(*args, **kwargs)


def configKwConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)
    return cls(*args, config=config, **_kwargs)


def configPosConstructor(cls, *args, **kwargs):
    assert hasattr(cls, "Config")
    _kwargs = {}
    if "level" in kwargs:
        _kwargs["level"] = kwargs.pop("level")
    config = cls.Config()
    for k, v in kwargs.items():
        setattr(config, k, v)

    return cls(config, *args, **_kwargs)


@pytest.fixture(params=[configPosConstructor, configKwConstructor, kwargsConstructor])
def conf_const(request):
    return request.param


@pytest.fixture
def rng():
    return acts.examples.RandomNumbers(seed=42)


@pytest.fixture
def basic_prop_seq(rng):
    def _basic_prop_seq_factory(geo, s=None):
        if s is None:
            s = acts.examples.Sequencer(events=10, numThreads=1)

        addParticleGun(
            s,
            ParticleConfig(num=10, pdg=acts.PdgParticle.eMuon, randomizeCharge=True),
            EtaConfig(-4.0, 4.0),
            rnd=rng,
        )

        trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
            level=acts.logging.WARNING,
            inputParticles="particles_generated",
            outputTrackParameters="params_particles_generated",
        )
        s.addAlgorithm(trkParamExtractor)

        nav = acts.Navigator(trackingGeometry=geo)
        stepper = acts.StraightLineStepper()

        prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

        alg = acts.examples.PropagationAlgorithm(
            level=acts.logging.WARNING,
            propagatorImpl=prop,
            sterileLogger=False,
            inputTrackParameters="params_particles_generated",
            outputSummaryCollection="propagation_summary",
        )
        s.addAlgorithm(alg)

        return s, alg

    return _basic_prop_seq_factory


@pytest.fixture
def trk_geo():
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    yield trackingGeometry


DetectorConfig = namedtuple(
    "DetectorConfig",
    [
        # Tests use this only to scope teardown (``with cfg.detector:``). The
        # detectors below are shared for the whole session, so nothing should
        # release them -- this is a no-op context manager and the real detector
        # is kept alive in _detector_config_cache.
        "detector",
        "trackingGeometry",
        "decorators",
        "geometrySelection",
        "digiConfigFile",
        "name",
    ],
)


_shared_detectors = {}


def _shared_detector(key: str, factory):
    """Build a detector once and reuse it for the rest of the session.

    Constructing the ODD takes ~1s and permanently retains ~16MB: the Acts
    detector elements are installed as extensions on the dd4hep DetElements
    (see DD4hepLayerBuilder::createSensitiveSurface) and so outlive the
    TrackingGeometry they were built for. The suite asks for a detector dozens
    of times, which is where a large part of the pytest process footprint comes
    from.

    Sharing is safe: TrackingGeometry is handed out as
    ``shared_ptr<const TrackingGeometry>`` and neither GenericDetector nor the
    DD4hep detectors populate any context decorators. Note this also means one
    dd4hep detector still backs exactly one TrackingGeometry, which the
    DetElement extension mechanism requires -- it has a single slot per element
    and throws on a second conversion.
    """
    if key not in _shared_detectors:
        _shared_detectors[key] = factory()
    return _shared_detectors[key]


def _get_generic_detector_config(srcdir: Path) -> DetectorConfig:
    detector = _shared_detector("generic", acts.examples.GenericDetector)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    return DetectorConfig(
        contextlib.nullcontext(),
        trackingGeometry,
        decorators,
        geometrySelection=(srcdir / "Examples/Configs/generic-seeding-config.json"),
        digiConfigFile=(srcdir / "Examples/Configs/generic-digi-smearing-config.json"),
        name="generic",
    )


def _get_odd_detector_config(srcdir: Path) -> DetectorConfig:
    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep not set up")

    detector = _shared_detector("odd", getOpenDataDetector)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
    return DetectorConfig(
        contextlib.nullcontext(),
        trackingGeometry,
        decorators,
        digiConfigFile=(srcdir / "Examples/Configs/odd-digi-smearing-config.json"),
        geometrySelection=(srcdir / "Examples/Configs/odd-seeding-config.json"),
        name="odd",
    )


def _srcdir() -> Path:
    return Path(__file__).resolve().parent.parent.parent.parent


@pytest.fixture
def odd_detector():
    """The shared session-lifetime default-configuration ODD (requires DD4hep).

    Prefer this over calling getOpenDataDetector() in a test. Do not use it as
    a context manager -- it is shared, so nothing may release it.
    """
    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep not set up")
    return _shared_detector("odd", getOpenDataDetector)


@pytest.fixture
def odd_detector_gen3():
    """The shared session-lifetime default-configuration Gen3 ODD."""
    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep not set up")
    return _shared_detector("odd-gen3", lambda: getOpenDataDetector(gen3=True))


@pytest.fixture
def generic_detector_config():
    """Detector config for the generic detector only."""
    return _get_generic_detector_config(_srcdir())


@pytest.fixture
def odd_detector_config():
    """Detector config for the Open Data Detector only (requires DD4hep)."""
    return _get_odd_detector_config(_srcdir())


@pytest.fixture(params=["generic", pytest.param("odd", marks=pytest.mark.odd)])
def detector_config(request):
    """Parametrized fixture that runs tests with both generic and ODD detectors."""
    srcdir = _srcdir()
    if request.param == "generic":
        return _get_generic_detector_config(srcdir)
    elif request.param == "odd":
        return _get_odd_detector_config(srcdir)
    else:
        raise ValueError(f"Invalid detector {request.param}")


@pytest.fixture
def ptcl_gun(rng):
    def _factory(s):
        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=2,
                    ),
                )
            ],
            outputEvent="particle_gun_event",
            randomNumbers=rng,
        )

        s.addReader(evGen)

        hepmc3Converter = acts.examples.hepmc3.HepMC3InputConverter(
            level=acts.logging.INFO,
            inputEvent=evGen.config.outputEvent,
            outputParticles="particles_generated",
            outputVertices="vertices_input",
        )
        s.addAlgorithm(hepmc3Converter)

        return evGen, hepmc3Converter

    return _factory


@pytest.fixture
def fatras(ptcl_gun, trk_geo, rng):
    def _factory(s):
        evGen, h3conv = ptcl_gun(s)

        field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
        simAlg = acts.examples.FatrasSimulation(
            level=acts.logging.INFO,
            inputParticles=h3conv.config.outputParticles,
            outputParticles="particles_simulated",
            outputSimHits="simhits",
            randomNumbers=rng,
            trackingGeometry=trk_geo,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=False,
            emEnergyLossIonisation=False,
            emEnergyLossRadiation=False,
            emPhotonConversion=False,
        )

        s.addAlgorithm(simAlg)

        # Digitization
        from acts.examples import json

        digiCfg = acts.examples.DigitizationAlgorithm.Config(
            digitizationConfigs=acts.examples.json.readDigiConfigFromJson(
                str(
                    Path(__file__).parent.parent.parent.parent
                    / "Examples/Configs/generic-digi-smearing-config.json"
                )
            ),
            surfaceByIdentifier=trk_geo.geoIdSurfaceMap(),
            randomNumbers=rng,
            inputSimHits=simAlg.config.outputSimHits,
        )
        digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)

        s.addAlgorithm(digiAlg)

        return evGen, simAlg, digiAlg

    return _factory


def _do_material_recording(d: Path):
    from material_recording import runMaterialRecording

    s = acts.examples.Sequencer(events=2, numThreads=1)

    with getOpenDataDetector() as detector:
        runMaterialRecording(
            detector=detector,
            s=s,
            tracksPerEvent=1000,
            outputFileBase=d / "geant4_material_tracks",
        )
        s.run()


@pytest.fixture(scope="session")
def material_recording_session(tmp_path_factory):
    tmp_dir = tmp_path_factory.getbasetemp().parent
    d = Path(tmp_dir) / "material_recording"

    if not helpers.geant4Enabled:
        pytest.skip("Geantino recording requested, but Geant4 is not set up")

    if not helpers.dd4hepEnabled:
        pytest.skip("DD4hep recording requested, but DD4hep is not set up")

    with filelock.FileLock(str(d) + ".lock"):

        if not d.exists():
            # Record into a scratch directory and only move it into place once
            # the child has succeeded. `d` lives above the per-run basetemp and
            # so outlives the session: publishing it up front means a failed
            # recording leaves an empty directory that every later run happily
            # reuses, turning a setup error into confusing assertion failures.
            staging = d.with_name(d.name + f".incomplete.{os.getpid()}")
            shutil.rmtree(staging, ignore_errors=True)
            staging.mkdir(parents=True)

            # explicitly ask for "spawn" as CI failures were observed with "fork"
            spawn_context = multiprocessing.get_context("spawn")
            p = spawn_context.Process(target=_do_material_recording, args=(staging,))
            p.start()
            p.join()
            if p.exitcode != 0:
                shutil.rmtree(staging, ignore_errors=True)
                raise RuntimeError("Failure to execute material recording")

            staging.rename(d)

        return Path(d)


@pytest.fixture
def material_recording(material_recording_session: Path, tmp_path: Path):
    target = tmp_path / material_recording_session.name
    shutil.copytree(material_recording_session, target)
    return target
