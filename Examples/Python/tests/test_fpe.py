import sys

import pytest

import acts
import acts.examples

pytestmark = pytest.mark.skipif(
    sys.platform != "linux", reason="FPE monitoring currently only suported on Linux"
)


class FpeMaker(acts.examples.IAlgorithm):
    def __init__(self, name):
        acts.examples.IAlgorithm.__init__(self, name, acts.logging.INFO)

    def execute(self, context):
        i = context.eventNumber % 4

        if i == 0 or i == 1:
            acts.FpeMonitor._trigger_divbyzero()
        elif i == 2:
            acts.FpeMonitor._trigger_overflow()
        elif i == 3:
            acts.FpeMonitor._trigger_invalid()

        return acts.examples.ProcessCode.SUCCESS


class FuncAlg(acts.examples.IAlgorithm):
    def __init__(self, name, func):
        acts.examples.IAlgorithm.__init__(self, name, acts.logging.INFO)
        self.func = func

    def execute(self, context):
        self.func(context)
        return acts.examples.ProcessCode.SUCCESS


@pytest.fixture(autouse=True)
def disable_log_threshold():
    prev = acts.logging.getFailureThreshold()
    acts.logging.setFailureThreshold(acts.logging.MAX)
    yield
    acts.logging.setFailureThreshold(prev)


def test_notrackfpe():
    s = acts.examples.Sequencer(
        events=3 * 100,
        trackFpes=False,
    )
    s.addAlgorithm(FpeMaker("FpeMaker"))

    s.run()

    res = s.fpeResult

    for x in acts.FpeType.values:
        assert res.count(x) == 0


_names = {
    acts.FpeType.FLTDIV: "DivByZero",
    acts.FpeType.FLTOVF: "Overflow",
    acts.FpeType.FLTINV: "Invalid",
}


_types = [
    pytest.param(acts.FpeType.FLTDIV, id="FLTDIV"),
    pytest.param(acts.FpeType.FLTOVF, id="FLTOVF"),
    pytest.param(acts.FpeType.FLTINV, id="FLTINV"),
]


@pytest.fixture(params=_types)
def fpe_type(request):
    yield request.param


def test_fpe_single_nofail(fpe_type):
    s = acts.examples.Sequencer(
        events=10,
    )

    s.addAlgorithm(
        FuncAlg(
            _names[fpe_type],
            lambda _: getattr(
                acts.FpeMonitor, f"_trigger_{_names[fpe_type].lower()}"
            )(),
        )
    )
    s.run()
    res = s.fpeResult
    for x in acts.FpeType.values:
        assert res.count(x) == (s.config.events if x == fpe_type else 0)


def test_fpe_single_fail(fpe_type):
    s = acts.examples.Sequencer(
        events=10,
        failOnFpe=True,
        numThreads=1,
    )

    s.addAlgorithm(
        FuncAlg(
            _names[fpe_type],
            lambda _: getattr(
                acts.FpeMonitor, f"_trigger_{_names[fpe_type].lower()}"
            )(),
        )
    )

    with pytest.raises(acts.FpeFailure):
        s.run()

    res = s.fpeResult
    for x in acts.FpeType.values:
        assert res.count(x) == (1 if x == fpe_type else 0)


def test_fpe_nocontext():
    class Alg(acts.examples.IAlgorithm):
        def __init__(self):
            acts.examples.IAlgorithm.__init__(self, "Alg", acts.logging.INFO)

        def execute(self, context):
            assert context.fpeMonitor is None
            return acts.examples.ProcessCode.SUCCESS

    s = acts.examples.Sequencer(
        events=10,
        trackFpes=False,
        numThreads=-1,
    )
    s.addAlgorithm(Alg())
    s.run()


def test_fpe_rearm(fpe_type):
    trigger = getattr(acts.FpeMonitor, f"_trigger_{_names[fpe_type].lower()}")

    class Alg(acts.examples.IAlgorithm):
        def __init__(self):
            acts.examples.IAlgorithm.__init__(self, "Alg", acts.logging.INFO)

        def execute(self, context):
            assert context.fpeMonitor is not None
            trigger()
            context.fpeMonitor.rearm()
            trigger()
            return acts.examples.ProcessCode.SUCCESS

    s = acts.examples.Sequencer(
        events=10,
        failOnFpe=False,
        numThreads=-1,
    )
    s.addAlgorithm(Alg())
    s.run()

    res = s.fpeResult
    for x in acts.FpeType.values:
        assert res.count(x) == (s.config.events * 2 if x == fpe_type else 0)


_locs = {
    acts.FpeType.FLTDIV: "acts/Examples/Python/src/ModuleEntry.cpp:85",
    acts.FpeType.FLTOVF: "acts/Examples/Python/src/ModuleEntry.cpp:91",
    acts.FpeType.FLTINV: "acts/Examples/Python/src/ModuleEntry.cpp:97",
}


def test_fpe_masking_single(fpe_type):

    trigger = getattr(acts.FpeMonitor, f"_trigger_{_names[fpe_type].lower()}")

    def func(context):
        trigger()
        assert context.fpeMonitor is not None
        context.fpeMonitor.rearm()
        trigger()

    # Mask, but it's below threshold

    s = acts.examples.Sequencer(
        events=10,
        numThreads=-1,
        failOnFpe=True,
        fpeMasks=[
            acts.examples.Sequencer.FpeMask(_locs[fpe_type], fpe_type, 1),
        ],
    )

    s.addAlgorithm(FuncAlg("Alg", func))

    with pytest.raises(acts.FpeFailure):
        s.run()

    # Mask

    s = acts.examples.Sequencer(
        events=10,
        numThreads=-1,
        failOnFpe=True,
        fpeMasks=[
            acts.examples.Sequencer.FpeMask(_locs[fpe_type], fpe_type, 3),
        ],
    )

    s.addAlgorithm(FuncAlg("Alg", func))

    s.run()

    res = s.fpeResult
    for x in acts.FpeType.values:
        assert res.count(x) == (s.config.events * 2 if x == fpe_type else 0)


def test_buffer_sufficient():
    s = acts.examples.Sequencer(
        events=10000,
    )

    s.addAlgorithm(FuncAlg("Invalid", lambda _: acts.FpeMonitor._trigger_invalid()))
    s.run()

    res = s.fpeResult
    for x in acts.FpeType.values:
        assert res.count(x) == (s.config.events if x == acts.FpeType.FLTINV else 0)
