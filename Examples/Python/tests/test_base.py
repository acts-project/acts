import time

import pytest

import acts

import acts.examples


def test_version():
    assert hasattr(acts, "__version__")
    assert hasattr(acts, "version")
    assert hasattr(acts.version, "major")
    assert hasattr(acts.version, "minor")
    assert hasattr(acts.version, "patch")
    assert hasattr(acts.version, "commit_hash")
    assert hasattr(acts.version, "commit_hash_short")


def test_logging():
    for l in ("VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"):
        assert hasattr(acts.logging, l)
        assert hasattr(acts.logging.Level, l)


def test_pgd_particle():
    assert len(acts.PdgParticle.__members__) == 19


def test_algebra():
    v3 = acts.Vector3(1, 2, 3)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2, 3, 4)
    with pytest.raises(TypeError):
        acts.Vector3(1, 2)

    v3 = acts.Vector3([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector3([1, 2])
    with pytest.raises(TypeError):
        acts.Vector3()

    v4 = acts.Vector4(1, 2, 3, 4)
    with pytest.raises(TypeError):
        v4 = acts.Vector4(1, 2, 3)
    v4 = acts.Vector4([1, 2, 3, 4])
    with pytest.raises(TypeError):
        acts.Vector4([1, 2, 3])
    with pytest.raises(TypeError):
        acts.Vector4()


def test_empty_sequencer(conf_const):
    s = acts.examples.Sequencer()
    with pytest.raises(RuntimeError):
        s.run()

    s = conf_const(acts.examples.Sequencer, events=1)
    s.run()


def test_sequencer_single_threaded(ptcl_gun, capfd):
    s = acts.examples.Sequencer(numThreads=1, events=2)
    ptcl_gun(s)
    s.run()
    cap = capfd.readouterr()
    assert cap.err == ""
    assert "Create Sequencer (single-threaded)" in cap.out
    assert "Processed 2 events" in cap.out
    print(cap.out)


def test_sequencer_multi_threaded(ptcl_gun, capfd):
    # This test can use 2 threads (for the 2 events),
    # but could be run single-threaded if threading is not available.
    s = acts.examples.Sequencer(numThreads=-1, events=2)
    ptcl_gun(s)
    s.run()
    cap = capfd.readouterr()
    assert cap.err == ""
    assert "Create Sequencer" in cap.out
    assert "Processed 2 events" in cap.out


class StallAlgorithm(acts.examples.IAlgorithm):
    sleep: float = 1

    def execute(self, ctx):

        if ctx.eventNumber == 50:
            print("BEGIN SLEEP")
            time.sleep(self.sleep)
            print("END OF SLEEP")
        else:
            time.sleep(0.01)

        return acts.examples.ProcessCode.SUCCESS


def test_sequencer_divergence():
    s = acts.examples.Sequencer(
        numThreads=5,
        events=100,
        logLevel=acts.logging.VERBOSE,
        maxInFlightRange=10,
        inFlightSyncTimeoutSeconds=0.5,
    )

    alg = StallAlgorithm(name="stall", level=acts.logging.INFO)
    alg.sleep = 1
    s.addAlgorithm(alg)

    with pytest.raises(RuntimeError) as excinfo:
        with acts.logging.ScopedFailureThreshold(acts.logging.MAX):
            s.run()

    assert "Timeout waiting for in-flight events" in str(excinfo.value)


def test_sequencer_divergence_recover():
    s = acts.examples.Sequencer(
        numThreads=5,
        events=100,
        logLevel=acts.logging.VERBOSE,
        maxInFlightRange=10,
        inFlightSyncTimeoutSeconds=2,
    )

    class StallAlgorithm(acts.examples.IAlgorithm):
        sleep: float = 1

        def execute(self, ctx):

            if ctx.eventNumber == 50:
                print("BEGIN SLEEP")
                time.sleep(self.sleep)
                print("END OF SLEEP")
            else:
                time.sleep(0.01)

            return acts.examples.ProcessCode.SUCCESS

    alg = StallAlgorithm(name="stall", level=acts.logging.INFO)
    alg.sleep = 1
    s.addAlgorithm(alg)

    s.run()


def test_random_number():
    rnd = acts.examples.RandomNumbers(seed=42)


def test_constructors():
    s1 = acts.examples.Sequencer()
    print(s1)
    s2 = acts.examples.Sequencer()
    print(s2)
