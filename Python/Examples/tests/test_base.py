import pytest

import acts

import acts.examples


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


def test_random_number():
    rnd = acts.examples.RandomNumbers(seed=42)


def test_constructors():
    s1 = acts.examples.Sequencer()
    print(s1)
    s2 = acts.examples.Sequencer()
    print(s2)
