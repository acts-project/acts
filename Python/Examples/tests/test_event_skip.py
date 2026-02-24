from pathlib import Path

import acts
import acts.examples


class PyAlgorithm(acts.examples.IAlgorithm):
    def __init__(self, name: str, output_file: Path, skip_even: bool):
        acts.examples.IAlgorithm.__init__(self, name, acts.logging.INFO)
        self.name = name
        self.output_file = output_file 
        self.skip_even = skip_even        

    def execute(self, context):

        with self.output_file.open("a", encoding="utf-8") as out:
            out.write(f"{context.eventNumber}\n")
            print(f"{self.name} - writing  event {context.eventNumber}")

        if context.eventNumber % 2 == 0:
            print(f"{self.name} - skipping event {context.eventNumber}")
            return acts.examples.ProcessCode.SKIP

        return acts.examples.ProcessCode.SUCCESS

def _read_sequence(path: Path) -> list[int]:
    with path.open(encoding="utf-8") as f:
        return [int(line.strip()) for line in f if line.strip()]


def test_event_skip(tmp_path):
    output_a = tmp_path / "algorithm_a.txt"
    output_b = tmp_path / "algorithm_b.txt"

    seq = acts.examples.Sequencer(events=20, numThreads=1)
    seq.addAlgorithm(PyAlgorithm("AlgorithmA", output_a, skip_even=True))
    seq.addAlgorithm(PyAlgorithm("AlgorithmB", output_b, skip_even=False))
    seq.run()

    sequence_a = _read_sequence(output_a)
    sequence_b = _read_sequence(output_b)

    print(sequence_a)
    print(sequence_b)

    assert sequence_a == list(range(20))
    assert sequence_b == list(range(1, 20, 2))
