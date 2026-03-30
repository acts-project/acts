import acts
import acts.examples


def test_hello_world_sequence():
    class WhiteBoardKeyInspector(acts.examples.IAlgorithm):
        def __init__(self):
            super().__init__(name="WhiteBoardKeyInspector", level=acts.logging.INFO)
            self.seenEvents = 0

        def execute(self, context):
            keys = context.eventStore.keys
            assert "random_data" in keys
            assert "copied_data" in keys
            self.seenEvents += 1
            return acts.examples.ProcessCode.SUCCESS

    eventsCount = 10
    logLevel = acts.logging.INFO

    rnd = acts.examples.RandomNumbers(seed=42)
    seq = acts.examples.Sequencer(events=eventsCount, numThreads=1, logLevel=logLevel)

    seq.addAlgorithm(acts.examples.HelloLoggerAlgorithm(level=logLevel))

    randomCfg = acts.examples.HelloRandomAlgorithm.Config()
    randomCfg.randomNumbers = rnd
    randomCfg.gaussParameters = (0.0, 2.5)
    randomCfg.uniformParameters = (-1.23, 4.25)
    randomCfg.gammaParameters = (1.0, 1.0)
    randomCfg.drawsPerEvent = 5000
    randomCfg.output = "random_data"
    seq.addAlgorithm(acts.examples.HelloRandomAlgorithm(randomCfg, level=logLevel))

    whiteboardCfg = acts.examples.HelloWhiteBoardAlgorithm.Config()
    whiteboardCfg.input = randomCfg.output
    whiteboardCfg.output = "copied_data"
    seq.addAlgorithm(
        acts.examples.HelloWhiteBoardAlgorithm(whiteboardCfg, level=logLevel)
    )

    inspector = WhiteBoardKeyInspector()
    seq.addAlgorithm(inspector)

    seq.run()

    assert inspector.seenEvents == eventsCount
