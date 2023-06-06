import acts
import acts.examples

class FpeMaker(acts.examples.IAlgorithm):

    def __init__(self, name):
        acts.examples.IAlgorithm.__init__(self, name, acts.logging.INFO)

    def execute(self, context):
        i = context.eventNumber % 4

        if i == 0 or i ==1:
            acts.FpeMonitor._trigger_divbyzero()
        elif i == 2:
            acts.FpeMonitor._trigger_overflow()
        elif i == 3:
            acts.FpeMonitor._trigger_invalid()

        return acts.examples.ProcessCode.SUCCESS

def test_divbyzero():
    s = acts.examples.Sequencer(events=3*100, numThreads=-1)

    s.addAlgorithm(FpeMaker("FpeMaker"))
    s.addAlgorithm(FpeMaker("FpeMaker2"))

    s.run()

