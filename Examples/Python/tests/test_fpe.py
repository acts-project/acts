import acts
import acts.examples

class FpeMaker(acts.examples.IAlgorithm):

    def __init__(self):
        acts.examples.IAlgorithm.__init__(self, "FpeMaker", acts.logging.INFO)

    def execute(self, context):
        i = context.eventNumber % 3

        if i == 0:
            acts.FpeMonitor._trigger_divbyzero()
        elif i == 1:
            acts.FpeMonitor._trigger_overflow()
        elif i == 2:
            acts.FpeMonitor._trigger_invalid()

        return acts.examples.ProcessCode.SUCCESS

def test_divbyzero():
    s = acts.examples.Sequencer(events=3*10, numThreads=-1)

    s.addAlgorithm(FpeMaker())

    s.run()

    #  with acts.FpeMonitor() as mon:
        #  print(mon)
        #  acts.FpeMonitor._trigger_divbyzero()
        #  print("post 1")
        #  acts.FpeMonitor._trigger_divbyzero()
        #  print("post 2")
        #  acts.FpeMonitor._trigger_overflow()
        #  print("post 3")
        #  acts.FpeMonitor._trigger_invalid()
        #  print("post 4")
        #  pass

        #  print(mon.result)
