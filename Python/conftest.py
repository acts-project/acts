import acts
import acts.examples

# The pytest suite drives whole Sequencer jobs whose loggers it cannot reach
# individually, so it arms the process-wide default instead.
acts.examples.setLogFailureThreshold(acts.logging.WARNING)
