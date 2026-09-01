import acts

# The pytest suite drives whole Sequencer jobs whose loggers it cannot reach
# individually, so it arms the process-wide default instead.
acts.logging.setFailureThreshold(acts.logging.WARNING)
