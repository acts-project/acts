from typing import Optional, Dict

import acts
from acts.examples import AliasAlgorithm


def addAlias(
    s: acts.examples.Sequencer,
    mapping: Dict[str, str],
    logLevel: Optional[acts.logging.Level] = None,
):
    """ """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    s.addAlgorithm(
        AliasAlgorithm(
            mapping=mapping,
            level=customLogLevel(),
        )
    )
