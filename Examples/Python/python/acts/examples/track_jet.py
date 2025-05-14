from typing import Optional, Union, Any, List
from pathlib import Path
from collections import namedtuple
from collections.abc import Iterable

import acts
from acts.examples import TruthJetAlgorithm

TruthJetConfig = namedtuple(
    "TruthJetConfig",
    ["inputTruthParticles", "outputJets"], 
    defaults = [None, None]
        )

def _getTruthJetKWargs(config: TruthJetConfig) -> dict:
    return {
        "inputTruthParticles": config.inputTruthParticles,
        "outputJets": config.outputJets
            }

def addTruthJetAlg(
    s: acts.examples.Sequencer,
    config: TruthJetConfig,
    loglevel: Optional[acts.logging.Level] = None) -> None:

    customLogLevel = acts.examples.defaultLogging(s, loglevel)
    truthJetAlg = acts.examples.TruthJetAlgorithm(
        **acts.examples.defaultKWArgs(**_getTruthJetKWargs(config)),
        level = customLogLevel(),
            )

    s.addAlgorithm(truthJetAlg)
