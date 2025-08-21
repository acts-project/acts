

#pragma once

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include <TFile.h>
#include <TTree.h>

namespace Acts {

void saveScoreMonitor(
    const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>&
        scoreMonitor,
    const std::string& monitorFile,
    const std::vector<std::string>& detectorNames = {} );

}
