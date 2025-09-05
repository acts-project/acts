// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/ScoreBasedAmbiguityResolution.hpp"

#include <TFile.h>
#include <TTree.h>

namespace Acts {

void saveScoreMonitor(
    const std::vector<Acts::ScoreBasedAmbiguityResolution::ScoreMonitor>&
        scoreMonitor,
    const std::string& monitorFile,
    const std::vector<std::string>& detectorNames = {},
    const std::vector<std::string>& optionalFunctions = {});

}
