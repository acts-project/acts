// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"

#include "Acts/Plugins/ExaTrkX/ExaTrkXTiming.hpp"

namespace Acts {

void ExaTrkXTrackFindingBase::getTracks(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates,
    LoggerWrapper logger) const {
  auto timeInfo = ExaTrkXTime{};
  getTracks(inputValues, spacepointIDs, trackCandidates, timeInfo, logger);
}

}  // namespace Acts
