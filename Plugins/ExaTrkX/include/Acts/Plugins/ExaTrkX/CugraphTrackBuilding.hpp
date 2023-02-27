// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace Acts {

class CugraphTrackBuilding final : public Acts::TrackBuildingBase {
 public:
  std::vector<std::vector<int>> operator()(std::any nodes, std::any edges,
                                           std::any edge_weights,
                                           std::vector<int> &spacepointIDs,
                                           const Logger &logger) override;
};

}  // namespace Acts
