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
  CugraphTrackBuilding(std::unique_ptr<const Logger> logger)
      : m_logger(std::move(logger)) {}

  std::vector<std::vector<int>> operator()(std::any nodes, std::any edges,
                                           std::any edge_weights,
                                           std::vector<int> &spacepointIDs,
                                           int deviceHint = -1) override;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }
};

}  // namespace Acts
