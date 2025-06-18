// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace Acts {

class CudaTrackBuilding final : public Acts::TrackBuildingBase {
 public:
  struct Config {
    bool useOneBlockImplementation = true;
    bool doJunctionRemoval = false;
  };

  CudaTrackBuilding(const Config &cfg, std::unique_ptr<const Logger> logger)
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  std::vector<std::vector<int>> operator()(
      PipelineTensors tensors, std::vector<int> &spacepointIDs,
      const ExecutionContext &execContext = {}) override;
  const Config &config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  const auto &logger() const { return *m_logger; }
};

}  // namespace Acts
