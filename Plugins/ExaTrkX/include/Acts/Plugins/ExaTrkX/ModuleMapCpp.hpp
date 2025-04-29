// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
#include <string>

namespace Acts {

namespace detail {
class GraphCreatorWrapperBase;
}

class ModuleMapCpp : public GraphConstructionBase {
 public:
  struct Config {
    std::string moduleMapPath;
    float rScale = 1.0;
    float phiScale = 1.0;
    float zScale = 1.0;
    float etaScale = 1.0;
    bool checkModuleConsistencyPerEvent = false;

    bool useGpu = false;
    int gpuDevice = 0;
    int gpuBlocks = 512;
  };

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<detail::GraphCreatorWrapperBase> m_graphCreator;
  std::vector<std::uint64_t> m_uniqueDoupletModuleIds;

  const auto &logger() const { return *m_logger; }

 public:
  ModuleMapCpp(const Config &cfg, std::unique_ptr<const Acts::Logger> logger);
  ~ModuleMapCpp();

  const auto &config() const { return m_cfg; }

  std::tuple<std::any, std::any, std::any> operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<std::uint64_t> &moduleIds,
      const ExecutionContext &execContext = {}) override;
};

}  // namespace Acts
