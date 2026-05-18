// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>
#include <unordered_map>

#include <traccc/edm/spacepoint_collection.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

namespace ActsExamples {

class ActsSpToTracccAlg final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSpacePoints = "spacepoints";
    std::string outputTracccSpacepoints = "acts-traccc-spacepoints";
  };

  explicit ActsSpToTracccAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

  mutable vecmem::host_memory_resource m_mr;

 private:
  Config m_cfg;

  ReadDataHandle<SpacePointContainer> m_inputSpacePoints{this,
                                                         "inputSpacePoints"};
  WriteDataHandle<traccc::edm::spacepoint_collection::host>
      m_outputTracccSpacepoints{this, "outputTracccSpacepoints"};
};

}  // namespace ActsExamples
