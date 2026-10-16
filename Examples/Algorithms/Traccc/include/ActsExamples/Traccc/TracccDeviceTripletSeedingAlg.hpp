// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"

#include <memory>
#include <optional>
#include <string>

#include <traccc/seeding/detail/seeding_config.hpp>

namespace ActsExamples {

/// Run the traccc triplet seeding on the device
class TracccDeviceTripletSeedingAlg final : public IAlgorithm {
 public:
  struct Config {
    /// The traccc detector, providing the backend
    std::shared_ptr<const TracccDetector> detector;
    /// Input device traccc spacepoints
    std::string inputSpacepoints = "traccc_device_spacepoints";
    /// Output device traccc seeds
    std::string outputSeeds = "traccc_device_seeds";
    /// Seed finder configuration, including the z component of the magnetic
    /// field (bFieldInZ) used by the seeding
    traccc::seedfinder_config seedFinder;
    /// Seed filter configuration
    traccc::seedfilter_config seedFilter;
  };

  explicit TracccDeviceTripletSeedingAlg(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::optional<
      TracccDeviceAlgorithm<traccc::device::triplet_seeding_algorithm>>
      m_algorithm;

  ReadDataHandle<TracccDeviceSpacepoints> m_inputSpacepoints{
      this, "InputSpacepoints"};
  WriteDataHandle<TracccDeviceSeeds> m_outputSeeds{this, "OutputSeeds"};
};

}  // namespace ActsExamples
