// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

namespace ActsExamples::Traccc {

struct TracccChainConfig {
  enum class Backend { CPU, CUDA };
  Backend backend = Backend::CPU;

  // Detector files
  std::string detectorFile;
  std::string digitizationFile;
  std::string conditionsFile;
  std::string gridFile;
  std::string materialFile;

  // Magnetic field file
  std::string magneticFieldFile;
  // z-component used by CPU track_params_estimation
  float bFieldInZ = 2.0f;

  // Algorithm configs
  traccc::clustering_config clusteringConfig;
  traccc::seedfinder_config seedfinderConfig;
  traccc::seedfilter_config seedfilterConfig;
  traccc::track_params_estimation_config trackParamsEstConfig;
  traccc::finding_config findingConfig;
  traccc::fitting_config fittingConfig;
  traccc::host::greedy_ambiguity_resolution_algorithm::config_type
      resolutionConfig;
};

struct EventResult {
  std::size_t n_measurements = 0;
  std::size_t n_spacepoints = 0;
  std::size_t n_seeds = 0;
  std::size_t n_track_candidates = 0;
  std::size_t n_fitted_tracks = 0;
};

class TracccChain {
 public:
  struct Impl;

  explicit TracccChain(const TracccChainConfig& cfg,
                       Acts::Logging::Level logLevel = Acts::Logging::INFO);
  ~TracccChain();

  TracccChain(TracccChain&&) noexcept;
  TracccChain& operator=(TracccChain&&) noexcept;

  EventResult operator()(
      traccc::edm::measurement_collection::host& measurements,
      traccc::edm::spacepoint_collection::host& spacepoints) const;

 private:
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples::Traccc
