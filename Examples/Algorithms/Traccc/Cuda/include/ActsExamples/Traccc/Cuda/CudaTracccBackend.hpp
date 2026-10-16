// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Traccc/TracccBackend.hpp"

#include <memory>

namespace ActsExamples {

/// traccc backend running on a CUDA device
///
/// The CUDA specific objects are kept out of this header, so that it can be
/// included without the CUDA toolkit headers.
class CudaTracccBackend final : public TracccBackend {
 public:
  struct Config {
    /// CUDA device to run on
    int device = 0;
  };

  explicit CudaTracccBackend(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);
  ~CudaTracccBackend() override;

  const Config& config() const { return m_cfg; }

  std::string name() const override { return "CUDA"; }

  const traccc::memory_resource& memoryResource() const override;
  vecmem::copy& copy() const override;
  void synchronize() const override;

  traccc::magnetic_field makeMagneticField(
      const traccc::magnetic_field& hostField) const override;

  TracccDeviceAlgorithm<MeasurementSortingAlgorithm> makeMeasurementSorting(
      std::unique_ptr<const Acts::Logger> logger) const override;

  TracccDeviceAlgorithm<
      traccc::device::silicon_pixel_spacepoint_formation_algorithm>
  makeSpacepointFormation(
      std::unique_ptr<const Acts::Logger> logger) const override;

  TracccDeviceAlgorithm<traccc::device::triplet_seeding_algorithm>
  makeTripletSeeding(const traccc::seedfinder_config& finderConfig,
                     const traccc::spacepoint_grid_config& gridConfig,
                     const traccc::seedfilter_config& filterConfig,
                     std::unique_ptr<const Acts::Logger> logger) const override;

  TracccDeviceAlgorithm<traccc::device::seed_parameter_estimation_algorithm>
  makeSeedParameterEstimation(
      const traccc::track_params_estimation_config& config,
      std::unique_ptr<const Acts::Logger> logger) const override;

  TracccDeviceAlgorithm<traccc::device::combinatorial_kalman_filter_algorithm>
  makeCombinatorialKalmanFilter(
      const traccc::finding_config& config,
      std::unique_ptr<const Acts::Logger> logger) const override;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  struct Impl;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
