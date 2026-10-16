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

#include <traccc/bfield/magnetic_field.hpp>
#include <traccc/edm/measurement_collection.hpp>
#include <traccc/finding/device/combinatorial_kalman_filter_algorithm.hpp>
#include <traccc/finding/finding_config.hpp>
#include <traccc/seeding/detail/seeding_config.hpp>
#include <traccc/seeding/detail/track_params_estimation_config.hpp>
#include <traccc/seeding/device/seed_parameter_estimation_algorithm.hpp>
#include <traccc/seeding/device/silicon_pixel_spacepoint_formation_algorithm.hpp>
#include <traccc/seeding/device/triplet_seeding_algorithm.hpp>
#include <traccc/utils/algorithm.hpp>
#include <traccc/utils/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace ActsExamples {

/// A backend specific traccc algorithm, accessed through its backend agnostic
/// base class, together with the copy object it was set up with.
template <typename algorithm_t>
class TracccDeviceAlgorithm {
 public:
  using algorithm_type = algorithm_t;

  TracccDeviceAlgorithm(vecmem::copy& copy,
                        std::shared_ptr<const algorithm_type> algorithm)
      : m_copy(&copy), m_algorithm(std::move(algorithm)) {}

  /// The copy object used by the algorithm
  vecmem::copy& copy() const { return *m_copy; }

  const algorithm_type& operator*() const { return *m_algorithm; }
  const algorithm_type* operator->() const { return m_algorithm.get(); }

 private:
  vecmem::copy* m_copy;
  std::shared_ptr<const algorithm_type> m_algorithm;
};

/// Interface for a traccc device backend (CUDA, HIP, SYCL, Alpaka, ...)
///
/// A backend owns the device memory resources, the copy object and the
/// stream/queue, and constructs the backend specific implementations of the
/// traccc device algorithms. The algorithms using the backend only ever see
/// the backend agnostic traccc base classes, so they do not depend on any
/// backend headers.
///
/// @note All algorithms share a single stream/queue, so the sequencer has to
///       run single threaded.
class TracccBackend {
 public:
  using MeasurementSortingAlgorithm =
      traccc::algorithm<traccc::edm::measurement_collection::buffer(
          const traccc::edm::measurement_collection::const_view&)>;

  virtual ~TracccBackend() = default;

  /// Human readable name of the backend
  virtual std::string name() const = 0;

  /// Device (main) and host-accessible memory resources
  virtual const traccc::memory_resource& memoryResource() const = 0;

  /// Copy object for host <-> device transfers on the backend stream/queue
  virtual vecmem::copy& copy() const = 0;

  /// Block until all work submitted to the backend stream/queue is done
  virtual void synchronize() const = 0;

  /// Make a copy of a host magnetic field usable on the device
  virtual traccc::magnetic_field makeMagneticField(
      const traccc::magnetic_field& hostField) const = 0;

  /// @name Algorithm factories
  /// @{
  virtual TracccDeviceAlgorithm<MeasurementSortingAlgorithm>
  makeMeasurementSorting(std::unique_ptr<const Acts::Logger> logger) const = 0;

  virtual TracccDeviceAlgorithm<
      traccc::device::silicon_pixel_spacepoint_formation_algorithm>
  makeSpacepointFormation(std::unique_ptr<const Acts::Logger> logger) const = 0;

  virtual TracccDeviceAlgorithm<traccc::device::triplet_seeding_algorithm>
  makeTripletSeeding(const traccc::seedfinder_config& finderConfig,
                     const traccc::spacepoint_grid_config& gridConfig,
                     const traccc::seedfilter_config& filterConfig,
                     std::unique_ptr<const Acts::Logger> logger) const = 0;

  virtual TracccDeviceAlgorithm<
      traccc::device::seed_parameter_estimation_algorithm>
  makeSeedParameterEstimation(
      const traccc::track_params_estimation_config& config,
      std::unique_ptr<const Acts::Logger> logger) const = 0;

  virtual TracccDeviceAlgorithm<
      traccc::device::combinatorial_kalman_filter_algorithm>
  makeCombinatorialKalmanFilter(
      const traccc::finding_config& config,
      std::unique_ptr<const Acts::Logger> logger) const = 0;
  /// @}
};

}  // namespace ActsExamples
