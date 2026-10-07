// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/abstract_awaitable.hpp"
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"

// System include(s).
#include <memory>

namespace traccc::device {

/// Base class for algorithms sorting reconstructed measurements.
class measurement_sorting_algorithm
    : public algorithm<edm::measurement_collection::buffer(
          const edm::measurement_collection::const_view&)>,
      public messaging,
      public algorithm_base,
      public virtual abstract_awaitable {
 public:
  /// Data passed to the backend-specific sorting implementation.
  struct measurement_sorting_kernel_payload {
    /// Number of measurements to sort.
    unsigned int n_measurements;
    /// Measurements to sort.
    const edm::measurement_collection::const_view& measurements;
    /// Output buffer receiving the sorted measurements.
    edm::measurement_collection::buffer& output;
  };

  /// Constructor for the measurement-sorting algorithm.
  ///
  /// @param mr The memory resource(s) to use
  /// @param copy The copy object to use
  /// @param logger The logger instance to use
  ///
  measurement_sorting_algorithm(
      const memory_resource& mr, const vecmem::copy& copy,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Sort the measurements by their identifier.
  [[nodiscard]] output_type operator()(
      const edm::measurement_collection::const_view& measurements)
      const override;

 protected:
  /// Launch the backend-specific measurement-sorting implementation.
  virtual void sorting_kernel(
      const measurement_sorting_kernel_payload& payload) const = 0;
};

}  // namespace traccc::device
