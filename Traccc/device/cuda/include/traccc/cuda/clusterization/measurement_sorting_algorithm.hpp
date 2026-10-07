// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>

namespace traccc::cuda {

/// Algorithm sorting the reconstructed measurements
///
/// The track finding algorithm expects measurements belonging to a single
/// detector module to be consecutive in memory. But
/// @c traccc::cuda::clusterization_algorithm does not (currently) produce the
/// measurements in such an ordered state. This is where this algorithm comes
/// to the rescue.
///
class measurement_sorting_algorithm
    : public algorithm<edm::measurement_collection::buffer(
          const edm::measurement_collection::const_view&)>,
      public messaging,
      public cuda::algorithm_base {
 public:
  /// Constructor for the algorithm
  ///
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use in the algorithm
  /// @param str The CUDA stream to schedule the measurement sorting in
  /// @param logger The logger to use in the algorithm
  ///
  measurement_sorting_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      const stream_wrapper& str,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Callable operator performing the sorting on a container
  ///
  /// @param measurements The measurements to sort
  ///
  [[nodiscard]] output_type operator()(
      const edm::measurement_collection::const_view& measurements)
      const override;

 private:
  /// The memory resource(s) to use
  traccc::memory_resource m_mr;
  /// Copy object to use in the algorithm
  std::reference_wrapper<const vecmem::copy> m_copy;
};  // class measurement_sorting_algorithm

}  // namespace traccc::cuda
