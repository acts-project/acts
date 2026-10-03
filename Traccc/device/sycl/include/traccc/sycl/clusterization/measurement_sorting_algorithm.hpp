// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>

namespace traccc::sycl {

/// Algorithm sorting the reconstructed measurements in their container
///
/// The track finding algorithm expects measurements belonging to a single
/// detector module to be consecutive in memory. But
/// @c traccc::sycl::clusterization_algorithm does not (currently) produce the
/// measurements in such an ordered state. This is where this algorithm comes
/// to the rescue.
///
class measurement_sorting_algorithm
    : public algorithm<edm::measurement_collection::buffer(
          const edm::measurement_collection::const_view&)>,
      public messaging {
 public:
  /// Constructor for the algorithm
  ///
  /// @param mr Unused, here for consistency of interface (see CUDA)
  /// @param copy The copy object to use in the algorithm
  /// @param queue Wrapper for the for the SYCL queue for kernel invocation
  ///
  measurement_sorting_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      queue_wrapper& queue,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Callable operator performing the sorting on a container
  ///
  /// @param measurements The measurements to sort
  ///
  [[nodiscard]] output_type operator()(
      const edm::measurement_collection::const_view& measurements)
      const override;

  /// Callable operator performing the sorting, consuming the input
  ///
  /// Unlike the view overload, this one blocks until the device work has
  /// finished, and then releases the input buffer.
  ///
  /// @param measurements The measurements to sort
  ///
  [[nodiscard]] output_type operator()(
      edm::measurement_collection::buffer&& measurements) const;

  /// A const buffer cannot be consumed. Without this overload, moving from
  /// one would silently select the view overload.
  output_type operator()(const edm::measurement_collection::buffer&&) const =
      delete;

 private:
  /// Memory resource(s) to use
  traccc::memory_resource m_mr;
  /// Copy object to use in the algorithm
  std::reference_wrapper<const vecmem::copy> m_copy;
  /// The SYCL queue to use
  std::reference_wrapper<queue_wrapper> m_queue;
};  // class measurement_sorting_algorithm

}  // namespace traccc::sycl
