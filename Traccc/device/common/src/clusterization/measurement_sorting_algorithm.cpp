// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/clusterization/device/measurement_sorting_algorithm.hpp"

// System include(s).
#include <utility>

namespace traccc::device {

measurement_sorting_algorithm::measurement_sorting_algorithm(
    const memory_resource& mr, const vecmem::copy& copy,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base{mr, copy} {}

measurement_sorting_algorithm::output_type
measurement_sorting_algorithm::operator()(
    const edm::measurement_collection::const_view& measurements) const {
  // Exit early if there are no measurements.
  if (measurements.capacity() == 0) {
    return {};
  }

  // Get the number of measurements.
  edm::measurement_collection::const_view::size_type n_measurements = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(measurements, *(mr().host));
    n_measurements = size.get();
  } else {
    n_measurements = copy().get_size(measurements);
  }

  // Create and resize the output buffer.
  output_type result{measurements.capacity(), mr().main,
                     vecmem::data::buffer_type::resizable};
  copy().setup(result)->ignore();
  if (n_measurements == 0) {
    return result;
  }
  copy()(measurements.size(), result.size())->ignore();

  sorting_kernel(
      {static_cast<unsigned int>(n_measurements), measurements, result});

  // Return the sorted buffer.
  return result;
}

}  // namespace traccc::device
