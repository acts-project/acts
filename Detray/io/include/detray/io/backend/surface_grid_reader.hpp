// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/grid_builder.hpp"
#include "detray/io/backend/detail/grid_reader.hpp"
#include "detray/io/frontend/payloads.hpp"

// System include(s)
#include <string_view>

namespace detray::io {

/// @brief Surface grid reader backend
///
/// Fills a @c detector_builder from a surface @c detector_grids_payload
template <typename detector_t,
          typename CAP = std::integral_constant<std::size_t, 9>,
          typename DIM = std::integral_constant<std::size_t, 2>>
class surface_grid_reader
    : public detail::grid_reader<detector_t, typename detector_t::surface_type,
                                 grid_builder, CAP, DIM> {
  using grid_reader_t =
      detail::grid_reader<detector_t, typename detector_t::surface_type,
                          grid_builder, CAP, DIM>;
  using base_type = grid_reader_t;

 public:
  /// Tag the reader as "surface_grids"
  static constexpr std::string_view s_tag = "surface_grids";

  /// Payload type that the reader processes
  using payload_type = detector_grids_payload<std::size_t, io::accel_id>;

  /// Same constructors for this class as for base_type
  using base_type::base_type;

  /// @returns the tag of the reader: "homogeneous_material"
  std::string_view tag() const override { return s_tag; }

  /// Convert the detector grids @param grids_data from their IO
  /// payload
  void from_payload(detector_builder<typename detector_t::metadata,
                                     volume_builder>& det_builder,
                    const detector_payload& det_data) const override {
    DETRAY_VERBOSE_HOST("Reading payload object...");

    if (!det_data.surface_grids.has_value()) {
      std::string err_str{"No data in surface grids payload"};
      DETRAY_FATAL_HOST(err_str);
      throw std::invalid_argument(err_str);
    }

    const payload_type& grids_data = *det_data.surface_grids;

    grid_reader_t::from_payload_impl(det_builder, grids_data);
  }
};

}  // namespace detray::io
