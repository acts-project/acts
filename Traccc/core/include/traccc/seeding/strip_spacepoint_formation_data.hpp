// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cstdint>

#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/container.hpp"

namespace traccc {

inline constexpr std::uint64_t invalid_strip_surface_link = UINT64_MAX;

/// Strip-family geometry, independent of detector region.
enum class strip_geometry_model : unsigned int { linear, radial };
enum class strip_local_frame : unsigned int { cartesian, polar };
enum class strip_pair_category : unsigned int { standard, overlap };
enum class strip_pairing_mode : unsigned int { difference, windows };
enum class strip_pairing_coordinate : unsigned int {
  local0,
  local1,
  strip_index
};

/// Map local[component] to scale * floor(local / pitch) + offset.
/// Truncate toward zero to obtain the integer strip index. The adapter selects
/// whether the [0, count) bounds apply before or after truncation.
struct strip_index_mapping {
  unsigned int component{0u};
  unsigned int count{0u};
  scalar pitch{0.f};
  scalar scale{1.f};
  scalar offset{0.f};
  bool check_before_truncation{true};
};

struct linear_strip_geometry {
  unsigned int n_rows{0u};
  unsigned int n_strips{0u};
  scalar pitch{0.f};
  scalar row_length{0.f};
  scalar row_coordinate{0.f};
};

/// One radial strip row, selected by the adapter. Radius bounds are in the
/// surface frame; strip angles and angular pitch are about the strip focus.
struct radial_strip_geometry {
  unsigned int n_strips{0u};
  scalar angular_pitch{0.f};
  scalar min_radius{0.f};
  scalar max_radius{0.f};
  scalar stereo_angle{0.f};
  scalar focal_radius{0.f};
  scalar frame_chord{0.f};
  scalar sin_stereo{0.f};
  scalar cos_stereo{1.f};
  strip_local_frame local_frame{strip_local_frame::cartesian};
};

/// Static geometry records must be strictly sorted by full surface_link.
/// The frame maps (u, v) to origin + u * local_u + v * local_v.
struct strip_measurement_surface_info {
  std::uint64_t surface_link{invalid_strip_surface_link};
  strip_geometry_model geometry_model{strip_geometry_model::linear};
  std::array<scalar, 3> origin{};
  std::array<scalar, 3> local_u{};
  std::array<scalar, 3> local_v{};
  strip_index_mapping index_mapping{};
  linear_strip_geometry linear{};
  radial_strip_geometry radial{};
  scalar spacepoint_variance_r{0.f};
  scalar spacepoint_variance_z{0.f};
};

using strip_measurement_surface_info_collection_types =
    collection_types<strip_measurement_surface_info>;

/// One directed surface pair. Sort by (reference_surface_link,
/// candidate_surface_link), with unique directed pairs and both links present
/// in the geometry collection. A reverse rule creates a separate pair.
/// All windows are inclusive. Experiment-specific cuts are adapter inputs.
struct strip_pairing_rule {
  std::uint64_t reference_surface_link{invalid_strip_surface_link};
  std::uint64_t candidate_surface_link{invalid_strip_surface_link};
  strip_pair_category category{strip_pair_category::standard};
  strip_pairing_mode mode{strip_pairing_mode::difference};
  strip_pairing_coordinate reference_coordinate{
      strip_pairing_coordinate::local0};
  strip_pairing_coordinate candidate_coordinate{
      strip_pairing_coordinate::local0};
  scalar difference_scale{1.f};
  scalar difference_min{0.f};
  scalar difference_max{0.f};
  scalar reference_min{0.f};
  scalar reference_max{0.f};
  scalar candidate_min{0.f};
  scalar candidate_max{0.f};
  /// Physical gap allowance and dimensionless endpoint allowance.
  scalar strip_length_gap_tolerance{0.f};
  scalar strip_length_tolerance{0.f};
};

using strip_pairing_rule_collection_types =
    collection_types<strip_pairing_rule>;

}  // namespace traccc
