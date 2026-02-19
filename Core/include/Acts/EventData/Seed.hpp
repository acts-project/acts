// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>

namespace Acts {

/// Seed built from N external space points.
template <typename external_space_point_t, std::size_t N = 3ul>
class Seed {
  static_assert(N >= 3ul);

 public:
  /// Type of the external space point
  using value_type = external_space_point_t;
  /// Number of space points in the seed
  static constexpr std::size_t DIM = N;

  /// Constructor from N space points
  /// @param points The space points to build the seed from
  template <typename... args_t>
    requires(sizeof...(args_t) == N) &&
            (std::same_as<external_space_point_t, args_t> && ...)
  explicit Seed(const args_t&... points);

  /// Set the z-coordinate of the vertex
  /// @param vertex The vertex z-coordinate
  void setVertexZ(float vertex);
  /// Set the quality of the seed
  /// @param seedQuality The seed quality value
  void setQuality(float seedQuality);

  /// Get the space points in the seed
  /// @return Array of pointers to the space points
  const std::array<const external_space_point_t*, N>& sp() const;
  /// Get the z-coordinate of the vertex
  /// @return The vertex z-coordinate
  float z() const;
  /// Get the quality of the seed
  /// @return The seed quality value
  float seedQuality() const;

 private:
  std::array<const external_space_point_t*, N> m_spacepoints{};
  float m_vertexZ{0.f};
  float m_seedQuality{-std::numeric_limits<float>::infinity()};
};

}  // namespace Acts

#include "Acts/EventData/Seed.ipp"
