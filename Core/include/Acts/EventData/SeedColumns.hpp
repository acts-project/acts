// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/EnumBitwiseOperators.hpp"

#include <cstdint>

namespace Acts {

/// Enumeration of available columns for space point data storage
enum class SeedColumns : std::uint32_t {
  None = 0,  ///< No columns

  SpacePointIndices =
      1 << 0,        ///< Indices of space points associated with the seed
  Quality = 1 << 1,  ///< Quality of the seed
  VertexZ = 1 << 2,  ///< Z coordinate of the vertex associated with

  /// All columns
  All = SpacePointIndices | Quality | VertexZ,
};

/// Enable bitwise operators for SeedColumns enum
ACTS_DEFINE_ENUM_BITWISE_OPERATORS(SeedColumns);

}  // namespace Acts
