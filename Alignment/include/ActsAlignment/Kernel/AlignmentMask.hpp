// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/EnumBitwiseOperators.hpp"

#include <limits>
#include <type_traits>

namespace ActsAlignment {

/// Collection of bit masks to steer which alignment degrees of freedom are
/// fitted and which stay fixed.
/// These mask values can be combined using binary operators, so
/// (AlignmentMask::Center0 | AlignmentMask::Center1) will instruct
/// alignment for translations along the local x and y axes.
/// The enum is used as a strong type wrapper around the bits to prevent
/// autoconversion from integer
enum struct AlignmentMask : std::uint8_t {
  None = 0,
  Center0 = 1 << 0,
  Center1 = 1 << 1,
  Center2 = 1 << 2,
  Rotation0 = 1 << 3,
  Rotation1 = 1 << 4,
  Rotation2 = 1 << 5,

  All = std::numeric_limits<std::uint8_t>::max(),  // should be all ones
};

ACTS_DEFINE_ENUM_BITWISE_OPERATORS(AlignmentMask)

}  // namespace ActsAlignment
