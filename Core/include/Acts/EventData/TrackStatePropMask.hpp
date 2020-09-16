// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <type_traits>

namespace Acts {

/// Collection of bit masks to enable steering which components of a track state
/// should be initialized, and which should be left invalid.
/// These mask values can be combined using binary operators, so
/// (TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian) will instruct
/// allocating storage for both predicted parameters (including covariance) and
/// a jacobian.
/// The enum is used as a strong type wrapper around the bits to prevent
/// autoconversion from integer
enum struct TrackStatePropMask : uint8_t {
  None = 0,
  Predicted = 1 << 0,
  Filtered = 1 << 1,
  Smoothed = 1 << 2,
  Jacobian = 1 << 3,

  Uncalibrated = 1 << 4,
  Calibrated = 1 << 5,

  All = std::numeric_limits<uint8_t>::max(),  // should be all ones
};

constexpr TrackStatePropMask operator|(TrackStatePropMask lhs,
                                       TrackStatePropMask rhs) {
  return static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) |
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
}

constexpr TrackStatePropMask operator&(TrackStatePropMask lhs,
                                       TrackStatePropMask rhs) {
  return static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) &
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
}

constexpr TrackStatePropMask operator^(TrackStatePropMask lhs,
                                       TrackStatePropMask rhs) {
  return static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) ^
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
}

constexpr TrackStatePropMask operator~(TrackStatePropMask op) {
  return static_cast<TrackStatePropMask>(
      ~static_cast<std::underlying_type<TrackStatePropMask>::type>(op));
}

constexpr TrackStatePropMask& operator|=(TrackStatePropMask& lhs,
                                         TrackStatePropMask rhs) {
  lhs = static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) |
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
  return lhs;
}

constexpr TrackStatePropMask& operator&=(TrackStatePropMask& lhs,
                                         TrackStatePropMask rhs) {
  lhs = static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) &
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
  return lhs;
}

constexpr TrackStatePropMask& operator^=(TrackStatePropMask& lhs,
                                         TrackStatePropMask rhs) {
  lhs = static_cast<TrackStatePropMask>(
      static_cast<std::underlying_type<TrackStatePropMask>::type>(lhs) ^
      static_cast<std::underlying_type<TrackStatePropMask>::type>(rhs));
  return lhs;
}

}  // namespace Acts