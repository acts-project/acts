// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <type_traits>

/// Collection of bit masks to enable steering which components of a track state
/// should be initialized, and which should be left invalid.
/// These mask values can be combined using binary operators, so
/// (TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian) will instruct
/// allocating storage for both predicted parameters (including covariance) and
/// a jacobian.
/// The enum is used as a strong type wrapper around the bits to prevent
/// autoconversion from integer
enum struct TrackStatePropMask : uint16_t {
  None = 0,
  BoundPredicted = 1 << 0,
  BoundFiltered = 1 << 1,
  BoundSmoothed = 1 << 2,

  FreePredicted = 1 << 3,
  FreeFiltered = 1 << 4,
  FreeSmoothed = 1 << 5,

  JacobianBoundToBound = 1 << 6,
  JacobianBoundToFree = 1 << 7,
  JacobianFreeToBound = 1 << 8,
  JacobianFreeToFree = 1 << 9,

  Uncalibrated = 1 << 10,
  Calibrated = 1 << 11,

  BoundAll = BoundPredicted + BoundFiltered + BoundSmoothed +
             JacobianBoundToBound + Uncalibrated + Calibrated,
  All = std::numeric_limits<uint16_t>::max(),  // should be all ones
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
