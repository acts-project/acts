// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <string_view>

namespace Acts {

namespace detail_tsp {
inline constexpr HashedString kPreviousKey = hashString("previous");
inline constexpr HashedString kChi2Key = hashString("chi2");
inline constexpr HashedString kPathLengthKey = hashString("pathLength");
inline constexpr HashedString kTypeFlagsKey = hashString("typeFlags");
inline constexpr HashedString kPredictedKey = hashString("predicted");
inline constexpr HashedString kFilteredKey = hashString("filtered");
inline constexpr HashedString kSmoothedKey = hashString("smoothed");
inline constexpr HashedString kJacobianKey = hashString("jacobian");
inline constexpr HashedString kProjectorKey = hashString("projector");
inline constexpr HashedString kUncalibratedKey =
    hashString("uncalibratedSourceLink");
inline constexpr HashedString kCalibratedKey = hashString("calibrated");
inline constexpr HashedString kCalibratedCovKey = hashString("calibratedCov");
inline constexpr HashedString kMeasDimKey = hashString("measdim");
inline constexpr HashedString kNextKey = hashString("next");

}  // namespace detail_tsp

namespace detail_tsp {

template <typename Derived, bool read_only>
class TrackStateProxyCommon {
 protected:
  using IndexType = Acts::TrackIndexType;
  // using ConstParameters = const_parameters_t;
  // using Parameters = parameters_t;
  // using ConstCovariance = const_covariance_t;
  // using Covariance = covariance_t;
  //
  // constexpr Derived& derived() { return static_cast<Derived&>(*this); }
  // constexpr const Derived& derived() const {
  //   return static_cast<const Derived&>(*this);
  // }

 public:
};

}  // namespace detail_tsp
}  // namespace Acts
