// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc::utils {

template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE std::array<float, 2u> to_float_array(
    const detray::dpoint2D<ALGEBRA_TYPE>& point) {
  return {static_cast<float>(getter::element(point, 0u)),
          static_cast<float>(getter::element(point, 1u))};
}

template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE std::array<float, 3u> to_float_array(
    const detray::dpoint3D<ALGEBRA_TYPE>& point) {
  return {static_cast<float>(getter::element(point, 0u)),
          static_cast<float>(getter::element(point, 1u)),
          static_cast<float>(getter::element(point, 2u))};
}

template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE detray::dpoint2D<ALGEBRA_TYPE> to_dpoint2D(
    const std::array<float, 2u>& arr) {
  detray::dpoint2D<ALGEBRA_TYPE> point;
  getter::element(point, 0u) =
      static_cast<typename ALGEBRA_TYPE::value_type>(arr[0u]);
  getter::element(point, 1u) =
      static_cast<typename ALGEBRA_TYPE::value_type>(arr[1u]);
  return point;
}

template <detray::concepts::algebra ALGEBRA_TYPE>
TRACCC_HOST_DEVICE detray::dpoint3D<ALGEBRA_TYPE> to_dpoint3D(
    const std::array<float, 3u>& arr) {
  detray::dpoint3D<ALGEBRA_TYPE> point;
  getter::element(point, 0u) =
      static_cast<typename ALGEBRA_TYPE::value_type>(arr[0u]);
  getter::element(point, 1u) =
      static_cast<typename ALGEBRA_TYPE::value_type>(arr[1u]);
  getter::element(point, 2u) =
      static_cast<typename ALGEBRA_TYPE::value_type>(arr[2u]);
  return point;
}

}  // namespace traccc::utils
