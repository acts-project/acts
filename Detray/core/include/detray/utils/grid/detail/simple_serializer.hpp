// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"

namespace detray {

/// @brief Transforms a local axis bin index to a global grid bin index and vice
/// versa
///
/// Serializers allow to create a memory local data layout if advantageous.
///
/// @note the bin indices are expected to start at zero.
template <std::size_t kDIM>
struct simple_serializer {};

/// @brief Simple serializer specialization for a single axis
template <>
struct simple_serializer<1> {
  /// @returns the axis local bin, which is also the global bin
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(
      multi_axis_t & /*axes*/, typename multi_axis_t::loc_bin_index mbin) const
      -> dindex {
    return mbin[0];
  }

  /// @returns the global bin, which is also the axis local bin
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(multi_axis_t & /*axes*/, dindex gbin) const
      -> typename multi_axis_t::loc_bin_index {
    return {gbin};
  }
};

/// @brief Simple serializer specialization for a 3D multi-axis
template <>
struct simple_serializer<2> {
  /// @brief Create a serial bin from a multi-bin - 3D
  ///
  /// @tparam multi_axis_t is the type of multi-dimensional axis
  ///
  /// @param axes contains all axes (multi-axis)
  /// @param mbin contains a bin index for every axis in the multi-axis.
  ///
  /// @returns a dindex for the bin data storage
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(
      multi_axis_t &axes, typename multi_axis_t::loc_bin_index mbin) const
      -> dindex {
    dindex offset{mbin[1] * axes.template get_axis<0>().nbins()};
    return offset + mbin[0];
  }

  /// @brief Create a bin tuple from a serialized bin - 3D
  ///
  /// @tparam multi_axis_t is the type of multi-dimensional axis
  ///
  /// @param axes contains all axes (multi-axis)
  /// @param gbin the global (serial) bin
  ///
  /// @return a 2-dimensional multi-bin
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes, dindex gbin) const ->
      typename multi_axis_t::loc_bin_index {
    dindex nbins_axis0 = axes.template get_axis<0>().nbins();

    dindex bin0{gbin % nbins_axis0};
    dindex bin1{gbin / nbins_axis0};

    return {bin0, bin1};
  }
};

/// @brief Simple serializer specialization for a 3D multi-axis
template <>
struct simple_serializer<3> {
  /// @brief Create a serial bin from a multi-bin - 3D
  ///
  /// @tparam multi_axis_t is the type of multi-dimensional axis
  ///
  /// @param axes contains all axes (multi-axis)
  /// @param mbin contains a bin index for every axis in the multi-axis.
  ///
  /// @returns a dindex for the bin data storage
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(
      multi_axis_t &axes, typename multi_axis_t::loc_bin_index mbin) const
      -> dindex {
    dindex nbins_axis0 = axes.template get_axis<0>().nbins();
    dindex nbins_axis1 = axes.template get_axis<1>().nbins();

    dindex offset1{mbin[1] * nbins_axis0};
    dindex offset2{mbin[2] * (nbins_axis0 * nbins_axis1)};

    return offset1 + offset2 + mbin[0];
  }

  /// @brief Create a bin tuple from a serialized bin - 3D
  ///
  /// @tparam multi_axis_t is the type of multi-dimensional axis
  ///
  /// @param axes contains all axes (multi-axis)
  /// @param gbin the global (serial) bin
  ///
  /// @return a 3-dimensional multi-bin
  template <typename multi_axis_t>
  DETRAY_HOST_DEVICE auto operator()(multi_axis_t &axes, dindex gbin) const ->
      typename multi_axis_t::loc_bin_index {
    dindex nbins_axis0 = axes.template get_axis<0>().nbins();
    dindex nbins_axis1 = axes.template get_axis<1>().nbins();

    dindex bin0{gbin % nbins_axis0};
    dindex bin1{(gbin / nbins_axis0) % nbins_axis1};
    dindex bin2{gbin / (nbins_axis0 * nbins_axis1)};

    return {bin0, bin1, bin2};
  }
};

}  // namespace detray
