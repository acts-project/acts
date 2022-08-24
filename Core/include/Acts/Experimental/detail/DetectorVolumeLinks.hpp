// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {
namespace Experimental {

class DetectorVolume;

namespace detail {

/// @brief a null volume link - explicitely
///
/// @note the method parameters are ignored
///
/// @return the registred volume
inline static const DetectorVolume* nullVolumeLink(
    const GeometryContext& /*ignored*/, const Vector3& /*ignored*/,
    const Vector3& /*ignored*/) {
  return nullptr;
}

/// @brief  The implementation of a single link,
/// The given volume is returned and geometry context,
/// position and direction ignored.
struct SingleLinkImpl {
  /// The single volume to point to
  const DetectorVolume* dVolume = nullptr;

  /// @brief  the function call to be connected with the Delegete
  ///
  /// @note the method parameters are ignored
  ///
  /// @return the registred volume
  const DetectorVolume* targetVolume(const GeometryContext& /*ignored*/,
                                     const Vector3& /*ignored*/,
                                     const Vector3& /*ignored*/) const {
    return dVolume;
  }
};

/// @brief  The implementation of a multi link relationship
/// based on 1D information
///
/// The given volume is returned and geometry context,
/// position and direction ignored.
///
/// @note this does not apply any transform operation
struct MultiLinkCast1DImpl {
  std::vector<const DetectorVolume*> dVolumes = {};
  std::vector<ActsScalar> cBoundaries = {};
  BinningValue bValue;

  /// @brief  the function call to be connected with the Delegete
  ///
  /// @param position is the 3D global position for this cast
  ///
  /// @note the geometry context is ingored
  /// @note the direction is ignored
  ///
  /// @return the registred volume
  const DetectorVolume* targetVolume(const GeometryContext& /*ignored*/,
                                     const Vector3& position,
                                     const Vector3& /*ignored*/) const {
    // cast out the value and, find lower bound and return
    ActsScalar cValue = VectorHelpers::cast(position, bValue);
    auto lb = std::lower_bound(cBoundaries.begin(), cBoundaries.end(), cValue);
    size_t b = static_cast<size_t>(std::distance(cBoundaries.begin(), lb) - 1u);

    return dVolumes[b];
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
