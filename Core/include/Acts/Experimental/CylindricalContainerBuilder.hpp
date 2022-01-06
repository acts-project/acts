// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class DetectorVolume;

using VolumeBuilder =
    std::function<std::vector<std::shared_ptr<DetectorVolume>>(
        const GeometryContext&)>;

// This is a mockup builder that simple provided pre-built volumes
struct CylindricalVolumeProvider {
  std::vector<std::shared_ptr<DetectorVolume>> volumes;

  // Return the pre-built vector and ignore the gometry context
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& /*ignored*/) {
    return volumes;
  }
};

struct CylindricalContainerBuilderR {
  /// Call operator to create a cylindrical container ordered in R
  ///
  /// @param gctx the current geometry context
  /// @param volumeBuilder the configured volume builder
  /// @param name the name of the new volume
  ///
  /// @return the new (cylindrical) container volume
  static std::shared_ptr<DetectorVolume> build(
      const GeometryContext& gctx, const VolumeBuilder& volumeBuilder,
      const std::string& name);
};

struct CylindricalContainerBuilderZ {
  /// Volume builder for the creation
  VolumeBuilder volumeBuilder;

  /// Name of the builder
  std::string name = "Unnamed";

  /// Call operator to create a cylindrical container ordered in Z
  ///
  /// @param gctx the current geometry context
  /// @param volumeBuilder the configured volume builder
  /// @param name the name of the new volume
  ///
  /// @return the new (cylindrical) container volume
  static std::shared_ptr<DetectorVolume> build(
      const GeometryContext& gctx, const VolumeBuilder& volumeBuilder,
      const std::string& name);
};

struct CylindricalContainerBuilderPhi {
  /// Call operator to create a cylindrical container ordered in Z
  ///
  /// @param gctx the current geometry context
  /// @param volumeBuilder the configured volume builder
  /// @param name the name of the new volume
  ///
  /// @return the new (cylindrical) container volume
  static std::shared_ptr<DetectorVolume> build(
      const GeometryContext& gctx, const VolumeBuilder& volumeBuilder,
      const std::string& name);
};

template <BinningValue kBINNED>
struct CylindricalContainerBuilder {
  const std::vector<BinningValue> allowedValues = {binR, binZ, binPhi};

  /// Volume builder for the creation
  VolumeBuilder volumeBuilder;

  /// Name of the builder
  std::string name = "Unnamed";

  /// Call operator to create a cylindrical container ordered
  /// eith  R, Z or phi
  ///
  /// @param gctx the current geometry context
  ///
  /// @note throws exception if a non-allowed value is used
  ///
  /// @return the new (cylindrical) container volume
  std::shared_ptr<DetectorVolume> operator()(
      const GeometryContext& gctx) noexcept(false) {
    if (std::find(allowedValues.begin(), allowedValues.end(), kBINNED) ==
        allowedValues.end()) {
      throw std::invalid_argument(
          "CylindricalVolumeBuilder can only be constructed with binning in r, "
          "z or phi.");
    }

    if (kBINNED == binR) {
      return CylindricalContainerBuilderR::build(gctx, volumeBuilder, name);
    } else if (kBINNED == binZ) {
      return CylindricalContainerBuilderZ::build(gctx, volumeBuilder, name);
    }
    return CylindricalContainerBuilderPhi::build(gctx, volumeBuilder, name);
  }
};

}  // namespace Acts
