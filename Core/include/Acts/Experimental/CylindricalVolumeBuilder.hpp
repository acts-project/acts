// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/CylindricalVolumeHelper.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>
#include <string>
#include <vector>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

template <BinningValue kBINNED = binValues>
struct CylindricalVolumeBuilder {
  /// Building configuration - internal blue prints
  std::vector<InternalBlueprint> internalBlueprints = {};

  /// Building configuration - external volume extent
  GeometricExtent volumeExtent = GeometricExtent{};
  
  /// Naming configuration - the base name of the volume
  std::string name = "Unnamed";
  /// Naming configuration - the suffix for gap volumes
  std::string suffixGapStr = "_r_g";
  /// Naming configuration - the suffix for structure volumes
  std::string suffixVolStr = "_r_l";

  /// Logging configuraiton - the log level for screen output
  Logging::Level logLevel = Logging::INFO;

  /// Call operator that creates a single cylindrical volume
  ///
  /// @param gctx the geometry context for this building call
  ///
  /// @note throws an exception if more than one Internal blueprints
  /// are provided - and if the internal blueprint does not fit into
  /// the external restriction
  ///
  /// @note throws an exception if any other binning structure than
  /// binR, binZ, binPhi or default = binValues is chosen
  ///
  /// @return the new (cylindrical) container volume
  std::vector<std::shared_ptr<DetectorVolume>> operator()(
      const GeometryContext& gctx) noexcept(false) {
    // Single volume case
    if (kBINNED == binValues) {
      // No internal blueprint is given
      if (internalBlueprints.size() == 0) {
        if (not volumeExtent.constrains(binR) and
            not volumeExtent.constrains(binZ)) {
          throw std::invalid_argument(
              "\n *** CylindricalVolumeBuilder: empty volume needs to restrict at "
              "least r "
              "and z.");
        }
        return {DetectorVolume::makeShared(
            CylindricalVolumeHelper::buildTransform(volumeExtent),
            std::move(CylindricalVolumeHelper::buildBounds(volumeExtent)), name)};
      }
      if (internalBlueprints.size() > 1) {
        throw std::invalid_argument(
            "\n *** CylindricalVolumeBuilder: requires at most one InternalBlueprint "
            "object.");
      }
      return {CylindricalVolumeHelper::buildVolume(gctx, internalBlueprints[0],
                                                   volumeExtent, name)};
    }
    // No blue prints but binning requested does not work
    if (internalBlueprints.empty()) {
      throw std::invalid_argument(
          "\n *** CylindricalVolumeBuilder: requires at least one InternalBlueprint "
          "object when binning is requested.");
    }

    std::vector<BinningValue> allowedValues = {binR, binPhi, binZ};
    if (std::find(allowedValues.begin(), allowedValues.end(), kBINNED) ==
        allowedValues.end()) {
      throw std::invalid_argument(
          "\n *** CylindricalVolumeBuilder: volumes must be created in R, z, or phi "
          "when binning is requested.");
    }
    // Remove the adapting value from the allowed ones
    allowedValues.erase(
        std::remove(allowedValues.begin(), allowedValues.end(), kBINNED),
        allowedValues.end());

    // Harmonize the blue prints and extents
    CylindricalVolumeHelper::checkExtents(volumeExtent, internalBlueprints);
    GeometricExtent cExtent = CylindricalVolumeHelper::harmonizeExtents(
        volumeExtent, internalBlueprints, allowedValues);

    // Built the volumes with harmonized z/phi values and adapting in R
    return CylindricalVolumeHelper::buildVolumes(
        gctx, internalBlueprints, cExtent, allowedValues, kBINNED, name,
        suffixVolStr, suffixGapStr, logLevel);
  }
};

}  // namespace Acts