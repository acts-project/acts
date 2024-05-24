// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/interface/IGeoShapeConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

class GeoFullPhysVol;
class GeoTube;

namespace Acts {

/// @brief The GeoTube converter
///
/// This is a dedicated converter for GeoTube shapes
class GeoTubeConverter : public IGeoShapeConverter {
 public:
  Surface::SurfaceType targetShape = Surface::SurfaceType::Straw;

  /// @brief Convert a GeoTube to a detector element and surface
  GeoTubeConverter() = default;
  /// Destructor
  ~GeoTubeConverter() override = default;

  /// @brief Convert a GeoFullPhysVol to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  ///
  /// @note the Result will be not ok if the provided GeoShape is not
  /// of the extpected type by the converter
  ///
  /// @return The detector element and surface
  Result<GeoModelSensitiveSurface> toSensitiveSurface(
      const GeoFullPhysVol& geoFPV) const override;

  /// @brief Convert a GeoFullPhysVol to a detector element and passive surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  ///
  /// @note the Result will be not ok if the provided GeoShape is not
  /// of the extpected type by the converter
  ///
  /// @return The representing surface
  Result<std::shared_ptr<Surface>> toPassiveSurface(
      const GeoFullPhysVol& geoFPV) const override;

 private:
  /// @brief Convert a GeoTube to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoTube The GeoTube shape to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  std::tuple<std::shared_ptr<GeoModelDetectorElement>, std::shared_ptr<Surface>>
  toSurface(const GeoFullPhysVol& geoFPV, const GeoTube& geoTube,
            const Transform3& absTransform, bool sensitive) const;
};

}  // namespace Acts
