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
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

class GeoFullPhysVol;
class GeoTrd;

namespace Acts {

class Surface;

/// @brief The GeoTrd converter
///
/// This is a dedicated converter for GeoTrd shapes
class GeoTrdConverter : public IGeoShapeConverter {
 public:
  /// @brief Convert a GeoTrd to a detector element and surface
  GeoTrdConverter() = default;
  /// Destructor
  ~GeoTrdConverter() override = default;

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
  /// @brief Convert a GeoTrd to a detector element and surface
  ///
  /// @param geoFPV The full physical volume to convert (contains shape)
  /// @param geoTrd The GeoTrd to convert
  /// @param absTransform from the GeoPhysVol
  /// @param bool sensitive
  ///
  /// @return The detector element and surface
  std::tuple<std::shared_ptr<GeoModelDetectorElement>, std::shared_ptr<Surface>>
  toSurface(const GeoFullPhysVol& geoFPV, const GeoTrd& geoTrd,
            const Transform3& absTransform, bool sensitive) const;
};

}  // namespace Acts
