// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BoundFactory.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"

#include <memory>
#include <tuple>

class GeoFullPhysVol;

namespace Acts {
class Surface;
}

namespace ActsPlugins {

/// @addtogroup geomodel_plugin
/// @{

/// @class IGeoShapeConverter
///
/// Interface for the conversion of GeoShapes to Acts surfaces
class IGeoShapeConverter {
 public:
  /// @brief Virtual destructor
  virtual ~IGeoShapeConverter() = default;

  /// @brief Convert a GeoShape into a sensitive surface with associated
  ///        GeoModelDetectorElement
  /// @param geoPV The physical volume to convert
  /// @param transform: Placement of the constructed detector element
  /// @param boundFactory: Reference to the bound factory to share equivalent bounds
  ///                      across multiple surfaces
  /// @return The detector element and surface
  virtual Acts::Result<GeoModelSensitiveSurface> toSensitiveSurface(
      PVConstLink geoPV, const Acts::Transform3& transform,
      Acts::SurfaceBoundFactory& boundFactory) const = 0;

  /// @brief Convert a GeoShape into a passive surface
  /// @param geoPV The physical volume to convert
  /// @param transform: Placement of the constructed detector element
  /// @param boundFactory: Reference to the bound factory to share equivalent bounds
  ///                      across multiple surfaces
  /// @return The detector element and surface
  virtual Acts::Result<std::shared_ptr<Acts::Surface>> toPassiveSurface(
      PVConstLink geoPV, const Acts::Transform3& transform,
      Acts::SurfaceBoundFactory& boundFactory) const = 0;
};

/// @}

}  // namespace ActsPlugins
