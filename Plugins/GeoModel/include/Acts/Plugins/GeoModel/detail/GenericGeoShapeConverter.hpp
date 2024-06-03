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
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Plugins/GeoModel/interface/IGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoFullPhysVol.h>

namespace Acts::detail {

template<typename Shape, typename Converter>
struct GenericGeoShapeConverter : public IGeoShapeConverter {

Acts::Result<Acts::GeoModelSensitiveSurface>
toSensitiveSurface(const GeoFullPhysVol& geoFPV) const override {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();
    auto concreteShape = dynamic_cast<const Shape*>(geoShape);
    if (concreteShape != nullptr) {
      return Result<GeoModelSensitiveSurface>::success(
          Converter{}(geoFPV, *concreteShape, transform, true));
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }
  return Result<GeoModelSensitiveSurface>::failure(
      GeoModelConversionError::MissingLogicalVolume);
}

Acts::Result<std::shared_ptr<Acts::Surface>>
toPassiveSurface(const GeoFullPhysVol& geoFPV) const override {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();

    auto concreteShape = dynamic_cast<const Shape*>(geoShape);
    if (concreteShape != nullptr) {
      // Conversion function call with sensitive = false
      auto [element, surface] = Converter{}(geoFPV, *concreteShape, transform, false);
      return Result<std::shared_ptr<Surface>>::success(surface);
    }
    return Result<std::shared_ptr<Surface>>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }
  return Result<std::shared_ptr<Surface>>::failure(
      GeoModelConversionError::MissingLogicalVolume);
}

};

}
