// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/BoundFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelConversionError.hpp"
#include "ActsPlugins/GeoModel/IGeoShapeConverter.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoTube.h>

namespace ActsPlugins::detail {

template <typename Shape, typename Converter>
struct GenericGeoShapeConverter : public IGeoShapeConverter {
  Acts::Result<GeoModelSensitiveSurface> toSensitiveSurface(
      PVConstLink geoPV, const Acts::Transform3& transform,
      Acts::SurfaceBoundFactory& boundFactory) const override {
    // Retrieve logical volume and absolute transform
    const GeoLogVol* logVol = geoPV->getLogVol();
    if (logVol != nullptr) {
      const GeoShape* geoShape = logVol->getShape();
      auto concreteShape = dynamic_cast<const Shape*>(geoShape);
      if (concreteShape != nullptr) {
        return Converter{}(geoPV, *concreteShape, transform, boundFactory,
                           true);
      }
      return Acts::Result<GeoModelSensitiveSurface>::failure(
          GeoModelConversionError::WrongShapeForConverter);
    }
    return Acts::Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::MissingLogicalVolume);
  }

  Acts::Result<std::shared_ptr<Acts::Surface>> toPassiveSurface(
      PVConstLink geoPV, const Acts::Transform3& transform,
      Acts::SurfaceBoundFactory& boundFactory) const override {
    // Retrieve logical volume and absolute transform
    const GeoLogVol* logVol = geoPV->getLogVol();
    if (logVol != nullptr) {
      const GeoShape* geoShape = logVol->getShape();

      auto concreteShape = dynamic_cast<const Shape*>(geoShape);
      if (concreteShape != nullptr) {
        // Conversion function call with sensitive = false
        auto res =
            Converter{}(geoPV, *concreteShape, transform, boundFactory, false);
        if (!res.ok()) {
          return res.error();
        }

        const auto& [el, surface] = res.value();

        return Acts::Result<std::shared_ptr<Acts::Surface>>::success(surface);
      }
      return Acts::Result<std::shared_ptr<Acts::Surface>>::failure(
          GeoModelConversionError::WrongShapeForConverter);
    }
    return Acts::Result<std::shared_ptr<Acts::Surface>>::failure(
        GeoModelConversionError::MissingLogicalVolume);
  }
};

}  // namespace ActsPlugins::detail
