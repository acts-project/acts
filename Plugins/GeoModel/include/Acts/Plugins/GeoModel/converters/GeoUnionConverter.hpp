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
#include "Acts/Plugins/GeoModel/converters/GeoShiftConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/LineBounds.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeUnion.h>

#include "GeoTrdConverter.hpp"
#include "GeoShiftConverter.hpp"

namespace Acts {

/// This is super hacky and ugly
/// TODO Please to not merge!!!
struct GeoUnionConverter : public IGeoShapeConverter {
  bool useA = true;

  auto convert(const GeoFullPhysVol& geoFPV, const GeoShape& shape, const Transform3& transform,
               bool sensitive) const -> Result<GeoModelSensitiveSurface> {
    if (auto trd = dynamic_cast<const GeoTrd*>(&shape); trd != nullptr) {
      return detail::GeoTrdConverter{}(geoFPV, *trd, transform, sensitive);
    }
    if (auto shift = dynamic_cast<const GeoShapeShift*>(&shape); shift != nullptr) {
      return detail::GeoShiftConverter{}(geoFPV, *shift, transform, sensitive);
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }

  Acts::Result<Acts::GeoModelSensitiveSurface> toSensitiveSurface(
      const GeoFullPhysVol& geoFPV) const override {
    // Retrieve logcal volume and absolute transform
    const GeoLogVol* logVol = geoFPV.getLogVol();
    const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
    if (logVol != nullptr) {
      const GeoShape* geoShape = logVol->getShape();

      // This should be only called when this is clear
      // As I said, super hacky
      auto concreteShape = static_cast<const GeoShapeUnion*>(geoShape);

      auto pick = useA ? concreteShape->getOpA() : concreteShape->getOpB();

      if (concreteShape != nullptr) {
        return convert(geoFPV, *pick, transform, true);
      }
      return Result<GeoModelSensitiveSurface>::failure(
          GeoModelConversionError::WrongShapeForConverter);
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::MissingLogicalVolume);
  }

  Acts::Result<std::shared_ptr<Acts::Surface>> toPassiveSurface(
      const GeoFullPhysVol&) const override {
    throw std::runtime_error("not implemented");
  }
};

}  // namespace Acts
