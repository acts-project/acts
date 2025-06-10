// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelToDetectorVolume.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"

#include <numbers>

#include <GeoModelHelpers/GeoShapeUtils.h>
#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoPcon.h>
#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoShapeSubtraction.h>
#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoTubs.h>

namespace Acts::GeoModel {

std::shared_ptr<Volume> convertVolume(const Transform3& trf,
                                      const GeoShape* shape,
                                      VolumeBoundFactory& boundFactory) {
  assert(shape);
  std::shared_ptr<VolumeBounds> bounds{};
  GeoTrf::Transform3D newTrf = trf;
  const ShapeType id = shape->typeID();
  if (id == GeoTube::getClassTypeID()) {
    const auto* tube = dynamic_cast<const GeoTube*>(shape);
    bounds = boundFactory.makeBounds<CylinderVolumeBounds>(
        tube->getRMin(), tube->getRMax(), tube->getZHalfLength());
  } else if (id == GeoTubs::getClassTypeID()) {
    const auto* tubs = dynamic_cast<const GeoTubs*>(shape);
    bounds = boundFactory.makeBounds<CylinderVolumeBounds>(
        tubs->getRMin(), tubs->getRMax(), tubs->getZHalfLength(),
        tubs->getDPhi() / 2);
    newTrf = trf * GeoTrf::RotateZ3D(tubs->getSPhi() + 0.5 * tubs->getDPhi());
  } else if (id == GeoBox::getClassTypeID()) {
    const auto* box = dynamic_cast<const GeoBox*>(shape);
    bounds = boundFactory.makeBounds<CuboidVolumeBounds>(
        box->getXHalfLength(), box->getYHalfLength(), box->getZHalfLength());
  } else if (id == GeoSimplePolygonBrep::getClassTypeID() ||
             /// Union is converted into box. Revise in the future
             id == GeoShapeUnion::getClassTypeID() ||
             /// Will change in future, get bounding box for now
             id == GeoPcon::getClassTypeID()) {
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    shape->extent(xmin, ymin, zmin, xmax, ymax, zmax);
    bounds = boundFactory.makeBounds<CuboidVolumeBounds>(
        (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
  } else if (id == GeoTrd::getClassTypeID()) {
    const auto* trd = dynamic_cast<const GeoTrd*>(shape);
    double x1 = trd->getXHalfLength1();
    double x2 = trd->getXHalfLength2();
    double y1 = trd->getYHalfLength1();
    double y2 = trd->getYHalfLength2();
    double z = trd->getZHalfLength();
    if (y1 == y2) {
      if (x1 <= x2) {
        // y axis in ACTS is z axis in geomodel
        bounds = boundFactory.makeBounds<TrapezoidVolumeBounds>(x1, x2, z, y1);
        constexpr double rotationAngle = std::numbers::pi / 2.;
        newTrf = trf * GeoTrf::RotateX3D(rotationAngle);
      } else {
        bounds = boundFactory.makeBounds<TrapezoidVolumeBounds>(x2, x1, z, y1);
        constexpr double rotationAngle = std::numbers::pi;
        newTrf = trf * GeoTrf::RotateY3D(rotationAngle) *
                 GeoTrf::RotateZ3D(rotationAngle);
      }
    } else if (x1 == x2) {
      if (y1 < y2) {
        bounds = boundFactory.makeBounds<TrapezoidVolumeBounds>(y1, y2, z, x1);
        auto rotationAngle = std::numbers::pi / 2.;
        newTrf = trf * GeoTrf::RotateZ3D(rotationAngle) *
                 GeoTrf::RotateX3D(rotationAngle);
      } else {
        bounds = boundFactory.makeBounds<TrapezoidVolumeBounds>(y2, y1, z, x1);
        auto rotationAngle = std::numbers::pi;
        newTrf = trf * GeoTrf::RotateX3D(rotationAngle) *
                 GeoTrf::RotateZ3D(rotationAngle / 2) *
                 GeoTrf::RotateX3D(rotationAngle / 2);
      }
    } else {
      throw std::runtime_error("convertVolume() - Translating the GeoTrd " +
                               printGeoShape(shape) + " to ACTS failed");
    }

  } else if (id == GeoShapeSubtraction::getClassTypeID()) {
    return convertVolume(newTrf, getOps(shape).first, boundFactory);
  } else if (id == GeoShapeShift::getClassTypeID()) {
    auto compressed = compressShift(shape);
    if (compressed->typeID() != GeoShapeShift::getClassTypeID()) {
      return convertVolume(newTrf, compressed, boundFactory);
    }
    const auto shiftShape =
        dynamic_pointer_cast<const GeoShapeShift>(compressed);
    const GeoShape* shapeOp = shiftShape->getOp();
    return convertVolume(newTrf * shiftShape->getX(), shapeOp, boundFactory);
  } else {
    throw std::runtime_error("Cannot convert " + printGeoShape(shape));
  }
  return std::make_shared<Volume>(newTrf, bounds);
}

std::shared_ptr<Experimental::DetectorVolume> convertDetectorVolume(
    const GeometryContext& context, Volume& vol, const std::string& name,
    const std::vector<std::shared_ptr<Surface>>& sensitives) {
  auto portalGenerator = Experimental::defaultPortalAndSubPortalGenerator();
  return Experimental::DetectorVolumeFactory::construct(
      portalGenerator, context, name, vol.transform(), vol.volumeBoundsPtr(),
      sensitives,
      std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>{},
      Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());
}

}  // namespace Acts::GeoModel
