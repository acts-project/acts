// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

std::shared_ptr<Experimental::DetectorVolume> convertVolume(
    const GeometryContext& context, const GeoShape& shape,
    const std::string& name, const GeoTrf::Transform3D& transform,
    const std::vector<GeoModelSensitiveSurface>& sensitives) {
  // dummy volume for conversion with surfaces
  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> a;

  // type conversion from GeoModelSensitiveSurface to Surface
  std::vector<std::shared_ptr<Surface>> sensSurfaces(sensitives.size());
  std::transform(sensitives.begin(), sensitives.end(), sensSurfaces.begin(),
                 [](const std::tuple<std::shared_ptr<GeoModelDetectorElement>,
                                     std::shared_ptr<Surface>>& t) {
                   return std::get<1>(t);
                 });
  auto portalGenerator = Experimental::defaultPortalAndSubPortalGenerator();
  if (shape.typeID() == GeoTube::getClassTypeID()) {
    const GeoTube* tube = dynamic_cast<const GeoTube*>(&shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tube->getRMin(), tube->getRMax(),
                                               tube->getZHalfLength());
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces());
  } else if (shape.typeID() == GeoTubs::getClassTypeID()) {
    const GeoTubs* tubs = dynamic_cast<const GeoTubs*>(&shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tubs->getRMin(), tubs->getRMax(),
                                               tubs->getZHalfLength(),
                                               tubs->getDPhi() / 2);
    GeoTrf::Transform3D newTransform =
        transform * GeoTrf::RotateZ3D(tubs->getSPhi() + 0.5 * tubs->getDPhi());
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, newTransform, bounds,
        Experimental::tryAllPortalsAndSurfaces());
  } else if (shape.typeID() == GeoBox::getClassTypeID()) {
    // TODO do the surfaces
    const GeoBox* box = dynamic_cast<const GeoBox*>(&shape);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(box->getXHalfLength(),
                                             box->getYHalfLength(),
                                             box->getZHalfLength());
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds, sensSurfaces, a,
        Experimental::tryNoVolumes(), Experimental::tryAllPortalsAndSurfaces());
  } else if (shape.typeID() == GeoSimplePolygonBrep::getClassTypeID()) {
    const GeoSimplePolygonBrep* brep =
        dynamic_cast<const GeoSimplePolygonBrep*>(&shape);
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    brep->extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces());
  } else if (shape.typeID() == GeoTrd::getClassTypeID()) {
    const GeoTrd* trd = dynamic_cast<const GeoTrd*>(&shape);
    double x1 = trd->getXHalfLength1();
    double x2 = trd->getXHalfLength2();
    double y1 = trd->getYHalfLength1();
    double y2 = trd->getYHalfLength2();
    double z = trd->getZHalfLength();

    if (y1 == y2) {
      if (x1 <= x2) {
        // y axis in ACTS is z axis in geomodel
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x1, x2, z, y1);
        constexpr double rotationAngle = M_PI / 2;
        GeoTrf::Transform3D newTransform =
            transform * GeoTrf::RotateX3D(rotationAngle);
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds, sensSurfaces,
            a, Experimental::tryNoVolumes(),
            Experimental::tryAllPortalsAndSurfaces());
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x2, x1, z, y1);
        constexpr double rotationAngle = M_PI;
        GeoTrf::Transform3D newTransform = transform *
                                           GeoTrf::RotateY3D(rotationAngle) *
                                           GeoTrf::RotateZ3D(rotationAngle);
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds, sensSurfaces,
            a, Experimental::tryNoVolumes(),
            Experimental::tryAllPortalsAndSurfaces());
      }
    } else if (x1 == x2) {
      if (y1 < y2) {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y1, y2, z, x1);
        auto rotationAngle = M_PI / 2;
        GeoTrf::Transform3D newTransform = transform *
                                           GeoTrf::RotateZ3D(rotationAngle) *
                                           GeoTrf::RotateX3D(rotationAngle);
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds, sensSurfaces,
            a, Experimental::tryNoVolumes(),
            Experimental::tryAllPortalsAndSurfaces());
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y2, y1, z, x1);
        auto rotationAngle = M_PI;
        GeoTrf::Transform3D newTransform =
            transform * GeoTrf::RotateX3D(rotationAngle) *
            GeoTrf::RotateZ3D(rotationAngle / 2) *
            GeoTrf::RotateX3D(rotationAngle / 2);
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds, sensSurfaces,
            a, Experimental::tryNoVolumes(),
            Experimental::tryAllPortalsAndSurfaces());
      }
    } else {
      throw std::runtime_error("FATAL: Translating GeoTrd to ACTS failed");
    }
  } else if (shape.typeID() == GeoShapeUnion::getClassTypeID()) {
    const GeoShapeUnion* unionShape =
        dynamic_cast<const GeoShapeUnion*>(&shape);
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    unionShape->extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces());
  } else if (shape.typeID() == GeoShapeSubtraction::getClassTypeID()) {
    // Go down the left side (opA) of the subtraction until we reach a normal
    // shape
    const GeoShapeSubtraction* subtractionShape =
        dynamic_cast<const GeoShapeSubtraction*>(&shape);
    const GeoShape* shapeA = subtractionShape->getOpA();
    return convertVolume(context, *shapeA, name, transform, sensitives);
  } else if (shape.typeID() == GeoPcon::getClassTypeID()) {
    // Will change in future, get bounding box for now
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    const GeoPcon* pcon = dynamic_cast<const GeoPcon*>(&shape);
    pcon->extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces());
  }
  if (shape.typeID() == GeoShapeShift::getClassTypeID()) {
    const GeoShapeShift* shiftShape =
        dynamic_cast<const GeoShapeShift*>(&shape);
    const GeoShape* shapeOp = shiftShape->getOp();
    GeoTrf::Transform3D newTransform = transform * shiftShape->getX();
    return convertVolume(context, *shapeOp, name, newTransform, sensitives);
  }
  throw std::runtime_error("FATAL: Unsupported GeoModel shape");
}

}  // namespace Acts::GeoModel
