// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelToDetVol.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
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

namespace Acts {
namespace GeoModel {
void convertVolume(
    const GeometryContext& context, const GeoShape& shape,
    const std::string& name, const Transform3& transform,
    std::vector<std::shared_ptr<Experimental::DetectorVolume>>& volumes) {
  auto portalGenerator = Experimental::defaultPortalAndSubPortalGenerator();
  if (shape.typeID() == GeoTube::getClassTypeID()) {
    const GeoTube& tube = static_cast<const GeoTube&>(shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tube.getRMin(), tube.getRMax(),
                                               tube.getZHalfLength());
    volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces()));
  } else if (shape.typeID() == GeoTubs::getClassTypeID()) {
    const GeoTubs& tubs = static_cast<const GeoTubs&>(shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tubs.getRMin(), tubs.getRMax(),
                                               tubs.getZHalfLength(),
                                               tubs.getDPhi() / 2);
    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    newTransform.translate(transform.translation());
    newTransform.rotate(transform.rotation() *
                        Eigen::AngleAxisd(tubs.getSPhi() + 0.5 * tubs.getDPhi(),
                                          Acts::Vector3::UnitZ()));
    volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, newTransform, bounds,
        Experimental::tryAllPortalsAndSurfaces()));
  } else if (shape.typeID() == GeoBox::getClassTypeID()) {
    const GeoBox& box = static_cast<const GeoBox&>(shape);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            box.getXHalfLength(), box.getYHalfLength(), box.getZHalfLength());
    volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces()));
  } else if (shape.typeID() == GeoSimplePolygonBrep::getClassTypeID()) {
    const GeoSimplePolygonBrep& brep =
        static_cast<const GeoSimplePolygonBrep&>(shape);
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    brep.extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
    volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces()));
  } else if (shape.typeID() == GeoTrd::getClassTypeID()) {
    const GeoTrd& trd = static_cast<const GeoTrd&>(shape);
    double x1 = trd.getXHalfLength1();
    double x2 = trd.getXHalfLength2();
    double y1 = trd.getYHalfLength1();
    double y2 = trd.getYHalfLength2();
    double z = trd.getZHalfLength();
    if (y1 == y2) {
      if (x1 <= x2) {
        // y axis in ACTS is z axis in geomodel
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x1, x2, z, y1);
        auto rotationAngle = M_PI / 2;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            transform.rotation() *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()));
        volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortalsAndSurfaces()));
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x2, x1, z, y1);
        auto rotationAngle = M_PI;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            transform.rotation() *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitY()) *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitZ()));
        volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortalsAndSurfaces()));
      }
    } else if (x1 == x2) {
      if (y1 < y2) {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y1, y2, z, x1);
        auto rotationAngle = M_PI / 2;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            transform.rotation() *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitZ()) *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()));
        volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortalsAndSurfaces()));
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y2, y1, z, x1);
        auto rotationAngle = M_PI;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            transform.rotation() *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()) *
            Eigen::AngleAxisd(rotationAngle / 2, Acts::Vector3::UnitZ()) *
            Eigen::AngleAxisd(rotationAngle / 2, Acts::Vector3::UnitX()));
        volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortalsAndSurfaces()));
      }
    } else {
      throw std::runtime_error("FATAL: Translating GeoTrd to ACTS failed");
    }
  } else if (shape.typeID() == GeoShapeUnion::getClassTypeID()) {
    const GeoShapeUnion& unionShape = static_cast<const GeoShapeUnion&>(shape);
    const GeoShape* shapeA = unionShape.getOpA();
    const GeoShape* shapeB = unionShape.getOpB();
    convertVolume(context, *shapeA, name + "A", transform, volumes);
    convertVolume(context, *shapeB, name + "B", transform, volumes);
  } else if (shape.typeID() == GeoShapeSubtraction::getClassTypeID()) {
    // Go down the left side (opA) of the subtraction until we reach a normal
    // shape
    const GeoShapeSubtraction& subtractionShape =
        static_cast<const GeoShapeSubtraction&>(shape);
    const GeoShape* shapeA = subtractionShape.getOpA();
    convertVolume(context, *shapeA, name, transform, volumes);
  } else if (shape.typeID() == GeoPcon::getClassTypeID()) {
    // Will change in future, get bounding box for now
    double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
    const GeoPcon& pcon = static_cast<const GeoPcon&>(shape);
    pcon.extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            (xmax - xmin) / 2, (ymax - ymin) / 2, (zmax - zmin) / 2);
    volumes.emplace_back(Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortalsAndSurfaces()));
  }
  if (shape.typeID() == GeoShapeShift::getClassTypeID()) {
    const GeoShapeShift& shiftShape = static_cast<const GeoShapeShift&>(shape);
    const GeoShape* shapeOp = shiftShape.getOp();
    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    newTransform.translate(transform.translation());
    newTransform.rotate(transform.rotation());
    newTransform.translate(shiftShape.getX().translation());
    newTransform.rotate(shiftShape.getX().rotation());
    convertVolume(context, *shapeOp, name, newTransform, volumes);
  }
}
}  // namespace GeoModel
}  // namespace Acts
