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
namespace GeoModelToDetVol {
std::shared_ptr<Experimental::DetectorVolume> convert(
    const GeometryContext& context, const GeoShape& shape,
    const std::string& name, const Transform3& transform) {
  auto portalGenerator = Experimental::defaultPortalAndSubPortalGenerator();
  if (shape.typeID() == GeoTube::getClassTypeID()) {
    const GeoTube& tube = dynamic_cast<const GeoTube&>(shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tube.getRMin(), tube.getRMax(),
                                               tube.getZHalfLength());
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortals());
  } else if (shape.typeID() == GeoTubs::getClassTypeID()) {
    const GeoTubs& tubs = dynamic_cast<const GeoTubs&>(shape);
    std::shared_ptr<CylinderVolumeBounds> bounds =
        std::make_shared<CylinderVolumeBounds>(tubs.getRMin(), tubs.getRMax(),
                                               tubs.getZHalfLength(),
                                               tubs.getDPhi() / 2);
    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    newTransform.translate(transform.translation());
    newTransform.rotate(Eigen::AngleAxisd(tubs.getSPhi() + 0.5 * tubs.getDPhi(),
                                          Acts::Vector3::UnitZ()));
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, newTransform, bounds,
        Experimental::tryAllPortals());
  } else if (shape.typeID() == GeoBox::getClassTypeID()) {
    const GeoBox& box = dynamic_cast<const GeoBox&>(shape);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(
            box.getXHalfLength(), box.getYHalfLength(), box.getZHalfLength());
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortals());
  } else if (shape.typeID() == GeoSimplePolygonBrep::getClassTypeID()) {
    // Will change this in the future
    double xmin{0};
    double xmax{0};
    double ymin{0};
    double ymax{0};
    double zmin{0};
    double zmax{0};
    const GeoSimplePolygonBrep& brep =
        dynamic_cast<const GeoSimplePolygonBrep&>(shape);
    brep.extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(xmax - xmin, ymax - ymin,
                                             zmax - zmin);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortals());
  } else if (shape.typeID() == GeoTrd::getClassTypeID()) {
    const GeoTrd& trd = dynamic_cast<const GeoTrd&>(shape);
    float x1 = trd.getXHalfLength1();
    float x2 = trd.getXHalfLength2();
    float y1 = trd.getYHalfLength1();
    float y2 = trd.getYHalfLength2();
    float z = trd.getZHalfLength();
    if (y1 == y2) {
      if (x1 <= x2) {
        // y axis in ACTS is z axis in geomodel
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x1, x2, z, y1);
        auto rotationAngle = M_PI / 2;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()));
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortals());
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(x2, x1, z, y1);
        auto rotationAngle = M_PI;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitY()) *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitZ()));
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortals());
      }
    } else if (x1 == x2) {
      if (y1 < y2) {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y1, y2, z, x1);
        auto rotationAngle = M_PI / 2;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitZ()) *
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()));
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortals());
      } else {
        std::shared_ptr<TrapezoidVolumeBounds> bounds =
            std::make_shared<TrapezoidVolumeBounds>(y2, y1, z, x1);
        auto rotationAngle = M_PI;
        Acts::Transform3 newTransform = Acts::Transform3::Identity();
        newTransform.translate(transform.translation());
        newTransform.rotate(
            Eigen::AngleAxisd(rotationAngle, Acts::Vector3::UnitX()) *
            Eigen::AngleAxisd(rotationAngle / 2, Acts::Vector3::UnitZ()) *
            Eigen::AngleAxisd(rotationAngle / 2, Acts::Vector3::UnitX()));
        return Experimental::DetectorVolumeFactory::construct(
            portalGenerator, context, name, newTransform, bounds,
            Experimental::tryAllPortals());
      }
    } else {
      throw std::runtime_error("FATAL: Translating GeoTrd to ACTS failed");
    }
  } else if (shape.typeID() == GeoShapeUnion::getClassTypeID()) {
    // Get the bounding box of the union
    double xmin{0};
    double xmax{0};
    double ymin{0};
    double ymax{0};
    double zmin{0};
    double zmax{0};
    const GeoShapeUnion& unionShape = dynamic_cast<const GeoShapeUnion&>(shape);
    unionShape.extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(xmax - xmin, ymax - ymin,
                                             zmax - zmin);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortals());
  } else if (shape.typeID() == GeoShapeSubtraction::getClassTypeID()) {
    // Go down the left side (opA) of the subtraction until we reach a normal
    // shape
    const GeoShapeSubtraction& subtractionShape =
        dynamic_cast<const GeoShapeSubtraction&>(shape);
    const GeoShape* shapeA = subtractionShape.getOpA();
    return convert(context, *shapeA, name, transform);
  } else if (shape.typeID() == GeoPcon::getClassTypeID()) {
    // Will change in future, get bounding box for now
    double xmin{0};
    double xmax{0};
    double ymin{0};
    double ymax{0};
    double zmin{0};
    double zmax{0};
    const GeoPcon& pcon = dynamic_cast<const GeoPcon&>(shape);
    pcon.extent(xmin, ymin, zmin, xmax, ymax, zmax);
    std::shared_ptr<CuboidVolumeBounds> bounds =
        std::make_shared<CuboidVolumeBounds>(xmax - xmin, ymax - ymin,
                                             zmax - zmin);
    return Experimental::DetectorVolumeFactory::construct(
        portalGenerator, context, name, transform, bounds,
        Experimental::tryAllPortals());
  }
  if (shape.typeID() == GeoShapeShift::getClassTypeID()) {
    const GeoShapeShift& shiftShape = dynamic_cast<const GeoShapeShift&>(shape);
    const GeoShape* shapeOp = shiftShape.getOp();
    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    newTransform.translate(transform.translation());
    newTransform.rotate(transform.rotation());
    newTransform.translate(shiftShape.getX().translation());
    return convert(context, *shapeOp, name, newTransform);
  }
  throw std::runtime_error("Unknown shape type: " + shape.type());
}
}  // namespace GeoModelToDetVol
}  // namespace Acts
