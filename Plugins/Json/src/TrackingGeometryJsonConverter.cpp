// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/TrackingGeometryJsonConverter.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Plugins/Json/GeometryJsonKeys.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include <Acts/Surfaces/AnnulusBounds.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>

namespace {

using volumePointer = const Acts::TrackingVolume*;
using surfacePointer = const Acts::Surface*;

Acts::BinUtility defaultSurfaceBin(const Acts::Surface* surface) {
  Acts::BinUtility bUtility;

  // Check which type of bounds is associated to the surface
  const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
  const Acts::RadialBounds* radialBounds =
      dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
  const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
  const Acts::AnnulusBounds* annulusBounds =
      dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);
  const Acts::RectangleBounds* rectangleBounds =
      dynamic_cast<const Acts::RectangleBounds*>(&surfaceBounds);

  if (radialBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1,
        radialBounds->get(Acts::RadialBounds::eAveragePhi) -
            radialBounds->get(Acts::RadialBounds::eHalfPhiSector),
        radialBounds->get(Acts::RadialBounds::eAveragePhi) +
            radialBounds->get(Acts::RadialBounds::eHalfPhiSector),
        (radialBounds->get(Acts::RadialBounds::eHalfPhiSector) - M_PI) <
                Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility += Acts::BinUtility(1, radialBounds->rMin(), radialBounds->rMax(),
                                 Acts::open, Acts::binR);
    return bUtility;
  }
  if (cylinderBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1,
        cylinderBounds->get(Acts::CylinderBounds::eAveragePhi) -
            cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector),
        cylinderBounds->get(Acts::CylinderBounds::eAveragePhi) +
            cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector),
        (cylinderBounds->get(Acts::CylinderBounds::eHalfPhiSector) - M_PI) <
                Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility += Acts::BinUtility(
        1, -1 * cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
        cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ), Acts::open,
        Acts::binZ);
    return bUtility;
  }
  if (annulusBounds != nullptr) {
    bUtility +=
        Acts::BinUtility(1, annulusBounds->get(Acts::AnnulusBounds::eMinPhiRel),
                         annulusBounds->get(Acts::AnnulusBounds::eMaxPhiRel),
                         Acts::open, Acts::binPhi);
    bUtility += Acts::BinUtility(1, annulusBounds->rMin(),
                                 annulusBounds->rMax(), Acts::open, Acts::binR);
    return bUtility;
  }
  if (rectangleBounds != nullptr) {
    bUtility +=
        Acts::BinUtility(1, rectangleBounds->get(Acts::RectangleBounds::eMinX),
                         rectangleBounds->get(Acts::RectangleBounds::eMaxX),
                         Acts::open, Acts::binX);
    bUtility +=
        Acts::BinUtility(1, rectangleBounds->get(Acts::RectangleBounds::eMinY),
                         rectangleBounds->get(Acts::RectangleBounds::eMaxY),
                         Acts::open, Acts::binY);
    return bUtility;
  }
  return bUtility;
}

Acts::BinUtility defaultVolumeBin(const Acts::TrackingVolume* volume) {
  Acts::BinUtility bUtility;

  // Check which type of bound is associated to the volume
  auto cyBounds = dynamic_cast<const Acts::CylinderVolumeBounds*>(
      &(volume->volumeBounds()));
  auto cutcylBounds = dynamic_cast<const Acts::CutoutCylinderVolumeBounds*>(
      &(volume->volumeBounds()));
  auto cuBounds =
      dynamic_cast<const Acts::CuboidVolumeBounds*>(&(volume->volumeBounds()));

  if (cyBounds != nullptr) {
    bUtility +=
        Acts::BinUtility(1, cyBounds->get(Acts::CylinderVolumeBounds::eMinR),
                         cyBounds->get(Acts::CylinderVolumeBounds::eMaxR),
                         Acts::open, Acts::binR);
    bUtility += Acts::BinUtility(
        1, -cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector),
        cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector),
        (cyBounds->get(Acts::CylinderVolumeBounds::eHalfPhiSector) - M_PI) <
                Acts::s_epsilon
            ? Acts::closed
            : Acts::open,
        Acts::binPhi);
    bUtility += Acts::BinUtility(
        1, -cyBounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ),
        cyBounds->get(Acts::CylinderVolumeBounds::eHalfLengthZ), Acts::open,
        Acts::binZ);
    return bUtility;
  }
  if (cutcylBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1, cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eMinR),
        cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eMaxR), Acts::open,
        Acts::binR);
    bUtility += Acts::BinUtility(1, -M_PI, M_PI, Acts::closed, Acts::binPhi);
    bUtility += Acts::BinUtility(
        1, -cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eHalfLengthZ),
        cutcylBounds->get(Acts::CutoutCylinderVolumeBounds::eHalfLengthZ),
        Acts::open, Acts::binZ);
    return bUtility;
  } else if (cuBounds != nullptr) {
    bUtility += Acts::BinUtility(
        1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthX),
        cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthX), Acts::open,
        Acts::binX);
    bUtility += Acts::BinUtility(
        1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthY),
        cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthY), Acts::open,
        Acts::binY);
    bUtility += Acts::BinUtility(
        1, -cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthZ),
        cuBounds->get(Acts::CuboidVolumeBounds::eHalfLengthZ), Acts::open,
        Acts::binZ);
    return bUtility;
  }
  return bUtility;
}

}  // namespace

void Acts::to_json(nlohmann::json& j, const surfacePointer& t) {
  // Get the ID of the surface (redundant but help readability)
  std::ostringstream SurfaceID;
  SurfaceID << t->geometryId();
  j[Acts::jsonKey().geometryidkey] = SurfaceID.str();
  if (t->surfaceMaterial() == nullptr) {
    // Cast the surface bound to both disk and cylinder
    const Acts::SurfaceBounds& surfaceBounds = t->bounds();
    auto sTransform = t->transform(Acts::GeometryContext());
    j[Acts::jsonKey().typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    j[Acts::jsonKey().mapkey] = false;
    const Acts::BinUtility bUtility = defaultSurfaceBin(t);
    nlohmann::json jBin(bUtility);
    j[Acts::jsonKey().binkey] = jBin;

    const Acts::RadialBounds* radialBounds =
        dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
    if (radialBounds != nullptr) {
      j[Acts::jsonKey().surfacetypekey] = "Disk";
      j[Acts::jsonKey().surfacepositionkey] = sTransform.translation().z();
      j[Acts::jsonKey().surfacerangekey] = {radialBounds->rMin(),
                                            radialBounds->rMax()};
      return;
    }
    const Acts::CylinderBounds* cylinderBounds =
        dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
    if (cylinderBounds != nullptr) {
      j[Acts::jsonKey().surfacetypekey] = "Cylinder";
      j[Acts::jsonKey().surfacepositionkey] =
          cylinderBounds->get(Acts::CylinderBounds::eR);
      j[Acts::jsonKey().surfacerangekey] = {
          -1 * cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
          cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ)};
      return;
    }
    const Acts::AnnulusBounds* annulusBounds =
        dynamic_cast<const Acts::AnnulusBounds*>(&surfaceBounds);
    if (annulusBounds != nullptr) {
      j[Acts::jsonKey().surfacetypekey] = "Annulus";
      j[Acts::jsonKey().surfacepositionkey] = sTransform.translation().z();
      j[Acts::jsonKey().surfacerangekey] = {
          {annulusBounds->rMin(), annulusBounds->rMax()},
          {annulusBounds->phiMin(), annulusBounds->phiMax()}};
      return;
    }
  } else {
    to_json(j, t->surfaceMaterial());
  }
}

// Not intended for use, here for complitness sake
// to_json on surface is used as input for the material mapping
void Acts::from_json(const nlohmann::json& /*j*/, surfacePointer& /*t*/) {
  return;
}

void Acts::to_json(nlohmann::json& j, const volumePointer& t) {
  std::ostringstream volumeID;
  volumeID << t->geometryId();
  j[Acts::jsonKey().geometryidkey] = volumeID.str();
  j[Acts::jsonKey().namekey] = t->volumeName();
  if (t->volumeMaterial() == nullptr) {
    j[Acts::jsonKey().typekey] = "proto";
    // by default the protoMaterial is not used for mapping
    j[Acts::jsonKey().mapkey] = false;
    const Acts::BinUtility bUtility = defaultVolumeBin(t);
    nlohmann::json jBin(bUtility);
    j[Acts::jsonKey().binkey] = jBin;
  } else {
    to_json(j, t->volumeMaterial());
  }
  return;
}

// Not intended for use, here for complitness sake
// to_json on volume is used as input for the material mapping
void Acts::from_json(const nlohmann::json& /*j*/, volumePointer& /*t*/) {
  return;
}
