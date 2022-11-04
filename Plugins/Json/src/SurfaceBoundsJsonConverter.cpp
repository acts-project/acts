// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/SurfaceBoundsJsonConverter.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::SurfaceBounds& bounds) {
  j["type"] = boundTypes[bounds.type()];
  j["values"] = bounds.values();
}

Acts::SurfaceBounds::BoundsType Acts::boundsType(
    const std::string& bTypeString) {
  auto fType = std::find(boundTypes.begin(), boundTypes.end(), bTypeString);
  if (fType != boundTypes.end()) {
    return SurfaceBounds::BoundsType(std::distance(boundTypes.begin(), fType));
  }
  return SurfaceBounds::BoundsType::eOther;
}

std::shared_ptr<Acts::SurfaceBounds> Acts::createBounds(
    const SurfaceBounds::BoundsType& bType,
    const std::vector<Acts::ActsScalar>& bVector) {
  std::shared_ptr<SurfaceBounds> bounds = nullptr;
  // One could write this with a tuple iteration
  /// @TODO check that
  switch (bType) {
    case SurfaceBounds::BoundsType::eCone: {
      bounds = createBounds<Acts::ConeBounds>(bVector);
    } break;
    case SurfaceBounds::BoundsType::eCylinder: {
      bounds = createBounds<Acts::CylinderBounds>(bVector);
    } break; /*
     case SurfaceBounds::BoundsType::eDiamond: {
       return std::make_shared<Acts::DiamondBounds>();
     } break;
     */
    case SurfaceBounds::BoundsType::eDisc: {
      bounds = createBounds<Acts::RadialBounds>(bVector);
    } break;
    case SurfaceBounds::BoundsType::eEllipse: {
      bounds = createBounds<Acts::EllipseBounds>(bVector);
    } break;
    case SurfaceBounds::BoundsType::eLine: {
      bounds = createBounds<Acts::LineBounds>(bVector);
    } break;
    case SurfaceBounds::BoundsType::eRectangle: {
      bounds = createBounds<Acts::RectangleBounds>(bVector);
    } break;
    case SurfaceBounds::BoundsType::eTrapezoid: {
      bounds = createBounds<Acts::TrapezoidBounds>(bVector);
    } break;
    /*
    case SurfaceBounds::BoundsType::eTriangle: {
      return std::make_shared<Acts::TriangleBounds>();
    } break;
    */
    case SurfaceBounds::BoundsType::eDiscTrapezoid: {
      bounds = createBounds<Acts::DiscTrapezoidBounds>(bVector);
    } break;
    /*
    case SurfaceBounds::BoundsType::eConvexPolygon: {
      return std::make_shared<Acts::ConvexPolygonBounds>();
    } break;
    */
    case SurfaceBounds::BoundsType::eAnnulus: {
      bounds = createBounds<Acts::AnnulusBounds>(bVector);

    } break;
    default:
      break;
  }
  return bounds;
}
