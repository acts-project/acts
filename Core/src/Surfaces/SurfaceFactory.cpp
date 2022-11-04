// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceFactory.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

std::shared_ptr<Acts::Surface> Acts::SurfaceFactory::createSurface(
    const Transform3& transform, Surface::SurfaceType sType,
    SurfaceBounds::BoundsType bType, const std::vector<ActsScalar>& bValues) {
  std::shared_ptr<Surface> surface = nullptr;

  switch (sType) {
    case Surface::SurfaceType::Plane: {
      switch (bType) {
        case SurfaceBounds::BoundsType::eRectangle: {
          auto rBounds = createBounds<RectangleBounds>(bValues);
          surface = Surface::makeShared<PlaneSurface>(transform, rBounds);
        } break;
        case SurfaceBounds::BoundsType::eTrapezoid: {
          auto tBounds = createBounds<TrapezoidBounds>(bValues);
          surface = Surface::makeShared<PlaneSurface>(transform, tBounds);
        } break;
        case SurfaceBounds::BoundsType::eEllipse: {
          auto eBounds = createBounds<TrapezoidBounds>(bValues);
          surface = Surface::makeShared<PlaneSurface>(transform, eBounds);
        } break;
        case SurfaceBounds::BoundsType::eDiamond: {
          auto dBounds = createBounds<TrapezoidBounds>(bValues);
          surface = Surface::makeShared<PlaneSurface>(transform, dBounds);
        } break;
        case SurfaceBounds::BoundsType::eTriangle: {
          auto dBounds = createBounds<ConvexPolygonBounds<3u>>(bValues);
          surface = Surface::makeShared<PlaneSurface>(transform, dBounds);
        } break;
        default:
          break;  // todo checking and exception throwing
      };
    } break;
    case Surface::SurfaceType::Cone: {
      auto cBounds = createBounds<ConeBounds>(bValues);
      surface = Surface::makeShared<ConeSurface>(transform, cBounds);
    } break;
    case Surface::SurfaceType::Cylinder: {
      auto cBounds = createBounds<CylinderBounds>(bValues);
      surface = Surface::makeShared<CylinderSurface>(transform, cBounds);
    } break;
    case Surface::SurfaceType::Disc: {
      switch (bType) {
        case SurfaceBounds::BoundsType::eAnnulus: {
          auto aBounds = createBounds<AnnulusBounds>(bValues);
          surface = Surface::makeShared<DiscSurface>(transform, aBounds);
        } break;
        case SurfaceBounds::BoundsType::eDisc: {
          auto rBounds = createBounds<RadialBounds>(bValues);
          surface = Surface::makeShared<DiscSurface>(transform, rBounds);
        } break;
        case SurfaceBounds::BoundsType::eDiscTrapezoid: {
          auto dBounds = createBounds<DiscTrapezoidBounds>(bValues);
          surface = Surface::makeShared<DiscSurface>(transform, dBounds);
        } break;
        default:
          break;  // todo checking and exception throwing
      };
    } break;
    case Surface::SurfaceType::Straw: {
      auto lBounds = createBounds<LineBounds>(bValues);
      surface = Surface::makeShared<StrawSurface>(transform, lBounds);
    } break;
    case Surface::SurfaceType::Perigee: {
      surface = Surface::makeShared<PerigeeSurface>(transform.translation());
    } break;
    default:
      break;
  }
  return surface;
}

std::shared_ptr<Acts::SurfaceBounds> Acts::SurfaceFactory::createBounds(
    SurfaceBounds::BoundsType bType, const std::vector<ActsScalar>& bValues) {
  std::shared_ptr<SurfaceBounds> bounds = nullptr;
  // One could write this with a tuple iteration
  switch (bType) {
    case SurfaceBounds::BoundsType::eCone: {
      bounds = createBounds<Acts::ConeBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eCylinder: {
      bounds = createBounds<Acts::CylinderBounds>(bValues);
    } break; 
    case SurfaceBounds::BoundsType::eDisc: {
      bounds = createBounds<Acts::RadialBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eEllipse: {
      bounds = createBounds<Acts::EllipseBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eLine: {
      bounds = createBounds<Acts::LineBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eRectangle: {
      bounds = createBounds<Acts::RectangleBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eTrapezoid: {
      bounds = createBounds<Acts::TrapezoidBounds>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eTriangle: {
      bounds = createBounds<Acts::ConvexPolygonBounds<3u>>(bValues);
    } break;
    case SurfaceBounds::BoundsType::eDiscTrapezoid: {
      bounds = createBounds<Acts::DiscTrapezoidBounds>(bValues);
    } break;
    /*
    case SurfaceBounds::BoundsType::eConvexPolygon: {
      return std::make_shared<Acts::ConvexPolygonBounds>();
    } break;
    */
    case SurfaceBounds::BoundsType::eAnnulus: {
      bounds = createBounds<Acts::AnnulusBounds>(bValues);

    } break;
    default:
      break;
  }
  return bounds;
}