// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <stdexcept>

namespace Acts {

Portal::Portal(std::shared_ptr<RegularSurface> surface)
    : m_surface(std::move(surface)) {
  throw_assert(m_surface, "Portal surface is nullptr");
}

const Volume* Portal::resolveVolume(const GeometryContext& gctx,
                                    const Vector3& position,
                                    const Vector3& direction) const {
  const Vector3 normal = m_surface->normal(gctx, position);
  Direction side = Direction::fromScalar(normal.dot(direction));

  const std::unique_ptr<PortalLinkBase>& link =
      side == Direction::AlongNormal ? m_alongNormal : m_oppositeNormal;

  return nullptr;
  if (link == nullptr) {
    // no link is attached in this direction => this is the end of the world as
    // we know it. (i feel fine)
    return nullptr;
  } else {
    // return link->resolveVolume(position);
  }
}

// MARK: - PortalLinkBase

std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link) {
  link.toStream(os);
  return os;
}

// MARK: - GridPortalLinks

void GridPortalLink::checkConsistency(const CylinderSurface& cyl) const {
  if (cyl.bounds().get(CylinderBounds::eAveragePhi) != 0) {
    throw std::invalid_argument(
        "GridPortalLink: CylinderBounds: only average phi == 0 is "
        "supported. Rotate the cylinder surface.");
  };

  constexpr auto tolerance = s_onSurfaceTolerance;
  auto same = [](auto a, auto b) { return std::abs(a - b) < tolerance; };

  auto checkZ = [&](const IAxis& axis) {
    ActsScalar hlZ = cyl.bounds().get(CylinderBounds::eHalfLengthZ);
    if (!same(axis.getMin(), -hlZ) || !same(axis.getMax(), hlZ)) {
      throw std::invalid_argument(
          "GridPortalLink: CylinderBounds: invalid length setup.");
    }
  };
  auto checkRPhi = [&](const IAxis& axis) {
    ActsScalar hlPhi = cyl.bounds().get(CylinderBounds::eHalfPhiSector);
    ActsScalar r = cyl.bounds().get(CylinderBounds::eR);
    ActsScalar hlRPhi = r * hlPhi;

    if (!same(axis.getMin(), -hlRPhi) || !same(axis.getMax(), hlRPhi)) {
      throw std::invalid_argument(
          "GridPortalLink: CylinderBounds: invalid phi sector setup: axes "
          "don't match bounds");
    }

    // If full cylinder, make sure axis wraps around
    if (same(hlPhi, M_PI)) {
      if (axis.getBoundaryType() != AxisBoundaryType::Closed) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid phi sector setup: "
            "axis is not closed.");
      }
    } else {
      if (axis.getBoundaryType() != AxisBoundaryType::Bound) {
        throw std::invalid_argument(
            "GridPortalLink: CylinderBounds: invalid phi sector setup: "
            "axis is not bound.");
      }
    }
  };

  if (dim() == 1) {
    const IAxis& axisLoc0 = *grid().axes().front();
    if (direction() == BinningValue::binRPhi) {
      checkRPhi(axisLoc0);
    } else {
      checkZ(axisLoc0);
    }
  } else {  // DIM == 2
    const auto& axisLoc0 = *grid().axes().front();
    const auto& axisLoc1 = *grid().axes().back();
    checkRPhi(axisLoc0);
    checkZ(axisLoc1);
  }
}

void GridPortalLink::checkConsistency(const DiscSurface& disc) const {
  constexpr auto tolerance = s_onSurfaceTolerance;
  auto same = [](auto a, auto b) { return std::abs(a - b) < tolerance; };

  const auto* bounds = dynamic_cast<const RadialBounds*>(&disc.bounds());
  if (bounds == nullptr) {
    throw std::invalid_argument(
        "GridPortalLink: DiscBounds: invalid bounds type.");
  }

  if (bounds->get(RadialBounds::eAveragePhi) != 0) {
    throw std::invalid_argument(
        "GridPortalLink: DiscBounds: only average phi == 0 is supported. "
        "Rotate the disc surface.");
  }

  auto checkR = [&](const IAxis& axis) {
    ActsScalar minR = bounds->get(RadialBounds::eMinR);
    ActsScalar maxR = bounds->get(RadialBounds::eMaxR);
    if (!same(axis.getMin(), minR) || !same(axis.getMax(), maxR)) {
      throw std::invalid_argument(
          "GridPortalLink: DiscBounds: invalid radius setup.");
    }
  };

  auto checkPhi = [&](const IAxis& axis) {
    ActsScalar hlPhi = bounds->get(RadialBounds::eHalfPhiSector);
    if (!same(axis.getMin(), -hlPhi) || !same(axis.getMax(), hlPhi)) {
      throw std::invalid_argument(
          "GridPortalLink: DiscBounds: invalid phi sector setup.");
    }
    // If full disc, make sure axis wraps around
    if (same(hlPhi, M_PI)) {
      if (axis.getBoundaryType() != Acts::AxisBoundaryType::Closed) {
        throw std::invalid_argument(
            "GridPortalLink: DiscBounds: invalid phi sector setup: axis is "
            "not closed.");
      }
    } else {
      if (axis.getBoundaryType() != Acts::AxisBoundaryType::Bound) {
        throw std::invalid_argument(
            "GridPortalLink: DiscBounds: invalid phi sector setup: axis "
            "is not bound.");
      }
    }
  };

  if (dim() == 1) {
    const IAxis& axisLoc0 = *grid().axes().front();
    if (direction() == BinningValue::binR) {
      checkR(axisLoc0);
    } else {
      checkPhi(axisLoc0);
    }
  } else {  // DIM == 2
    const auto& axisLoc0 = *grid().axes().front();
    const auto& axisLoc1 = *grid().axes().back();
    checkR(axisLoc0);
    checkPhi(axisLoc1);
  }
}

std::unique_ptr<GridPortalLink> GridPortalLink::extendTo2D(
    const std::shared_ptr<CylinderSurface>& surface) const {
  assert(dim() == 1);
  if (direction() == BinningValue::binRPhi) {
    const auto& axisRPhi = *grid().axes().front();
    // 1D direction is binRPhi, so add a Z axis
    ActsScalar hlZ = surface->bounds().get(CylinderBounds::eHalfLengthZ);

    Axis axisZ{AxisBound, -hlZ, hlZ, 1};
    return axisRPhi.visit([&](const auto& axis) {
      return GridPortalLink::make(surface, Axis{axis}, std::move(axisZ));
    });
  } else {
    const auto& axisZ = *grid().axes().front();
    // 1D direction is binZ, so add an rPhi axis
    ActsScalar r = surface->bounds().get(CylinderBounds::eR);
    ActsScalar hlPhi = surface->bounds().get(CylinderBounds::eHalfPhiSector);
    ActsScalar hlRPhi = r * hlPhi;

    auto axis = [&](auto bdt) {
      return axisZ.visit([&](const auto& concreteAxis) {
        return GridPortalLink::make(surface, Axis{bdt, -hlRPhi, hlRPhi, 1},
                                    Axis{concreteAxis});
      });
    };

    if (surface->bounds().coversFullAzimuth()) {
      return axis(AxisClosed);
    } else {
      return axis(AxisBound);
    }
  }
}

std::unique_ptr<GridPortalLink> GridPortalLink::extendTo2D(
    const std::shared_ptr<DiscSurface>& surface) const {
  assert(dim() == 1);

  const auto* bounds = dynamic_cast<const RadialBounds*>(&surface->bounds());
  if (bounds == nullptr) {
    throw std::invalid_argument(
        "GridPortalLink: DiscBounds: invalid bounds type.");
  }

  if (direction() == BinningValue::binR) {
    const auto& axisR = *grid().axes().front();
    // 1D direction is binR, so add a phi axis
    ActsScalar hlPhi = bounds->get(RadialBounds::eHalfPhiSector);

    auto axis = [&](auto bdt) {
      return axisR.visit([&](const auto& concreteAxis) {
        return GridPortalLink::make(surface, Axis{concreteAxis},
                                    Axis{bdt, -hlPhi, hlPhi, 1});
      });
    };

    if (bounds->coversFullAzimuth()) {
      return axis(AxisClosed);
    } else {
      return axis(AxisBound);
    }
  } else {
    const auto& axisPhi = *grid().axes().front();
    // 1D direction is binPhi, so add an R axis
    ActsScalar rMin = bounds->get(RadialBounds::eMinR);
    ActsScalar rMax = bounds->get(RadialBounds::eMaxR);

    return axisPhi.visit([&](const auto& axis) {
      return GridPortalLink::make(surface, Axis{AxisBound, rMin, rMax, 1},
                                  Axis{axis});
    });
  }
}

}  // namespace Acts
