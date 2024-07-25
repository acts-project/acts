// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GridPortalLink.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

namespace Acts {

std::unique_ptr<GridPortalLink> GridPortalLink::make(
    const std::shared_ptr<RegularSurface>& surface,
    const TrackingVolume& volume, BinningValue direction) {
  std::unique_ptr<GridPortalLink> grid;

  if (const auto* cylinder =
          dynamic_cast<const CylinderSurface*>(surface.get());
      cylinder != nullptr) {
    if (direction == BinningValue::binRPhi) {
      ActsScalar r = cylinder->bounds().get(CylinderBounds::eR);
      if (cylinder->bounds().coversFullAzimuth()) {
        grid = GridPortalLink::make(surface, direction,
                                    Axis{AxisClosed, -M_PI * r, M_PI * r, 1});
      } else {
        ActsScalar hlPhi =
            cylinder->bounds().get(CylinderBounds::eHalfPhiSector);

        grid = GridPortalLink::make(surface, direction,
                                    Axis{AxisBound, -hlPhi * r, hlPhi * r, 1});
      }
    } else if (direction == BinningValue::binZ) {
      ActsScalar hlZ = cylinder->bounds().get(CylinderBounds::eHalfLengthZ);
      grid = GridPortalLink::make(surface, direction,
                                  Axis{AxisBound, -hlZ, hlZ, 1});
    } else {
      throw std::invalid_argument{"Invalid binning direction"};
    }
  } else if (const auto* disc = dynamic_cast<const DiscSurface*>(surface.get());
             disc != nullptr) {
    const auto& bounds = dynamic_cast<const RadialBounds&>(disc->bounds());
    if (direction == BinningValue::binR) {
      ActsScalar minR = bounds.get(RadialBounds::eMinR);
      ActsScalar maxR = bounds.get(RadialBounds::eMaxR);
      grid = GridPortalLink::make(surface, direction,
                                  Axis{AxisBound, minR, maxR, 1});
    } else if (direction == BinningValue::binPhi) {
      if (bounds.coversFullAzimuth()) {
        grid = GridPortalLink::make(surface, direction,
                                    Axis{AxisClosed, -M_PI, M_PI, 1});
      } else {
        ActsScalar hlPhi = bounds.get(RadialBounds::eHalfPhiSector);
        grid = GridPortalLink::make(surface, direction,
                                    Axis{AxisBound, -hlPhi, hlPhi, 1});
      }
    } else {
      throw std::invalid_argument{"Invalid binning direction"};
    }
  } else if (const auto* plane =
                 dynamic_cast<const PlaneSurface*>(surface.get());
             plane != nullptr) {
    throw std::invalid_argument{"Plane surface is not implemented yet"};

  } else {
    throw std::invalid_argument{"Surface type is not supported"};
  }

  assert(grid != nullptr);
  grid->setVolume(&volume);

  return grid;
}

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

void GridPortalLink::printContents(std::ostream& os) const {
  std::size_t dim = grid().axes().size();
  os << "----- GRID " << dim << "d -----" << std::endl;
  os << grid() << " along " << direction() << std::endl;

  std::string loc0;
  std::string loc1;

  bool flipped = false;

  if (surface().type() == Surface::Cylinder) {
    loc0 = "rPhi";
    loc1 = "z";
    flipped = direction() != BinningValue::binRPhi;
  } else if (surface().type() == Surface::Disc) {
    loc0 = "r";
    loc1 = "phi";
    flipped = direction() != BinningValue::binR;
  } else if (surface().type() == Surface::Plane) {
    loc0 = "x";
    loc1 = "y";
    flipped = direction() != BinningValue::binX;
  } else {
    throw std::invalid_argument{"Unsupported surface type"};
  }

  if (dim == 1) {
    auto loc = numLocalBins();

    if (flipped) {
      os << std::format("{} >       i=0 ", loc1);
      for (std::size_t i = 1; i <= loc.at(0) + 1; i++) {
        os << std::format("{:>11} ", std::format("i={}", i));
      }
      os << std::endl;

      os << "    ";
      for (std::size_t i = 0; i <= loc.at(0) + 1; i++) {
        const void* v = atLocalBins({i});
        os << std::format("{:11}", v) << " ";
      }
      os << std::endl;

    } else {
      os << std::format("v {}", loc0) << std::endl;
      for (std::size_t i = 0; i <= loc.at(0) + 1; i++) {
        os << "i=" << i << " ";
        const void* v = atLocalBins({i});
        os << std::format("{:11}", v) << " ";
        os << std::endl;
      }
    }

  } else {
    auto loc = numLocalBins();
    os << std::format("v {}|{} >   j=0 ", loc0, loc1, 0);
    for (std::size_t j = 1; j <= loc.at(1) + 1; j++) {
      os << std::format("{:>11} ", std::format("j={}", j));
    }
    os << std::endl;
    for (std::size_t i = 0; i <= loc.at(0) + 1; i++) {
      os << "i=" << i << " ";
      for (std::size_t j = 0; j <= loc.at(1) + 1; j++) {
        const void* v = atLocalBins({i, j});
        os << std::format("{:11}", v) << " ";
      }
      os << std::endl;
    }
  }
}

}  // namespace Acts
