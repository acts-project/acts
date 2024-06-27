// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisFwd.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <stdexcept>
#include <type_traits>

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

std::ostream& operator<<(std::ostream& os,
                         GridPortalLink::Direction direction) {
  switch (direction) {
    case GridPortalLink::Direction::loc0:
      os << "loc0";
      break;
    case GridPortalLink::Direction::loc1:
      os << "loc1";
      break;
  }
  return os;
}

// MARK: - PortalLinkBase

std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link) {
  link.toStream(os);
  return os;
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
    const GeometryContext& gctx, const PortalLinkBase& other,
    BinningValue direction, const Logger& logger) const {
  ACTS_DEBUG("Merging tro portals");

  ACTS_VERBOSE(" - this: " << *this);
  ACTS_VERBOSE(" - other: " << other);

  const auto& surfaceA = this->surface();
  const auto& surfaceB = other.surface();

  throw_assert(&surfaceA != &surfaceB,
               "Cannot merge portals to the same surface");

  throw_assert(surfaceA.type() == surfaceB.type(),
               "Cannot merge portals of different surface types");

  throw_assert(surfaceA.bounds().type() == surfaceB.bounds().type(),
               "Cannot merge portals of different surface bounds");

  if (const auto* cylA = dynamic_cast<const CylinderSurface*>(&surfaceA);
      cylA != nullptr) {
    const auto* cylB = dynamic_cast<const CylinderSurface*>(&surfaceB);
    throw_assert(cylB != nullptr,
                 "Cannot merge CylinderSurface with "
                 "non-CylinderSurface");
    throw_assert(
        direction == binZ || direction == binRPhi,
        "Invalid binning direction: " + binningValueNames()[direction]);

    auto mergedSurface = cylA->mergedWith(gctx, *cylB, direction);

    // mergeImpl(*cylB, other, direction, logger);
  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&surfaceA);
             discA != nullptr) {
    const auto* discB = dynamic_cast<const DiscSurface*>(&surfaceB);
    throw_assert(discB != nullptr,
                 "Cannot merge DiscSurface with non-DiscSurface");
    throw_assert(
        direction == binR || direction == binPhi,
        "Invalid binning direction: " + binningValueNames()[direction]);

    auto mergedSurface = discA->mergedWith(gctx, *discB, direction);
  } else {
    throw std::logic_error{"Surface type is not supported"};
  }

  return nullptr;
}

#if 0
std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
    const GridPortalLink1& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink1 + GridPortalLink1 with offset: "
             << offset.transpose());

  const auto& a = *this;
  const auto& b = other;

  assert(a.grid().axes().size() == 1);
  assert(b.grid().axes().size() == 1);
  ACTS_VERBOSE("Axis counts are good");

  if (a.direction() != b.direction()) {
    // 1D axes are not aligned, we cannot merge them
    throw std::logic_error{"Cannot merge 1D grids with different directions"};
  }

  const auto commonDirection = a.direction();
  ACTS_VERBOSE("Directions are consistent: " << commonDirection);

  const IAxis& axisA = *a.grid().axes().at(0);
  const IAxis& axisB = *b.grid().axes().at(0);

  const auto bdtA = axisA.getBoundaryType();
  const auto bdtB = axisB.getBoundaryType();

  if (bdtA == AxisBoundaryType::Open || bdtB == AxisBoundaryType::Open) {
    // Rejecting Open axes outright, so we don't have to handle overflow bins
    // for now
    ACTS_ERROR("Cannot merge 1D grids with Open axes");
    throw std::logic_error{"Cannot merge 1D grids with Open axes"};
  }

  ACTS_VERBOSE("Boundary types are consistent: " << bdtA << " + " << bdtB)

  ACTS_VERBOSE("- grid offset is " << offset.transpose());
  PortalDirection direction = PortalDirection::loc0;
  if (offset[0] != 0 && offset[1] != 0) {
    ACTS_ERROR("Offset is not aligned with either loc0 or loc1");
    throw std::logic_error{"Cannot merge 1D grids with diagonal offsets"};
  } else if (offset[0] != 0) {
    direction = PortalDirection::loc0;
  } else if (offset[1] != 0) {
    direction = PortalDirection::loc1;
  } else {
    ACTS_ERROR("Offset is zero");
    throw std::logic_error{"Cannot merge 1D grids with zero offset"};
  }
  ACTS_VERBOSE("=> merging along " << direction);
  ActsScalar localOffset = offset[0] != 0 ? offset[0] : offset[1];

  if (commonDirection == direction) {
    ACTS_VERBOSE("Merging along the common direction (" << commonDirection
                                                        << ") was requested");
    // Merging axes along the single binning direction, so we extend binnings
    if (bdtA != AxisBoundaryType::Bound || bdtB != AxisBoundaryType::Bound) {
      // one of the axes is not bound, cannot merge
      ACTS_ERROR("Axes are not bound, refusing to merge them");
      throw std::logic_error{
          "Cannot merge 1D grids with axes != Bound along common direction"};
    }

    ACTS_VERBOSE("Are they both equidistant and have same bin width?");
    if (axisA.isEquidistant() == axisB.isEquidistant()) {
      ActsScalar binWidthA =
          (axisA.getMax() - axisA.getMin()) / axisA.getNBins();
      ActsScalar binWidthB =
          (axisB.getMax() - axisB.getMin()) / axisB.getNBins();

      ActsScalar aMin = axisA.getMin();
      ActsScalar aMax = axisA.getMax();
      ActsScalar bMin = axisB.getMin() + localOffset;
      ActsScalar bMax = axisB.getMax() + localOffset;

      ACTS_VERBOSE("  - axis a: [" << aMin << " -> " << aMax << "] with "
                                   << axisA.getNBins() << " bins of "
                                   << binWidthA);
      ACTS_VERBOSE("  - axis b: [" << bMin << " -> " << bMax << "] with "
                                   << axisB.getNBins() << " bins of "
                                   << binWidthB);

      if (binWidthA == binWidthB) {
        ACTS_VERBOSE("  => yes!");

        ACTS_VERBOSE("Do their edges line up?");

        constexpr auto tolerance = s_onSurfaceTolerance;

        if (std::abs(aMax - bMin) > tolerance) {
          ACTS_ERROR(
              "=> no! Axes edges overlap or have gaps, refusing to merge");
          throw std::logic_error{
              "Cannot merge 1D grids with non-matching "
              "axes along common direction"};
        }

        ACTS_VERBOSE("=> yes!")

        // ActsScalar width =
        //     axisA.getMax() - axisA.getMin() + axisB.getMax() -
        //     axisB.getMin();
        ActsScalar max = axisB.getMax() + localOffset;
        ACTS_VERBOSE("New axis will be [" << aMin << " -> " << bMax << "]");

        return GridPortalLink::make(
            Axis{AxisBound{}, aMin, bMax, axisA.getNBins() + axisB.getNBins()});
      } else {
        ACTS_VERBOSE("  => no!");
      }
    }
  } else {
    ACTS_VERBOSE("Merging across the common direction (" << commonDirection
                                                         << ") was requested");
    // Merging axes across the common direction, we don't care about the
    // boundary type, as long as they're the same
    throw std::domain_error{"NotImplemented"};
  }

  return nullptr;
}

  std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
      const GridPortalLink2& other, const Vector2& offset, const Logger& logger)
      const {
    ACTS_DEBUG("Merge GridPortalLink1 + GridPortalLink2 with offset: "
               << offset.transpose());
    return nullptr;
  }

  // MARK : -GridPortalLink2

  std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
      const PortalLinkBase& other, const Vector2& offset, const Logger& logger)
      const {
    ACTS_DEBUG("Merge GridPortalLink2 + PortalLinkBase with offset: "
               << offset.transpose());
    return nullptr;
  }

  std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
      const GridPortalLink2& other, const Vector2& offset, const Logger& logger)
      const {
    ACTS_DEBUG("Merge GridPortalLink2 + GridPortalLink2 with offset: "
               << offset.transpose());
    return nullptr;
  }

  std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
      const GridPortalLink1& other, const Vector2& offset, const Logger& logger)
      const {
    ACTS_DEBUG("Merge GridPortalLink2 + GridPortalLink1 with offset: "
               << offset.transpose());
    return nullptr;
  }
#endif

}  // namespace Acts
