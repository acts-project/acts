// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

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

std::ostream& operator<<(std::ostream& os, PortalDirection direction) {
  switch (direction) {
    case PortalDirection::loc0:
      os << "loc0";
      break;
    case PortalDirection::loc1:
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

std::unique_ptr<GridPortalLink1> PortalLinkBase::merge1d(
    const GridPortalLink1& a, const GridPortalLink1& b, const Vector2& offset,
    const Logger& logger) {
  return nullptr;
}

// MARK: - GridPortalLink1

std::unique_ptr<PortalLinkBase> GridPortalLink1::merge(
    const PortalLinkBase& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink1 + PortalLinkBase > " << offset.transpose());
  return other.merge(*this, offset, logger);
}

// std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
//     const PortalLinkBase& other, const Vector2& offset,
//     const Logger& logger) const {
// ACTS_DEBUG("Merge PortalLinkBase + PortalLinkBase with offset: "
//            << offset.transpose());
//
// if (const auto* a1 = dynamic_cast<const GridPortalLink1*>(this);
//     a1 != nullptr) {
//   if (const auto* b1 = dynamic_cast<const GridPortalLink1*>(&other);
//       b1 != nullptr) {
//     ACTS_DEBUG("Merging GridPortalLink1 + GridPortalLink1");
//     merge1d(*a1, *b1, offset, logger);
//
//   } else if (const auto* b2 = dynamic_cast<const GridPortalLink2*>(&other);
//              b2 != nullptr) {
//     ACTS_DEBUG("Merging GridPortalLink1 + GridPortalLink2");
//
//   } else {
//     ACTS_DEBUG("Merging GridPortalLink1 + unknown");
//   }
// } else if (const auto* a2 = dynamic_cast<const GridPortalLink2*>(this);
//            a2 != nullptr) {
//   if (const auto* b1 = dynamic_cast<const GridPortalLink1*>(&other);
//       b1 != nullptr) {
//     ACTS_DEBUG("Merging GridPortalLink2 + GridPortalLink1");
//
//   } else if (const auto* b2 = dynamic_cast<const GridPortalLink2*>(&other);
//              b2 != nullptr) {
//     ACTS_DEBUG("Merging GridPortalLink2 + GridPortalLink2");
//
//   } else {
//     ACTS_DEBUG("Merging GridPortalLink2 + unknown");
//   }
// } else {
//   ACTS_DEBUG("Merging unknown + unknown");
// }
//
//   return nullptr;
// }

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
    const GridPortalLink2& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink1 + GridPortalLink2 with offset: "
             << offset.transpose());
  return nullptr;
}

// MARK : -GridPortalLink2

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const PortalLinkBase& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink2 + PortalLinkBase with offset: "
             << offset.transpose());
  return nullptr;
}

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const GridPortalLink2& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink2 + GridPortalLink2 with offset: "
             << offset.transpose());
  return nullptr;
}

std::unique_ptr<PortalLinkBase> GridPortalLink2::merge(
    const GridPortalLink1& other, const Vector2& offset,
    const Logger& logger) const {
  ACTS_DEBUG("Merge GridPortalLink2 + GridPortalLink1 with offset: "
             << offset.transpose());
  return nullptr;
}

}  // namespace Acts
