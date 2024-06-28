// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisFwd.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
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

// MARK: - PortalLinkBase

std::ostream& operator<<(std::ostream& os, const PortalLinkBase& link) {
  link.toStream(os);
  return os;
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
    const GeometryContext& gctx, const PortalLinkBase& other,
    BinningValue direction, const Logger& logger) const {
  ACTS_DEBUG("Merging two arbitrary portals");

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

    return mergeImpl(gctx, other, surfaceA, surfaceB, direction, logger);

  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&surfaceA);
             discA != nullptr) {
    const auto* discB = dynamic_cast<const DiscSurface*>(&surfaceB);
    throw_assert(discB != nullptr,
                 "Cannot merge DiscSurface with non-DiscSurface");
    throw_assert(
        direction == binR || direction == binPhi,
        "Invalid binning direction: " + binningValueNames()[direction]);

    return mergeImpl(gctx, other, surfaceA, surfaceB, direction, logger);

  } else {
    throw std::logic_error{"Surface type is not supported"};
  }
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::mergeImpl(
    const GeometryContext& gctx, const PortalLinkBase& other,
    const RegularSurface& surfaceA, const RegularSurface& surfaceB,
    BinningValue direction, const Logger& logger) const {
  ACTS_VERBOSE("Binary portal merging");
  throw std::logic_error{"Not implemented"};
}

// MARK: - GridPortalLinks

// @TODO: Revisit if WARNING messages are the way to go here

namespace {

std::unique_ptr<PortalLinkBase> mergeGridPortals(
    const GeometryContext& gctx, const GridPortalLink* a,
    const GridPortalLink* b, const CylinderSurface* surfaceA,
    const CylinderSurface* surfaceB, BinningValue direction,
    const Logger& logger) {
  assert(surfaceA != nullptr);
  assert(surfaceB != nullptr);
  assert(a->dim() == 2 || a->dim() == 1);
  assert(a->dim() == b->dim());

  constexpr auto tolerance = s_onSurfaceTolerance;

  if (a->dim() == 1) {
    ACTS_VERBOSE("Merge two 1D GridPortalLinks on CylinderSurfaces in "
                 << binningValueNames()[direction]);
    if (a->direction() != b->direction()) {
      ACTS_WARNING("GridPortalLinks have different directions");
      return nullptr;
    }

    auto [mergedSurface, reversed] =
        surfaceA->mergedWith(gctx, *surfaceB, direction);

    ACTS_VERBOSE("Merged surface: " << mergedSurface->toStream(gctx));

    // Normalize ordering of grid portals and surfaces: a is always at lower
    // range than b
    if (reversed) {
      std::swap(surfaceA, surfaceB);
      std::swap(a, b);
    }

    if (direction == binZ) {
      ACTS_VERBOSE("Grids are binned along "
                   << binningValueNames()[a->direction()]);
      if (a->direction() == binZ) {
        ACTS_VERBOSE("=> colinear merge");

        const auto& axisA = *a->grid().axes().front();
        const auto& axisB = *b->grid().axes().front();

        if (axisA.getBoundaryType() != axisB.getBoundaryType()) {
          ACTS_WARNING("AxisBoundaryTypes are different");
          return nullptr;
        }

        if (axisA.getBoundaryType() != AxisBoundaryType::Bound) {
          ACTS_WARNING(
              "AxisBoundaryType is not Bound, cannot do colinear merge");
          return nullptr;
        }

        auto mergeVariable = [&, mergedSurface](const auto& axisA,
                                                const auto& axisB) {
          ActsScalar halfWidth = (axisA.getMax() - axisA.getMin() +
                                  axisB.getMax() - axisB.getMin()) /
                                 2.0;

          ActsScalar shift = axisA.getMax() - halfWidth;
          ACTS_VERBOSE("    ~> shift: " << shift);

          std::vector<ActsScalar> binEdges;

          binEdges.reserve(axisA.getNBins() + axisB.getNBins() + 1);
          auto edgesA = axisA.getBinEdges();
          std::transform(edgesA.begin(), edgesA.end(),
                         std::back_inserter(binEdges),
                         [&](ActsScalar edge) { return edge + shift; });

          ActsScalar stitchPoint = binEdges.back();
          auto edgesB = axisB.getBinEdges();
          std::transform(std::next(edgesB.begin()), edgesB.end(),
                         std::back_inserter(binEdges), [&](ActsScalar edge) {
                           return edge - axisB.getMin() + stitchPoint;
                         });

          Axis merged{AxisBound, std::move(binEdges)};
          ACTS_VERBOSE("    ~> merged axis: " << merged);

          return GridPortalLink::make(*mergedSurface, binZ, std::move(merged));
        };

        AxisType aType = axisA.getType();
        AxisType bType = axisB.getType();
        if (aType == AxisType::Equidistant && bType == AxisType::Equidistant) {
          ACTS_VERBOSE(
              "===> potentially equidistant merge: checking bin widths");

          ActsScalar binsWidthA =
              (axisA.getMax() - axisA.getMin()) / axisA.getNBins();
          ActsScalar binsWidthB =
              (axisB.getMax() - axisB.getMin()) / axisB.getNBins();

          ACTS_VERBOSE("  ~> binWidths: " << binsWidthA << " vs "
                                          << binsWidthB);

          if (std::abs(binsWidthA - binsWidthB) < tolerance) {
            ACTS_VERBOSE("==> binWidths same: " << binsWidthA);

            ActsScalar halfWidth = (axisA.getMax() - axisA.getMin() +
                                    axisB.getMax() - axisB.getMin()) /
                                   2.0;
            Axis merged{AxisBound, -halfWidth, halfWidth,
                        axisA.getNBins() + axisB.getNBins()};

            ACTS_VERBOSE("    ~> merged axis: " << merged);

            std::unique_ptr<PortalLinkBase> mergedPortalLink =
                GridPortalLink::make(*mergedSurface, binZ, std::move(merged));

            // @TODO: Sync bin contents
            return mergedPortalLink;

          } else {
            ACTS_VERBOSE("==> binWidths differ: " << binsWidthA << " vs "
                                                  << binsWidthB
                                                  << " ~> variable merge");

            std::unique_ptr<PortalLinkBase> mergedPortalLink =
                mergeVariable(axisA, axisB);
            // @TODO: Sync bin contents
            return mergedPortalLink;
          }

        } else if (aType == AxisType::Variable && bType == AxisType::Variable) {
          ACTS_VERBOSE("===> variable merge");
          std::unique_ptr<PortalLinkBase> mergedPortalLink =
              mergeVariable(axisA, axisB);
          // @TODO: Sync bin contents
          return mergedPortalLink;
        } else if (aType == AxisType::Equidistant &&
                   bType == AxisType::Variable) {
          ACTS_WARNING("===> mixed merged");
          std::unique_ptr<PortalLinkBase> mergedPortalLink =
              mergeVariable(axisA, axisB);
          // @TODO: Sync bin contents
          return mergedPortalLink;
        } else {
          ACTS_WARNING("===> mixed merged");
          std::unique_ptr<PortalLinkBase> mergedPortalLink =
              mergeVariable(axisA, axisB);
          // @TODO: Sync bin contents
          return mergedPortalLink;
        }

        return nullptr;

      } else {
        ACTS_VERBOSE("=> perpendicular merge");
        return nullptr;
      }

      // AxisBoundaryType aBoundaryType =
      //     a.grid().axes().front()->getBoundaryType();
      // AxisBoundaryType bBoundaryType =
      //     b.grid().axes().front()->getBoundaryType();
      // ACTS_VERBOSE("AxisBoundaryTypes are:");
      // ACTS_VERBOSE(" - a: " << aBoundaryType);
      // ACTS_VERBOSE(" - b: " << bBoundaryType);
      //
      // if(aBoundaryType != bBoundaryType) {
      //   ACTS_WARNING("AxisBoundaryTypes are different");
      //   return nullptr;
      // }
      // else if(aBoundaryType )
      //
      // }
    } else if (direction == binRPhi) {
      ACTS_VERBOSE("BINRPHI");
      // Linear merge along rphi will NOT wrap around (doesn't make sense)
      // Cross merge might have wrap around
      throw std::logic_error{"Not implemented"};
    } else {
      ACTS_ERROR(
          "Invalid binning direction: " << binningValueNames()[direction]);
      throw std::invalid_argument{"Invalid binning direction"};
    }
  } else {
    ACTS_WARNING("2D grid merging is not implemented");
    return nullptr;
  }
}

std::unique_ptr<PortalLinkBase> mergeGridPortals(
    const GeometryContext& gctx, const GridPortalLink* a,
    const GridPortalLink* b, const DiscSurface* surfaceA,
    const DiscSurface* surfaceB, BinningValue direction, const Logger& logger) {
  assert(surfaceA != nullptr);
  assert(surfaceB != nullptr);
  ACTS_WARNING("Disc grid portal merging is not implemented");
  return nullptr;
}

std::unique_ptr<PortalLinkBase> mergeGridPortals(const GeometryContext& gctx,
                                                 const GridPortalLink* a,
                                                 const GridPortalLink* b,
                                                 const RegularSurface& surfaceA,
                                                 const RegularSurface& surfaceB,
                                                 BinningValue direction,
                                                 const Logger& logger) {
  assert(a->dim() == 2 || a->dim() == 1);
  assert(b -.dim() == 2 || b -.dim() == 1);

  if (a->dim() < b->dim()) {
    return mergeGridPortals(gctx, b, a, surfaceB, surfaceA, direction, logger);
  }

  ACTS_VERBOSE("Merging GridPortalLinks along "
               << binningValueNames()[direction] << ":");
  ACTS_VERBOSE(" - a: " << a->grid()
                        << " along: " << binningValueNames()[a->direction()]);
  ACTS_VERBOSE(" - b: " << b->grid()
                        << " along: " << binningValueNames()[b->direction()]);

  if (a->dim() == b->dim()) {
    ACTS_VERBOSE("Grid both have same dimension: " << a->dim());

    if (const auto* cylinder = dynamic_cast<const CylinderSurface*>(&surfaceA);
        cylinder != nullptr) {
      return mergeGridPortals(gctx, a, b, cylinder,
                              &dynamic_cast<const CylinderSurface&>(surfaceB),
                              direction, logger);
    } else if (const auto* disc = dynamic_cast<const DiscSurface*>(&surfaceA);
               disc != nullptr) {
      return mergeGridPortals(gctx, a, b, disc,
                              &dynamic_cast<const DiscSurface&>(surfaceB),
                              direction, logger);
    } else {
      ACTS_VERBOSE("Surface type is not supported here, falling back");
      return nullptr;
    }
  } else {
    ACTS_VERBOSE("Grids have different dimension, falling back");
    return nullptr;
  }
}
}  // namespace

std::unique_ptr<PortalLinkBase> GridPortalLink::mergeImpl(
    const GeometryContext& gctx, const PortalLinkBase& other,
    const RegularSurface& surfaceA, const RegularSurface& surfaceB,
    BinningValue direction, const Logger& logger) const {
  ACTS_VERBOSE("this: GridPortalLink<" << dim() << ">");
  ACTS_VERBOSE("Is other also GridPortalLink?");
  if (const auto* gridPortalLink =
          dynamic_cast<const GridPortalLink*>(&other)) {
    ACTS_VERBOSE("-> yes!");
    auto merged = mergeGridPortals(gctx, this, gridPortalLink, surfaceA,
                                   surfaceB, direction, logger);

    if (merged != nullptr) {
      return merged;
    }
    ACTS_VERBOSE("Grid merging failed, falling back to binary merging");
  } else {
    ACTS_VERBOSE("-> no! Falling back to binary merging");
  }
  return PortalLinkBase::mergeImpl(gctx, other, surfaceA, surfaceB, direction,
                                   logger);
}

}  // namespace Acts
