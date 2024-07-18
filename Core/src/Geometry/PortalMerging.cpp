// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

#include <type_traits>

namespace Acts {

namespace {

struct NoOtherAxis : public std::false_type {};

template <typename axis_t>
struct OtherAxis : public std::true_type {
  using AxisType = std::decay_t<axis_t>;
  const AxisType& m_axis;

  OtherAxis(const AxisType& axis) : m_axis(axis) {}
};

template <typename axis_t>
OtherAxis(axis_t&&) -> OtherAxis<std::decay_t<axis_t>>;

template <typename axis_t>
struct PrependAxis : public OtherAxis<std::decay_t<axis_t>> {
  static constexpr bool prepend = true;
};

template <typename axis_t>
PrependAxis(axis_t&&) -> PrependAxis<std::decay_t<axis_t>>;

template <typename axis_t>
struct AppendAxis : public OtherAxis<std::decay_t<axis_t>> {
  static constexpr bool prepend = false;
};

template <typename axis_t>
AppendAxis(axis_t&&) -> AppendAxis<std::decay_t<axis_t>>;

template <PortalSurfaceConcept surface_t, typename... Args,
          typename other_axis_t>
std::unique_ptr<GridPortalLink> makeGrid(const surface_t& surface,
                                         BinningValue direction,
                                         const Logger& logger,
                                         std::tuple<Args...> args,
                                         const other_axis_t& otherAxis) {
  // @TODO: PlaneSurface support

  ACTS_VERBOSE("Make resulting merged grid");

  // This is to make it possible to construct Axis with the tuple arguments
  auto axisFactory = [](auto&&... axisArgs) {
    return Axis{std::forward<decltype(axisArgs)>(axisArgs)...};
  };

  // Avoid copy-pasting identical code twice below
  auto makeGrid = [&](auto boundaryType) {
    auto axisArgs = std::tuple_cat(std::tuple{boundaryType}, std::move(args));
    auto merged = std::apply(axisFactory, std::move(axisArgs));

    if constexpr (!std::decay_t<other_axis_t>::value) {
      // No other axis
      ACTS_VERBOSE("    ~> merged axis: " << merged);
      return GridPortalLink::make(surface, direction, std::move(merged));
    } else if constexpr (other_axis_t::prepend) {
      // Prepend other axis
      ACTS_VERBOSE("    ~> other axis (prepend): " << otherAxis.m_axis);
      ACTS_VERBOSE("    ~> merged axis: " << merged);
      return GridPortalLink::make(surface, Axis{otherAxis.m_axis},
                                  std::move(merged));
    } else {
      // Append other axis
      ACTS_VERBOSE("    ~> merged axis: " << merged);
      ACTS_VERBOSE("    ~> other axis (append): " << otherAxis.m_axis);
      return GridPortalLink::make(surface, std::move(merged),
                                  Axis{otherAxis.m_axis});
    }
  };

  // Check if we're in the cylinder or disc case, and the resulting bounds wrap
  // around and should have closed binning
  if (direction == BinningValue::binPhi || direction == BinningValue::binRPhi) {
    if constexpr (std::is_same_v<CylinderSurface, surface_t>) {
      if (surface.bounds().coversFullAzimuth()) {
        return makeGrid(AxisClosed);
      }
    } else if (std::is_same_v<DiscSurface, surface_t>) {
      if (dynamic_cast<const RadialBounds&>(surface.bounds())
              .coversFullAzimuth()) {
        return makeGrid(AxisClosed);
      }
    }
  }

  return makeGrid(AxisBound);
}

template <typename other_axis_t>
std::unique_ptr<GridPortalLink> mergeVariable(
    const PortalSurfaceConcept auto& mergedSurface, const IAxis& axisA,
    const IAxis& axisB, ActsScalar /*tolerance*/, BinningValue direction,
    const Logger& logger, const other_axis_t& otherAxis) {
  ActsScalar halfWidth =
      (axisA.getMax() - axisA.getMin() + axisB.getMax() - axisB.getMin()) / 2.0;

  ActsScalar shift = axisA.getMax() - halfWidth;
  ACTS_VERBOSE("    ~> shift: " << shift);

  std::vector<ActsScalar> binEdges;

  binEdges.reserve(axisA.getNBins() + axisB.getNBins() + 1);
  auto edgesA = axisA.getBinEdges();
  std::transform(edgesA.begin(), edgesA.end(), std::back_inserter(binEdges),
                 [&](ActsScalar edge) { return edge + shift; });

  ActsScalar stitchPoint = binEdges.back();
  auto edgesB = axisB.getBinEdges();
  std::transform(
      std::next(edgesB.begin()), edgesB.end(), std::back_inserter(binEdges),
      [&](ActsScalar edge) { return edge - axisB.getMin() + stitchPoint; });

  return makeGrid(mergedSurface, direction, logger,
                  std::tuple{std::move(binEdges)}, otherAxis);
}

template <typename other_axis_t>
std::unique_ptr<GridPortalLink> mergeEquidistant(
    const PortalSurfaceConcept auto& mergedSurface, const IAxis& axisA,
    const IAxis& axisB, ActsScalar tolerance, BinningValue direction,
    const Logger& logger, other_axis_t otherAxis) {
  ACTS_VERBOSE("===> potentially equidistant merge: checking bin widths");

  ActsScalar binsWidthA = (axisA.getMax() - axisA.getMin()) / axisA.getNBins();
  ActsScalar binsWidthB = (axisB.getMax() - axisB.getMin()) / axisB.getNBins();

  ACTS_VERBOSE("  ~> binWidths: " << binsWidthA << " vs " << binsWidthB);

  if (std::abs(binsWidthA - binsWidthB) < tolerance) {
    ACTS_VERBOSE("==> binWidths same: " << binsWidthA);

    ActsScalar halfWidth =
        (axisA.getMax() - axisA.getMin() + axisB.getMax() - axisB.getMin()) /
        2.0;

    return makeGrid(
        mergedSurface, direction, logger,
        std::tuple{-halfWidth, halfWidth, axisA.getNBins() + axisB.getNBins()},
        otherAxis);

  } else {
    ACTS_VERBOSE("==> binWidths differ: " << binsWidthA << " vs " << binsWidthB
                                          << " ~> variable merge");

    std::unique_ptr<GridPortalLink> mergedPortalLink =
        mergeVariable(mergedSurface, axisA, axisB, tolerance, direction, logger,
                      other_axis_t(otherAxis));
    return mergedPortalLink;
  }
}

template <typename other_axis_t>
std::unique_ptr<GridPortalLink> colinearMerge(
    const PortalSurfaceConcept auto& mergedSurface, const IAxis& axisA,
    const IAxis& axisB, ActsScalar tolerance, BinningValue direction,
    const Logger& logger, const other_axis_t& otherAxis) {
  AxisType aType = axisA.getType();
  AxisType bType = axisB.getType();
  if (axisA.getBoundaryType() != axisB.getBoundaryType()) {
    ACTS_WARNING("AxisBoundaryTypes are different");
    return nullptr;
  }

  if (axisA.getBoundaryType() != AxisBoundaryType::Bound) {
    ACTS_WARNING("AxisBoundaryType is not Bound, cannot do colinear merge");
    return nullptr;
  }

  // convenience only
  auto mergeVariableLocal = [&] {
    return mergeVariable(mergedSurface, axisA, axisB, tolerance, direction,
                         logger, otherAxis);
  };

  if (aType == AxisType::Equidistant && bType == AxisType::Equidistant) {
    auto mergedPortalLink =
        mergeEquidistant(mergedSurface, axisA, axisB, tolerance, direction,
                         logger, other_axis_t(otherAxis));
    // @TODO: Sync bin contents
    return mergedPortalLink;
  } else if (aType == AxisType::Variable && bType == AxisType::Variable) {
    ACTS_VERBOSE("===> variable merge");
    auto mergedPortalLink = mergeVariableLocal();
    // @TODO: Sync bin contents
    return mergedPortalLink;
  } else if (aType == AxisType::Equidistant && bType == AxisType::Variable) {
    ACTS_WARNING("===> mixed merged");
    auto mergedPortalLink = mergeVariableLocal();
    // @TODO: Sync bin contents
    return mergedPortalLink;
  } else {
    ACTS_WARNING("===> mixed merged");
    auto mergedPortalLink = mergeVariableLocal();
    // @TODO: Sync bin contents
    return mergedPortalLink;
  }
}

std::unique_ptr<PortalLinkBase> mergeGridPortals(
    const GridPortalLink* a, const GridPortalLink* b,
    const CylinderSurface* surfaceA, const CylinderSurface* surfaceB,
    BinningValue direction, const Logger& logger) {
  assert(surfaceA != nullptr);
  assert(surfaceB != nullptr);
  assert(a->dim() == 2 || a->dim() == 1);
  assert(a->dim() == b->dim());

  constexpr auto tolerance = s_onSurfaceTolerance;

  auto [mergedSurface, reversed] =
      surfaceA->mergedWith(*surfaceB, direction, true, logger);
  ACTS_VERBOSE("Merged surface: " << *mergedSurface);

  // Normalize ordering of grid portals and surfaces: a is always at lower
  // range than b
  if (reversed) {
    std::swap(surfaceA, surfaceB);
    std::swap(a, b);
  }

  if (a->dim() == 1) {
    ACTS_VERBOSE("Merge two 1D GridPortalLinks on CylinderSurfaces in "
                 << direction);
    if (a->direction() != b->direction()) {
      ACTS_WARNING("GridPortalLinks have different directions");
      return nullptr;
    }

    const auto& axisA = *a->grid().axes().front();
    const auto& axisB = *b->grid().axes().front();

    if (direction == BinningValue::binZ) {
      ACTS_VERBOSE("Grids are binned along " << a->direction());
      if (a->direction() == BinningValue::binZ) {
        ACTS_VERBOSE("=> colinear merge");

        return colinearMerge(*mergedSurface, axisA, axisB, tolerance,
                             BinningValue::binZ, logger, NoOtherAxis{});

      } else {
        ACTS_VERBOSE("=> perpendicular merge");
        auto a2D = a->make2DGrid();
        auto b2D = b->make2DGrid();
        assert(a2D != nullptr);
        assert(b2D != nullptr);
        return mergeGridPortals(a2D.get(), b2D.get(), surfaceA, surfaceB,
                                direction, logger);
      }

    } else if (direction == BinningValue::binRPhi) {
      ACTS_VERBOSE("Grids are binned along " << a->direction());
      if (a->direction() == BinningValue::binRPhi) {
        ACTS_VERBOSE("=> colinear merge");

        return colinearMerge(*mergedSurface, axisA, axisB, tolerance,
                             BinningValue::binRPhi, logger, NoOtherAxis{});

      } else {
        ACTS_VERBOSE("=> perpendicular merge");

        auto a2D = a->make2DGrid();
        auto b2D = b->make2DGrid();
        assert(a2D != nullptr);
        assert(b2D != nullptr);
        return mergeGridPortals(a2D.get(), b2D.get(), surfaceA, surfaceB,
                                direction, logger);
      }
    } else {
      ACTS_ERROR("Invalid binning direction: " << direction);
      throw std::invalid_argument{"Invalid binning direction"};
    }
  } else {
    ACTS_VERBOSE("Merging two 2D GridPortalLinks on CylinderSurfaces in "
                 << direction << " direction");

    const auto& rPhiAxisA = *a->grid().axes().front();
    const auto& zAxisA = *a->grid().axes().back();
    const auto& rPhiAxisB = *b->grid().axes().front();
    const auto& zAxisB = *b->grid().axes().back();

    if (direction == BinningValue::binZ) {
      ACTS_VERBOSE("=> colinear merge along z");
      ACTS_VERBOSE("--> Checking if rPhi axes are identical");

      if (rPhiAxisA != rPhiAxisB) {
        ACTS_WARNING(
            "    ~> rPhi axes are not identical, falling back to binary "
            "merging");
        return nullptr;
      }
      ACTS_VERBOSE("    ~> they are!");

      return rPhiAxisA.visit(
          [&, mergedSurface](auto axis) -> std::unique_ptr<GridPortalLink> {
            ACTS_VERBOSE("    ~> rPhi axis: " << axis);
            return colinearMerge(*mergedSurface, zAxisA, zAxisB, tolerance,
                                 BinningValue::binZ, logger, PrependAxis{axis});
          });
    } else if (direction == BinningValue::binRPhi) {
      ACTS_VERBOSE("=> colinear merge along rPhi");
      ACTS_VERBOSE("--> Checking if z axes are identical");

      if (zAxisA != zAxisB) {
        ACTS_WARNING(
            "    ~> z axes are not identical, falling back to binary merging");
        return nullptr;
      }
      ACTS_VERBOSE("    ~> they are!");

      return zAxisA.visit([&, mergedSurface](
                              auto axis) -> std::unique_ptr<GridPortalLink> {
        ACTS_VERBOSE("    ~> z axis: " << axis);
        return colinearMerge(*mergedSurface, rPhiAxisA, rPhiAxisB, tolerance,
                             BinningValue::binRPhi, logger, AppendAxis{axis});
      });
    } else {
      ACTS_ERROR("Invalid binning direction: " << a->direction());
      throw std::invalid_argument{"Invalid binning direction"};
    }
  }
}

std::unique_ptr<PortalLinkBase> mergeGridPortals(const GridPortalLink* a,
                                                 const GridPortalLink* b,
                                                 const DiscSurface* surfaceA,
                                                 const DiscSurface* surfaceB,
                                                 BinningValue direction,
                                                 const Logger& logger) {
  (void)a;
  (void)b;
  (void)surfaceA;
  (void)surfaceB;
  (void)direction;
  assert(surfaceA != nullptr);
  assert(surfaceB != nullptr);
  ACTS_WARNING("Disc grid portal merging is not implemented");
  return nullptr;
}

std::unique_ptr<PortalLinkBase> mergeGridPortals(const GridPortalLink* a,
                                                 const GridPortalLink* b,
                                                 const RegularSurface& surfaceA,
                                                 const RegularSurface& surfaceB,
                                                 BinningValue direction,
                                                 const Logger& logger) {
  assert(a->dim() == 2 || a->dim() == 1);
  assert(b->dim() == 2 || b->dim() == 1);

  if (a->dim() < b->dim()) {
    return mergeGridPortals(b, a, surfaceB, surfaceA, direction, logger);
  }

  ACTS_VERBOSE("Merging GridPortalLinks along " << direction << ":");
  ACTS_VERBOSE(" - a: " << a->grid() << " along: " << a->direction());
  ACTS_VERBOSE(" - b: " << b->grid() << " along: " << b->direction());

  if (a->dim() == b->dim()) {
    ACTS_VERBOSE("Grid both have same dimension: " << a->dim());

    if (const auto* cylinder = dynamic_cast<const CylinderSurface*>(&surfaceA);
        cylinder != nullptr) {
      return mergeGridPortals(a, b, cylinder,
                              &dynamic_cast<const CylinderSurface&>(surfaceB),
                              direction, logger);
    } else if (const auto* disc = dynamic_cast<const DiscSurface*>(&surfaceA);
               disc != nullptr) {
      return mergeGridPortals(a, b, disc,
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

std::unique_ptr<PortalLinkBase> PortalLinkBase::merge(
    const GeometryContext& /*gctx*/, const PortalLinkBase& other,
    BinningValue direction, const Logger& logger) const {
  ACTS_DEBUG("Merging two arbitrary portals");

  ACTS_VERBOSE(" - this:  " << *this);
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
        direction == BinningValue::binZ || direction == BinningValue::binRPhi,
        "Invalid binning direction: " + binningValueName(direction));

    return mergeImpl(other, surfaceA, surfaceB, direction, logger);

  } else if (const auto* discA = dynamic_cast<const DiscSurface*>(&surfaceA);
             discA != nullptr) {
    const auto* discB = dynamic_cast<const DiscSurface*>(&surfaceB);
    throw_assert(discB != nullptr,
                 "Cannot merge DiscSurface with non-DiscSurface");
    throw_assert(
        direction == BinningValue::binR || direction == BinningValue::binPhi,
        "Invalid binning direction: " + binningValueName(direction));

    return mergeImpl(other, surfaceA, surfaceB, direction, logger);

  } else {
    throw std::logic_error{"Surface type is not supported"};
  }
}

std::unique_ptr<PortalLinkBase> PortalLinkBase::mergeImpl(
    const PortalLinkBase& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  (void)other;
  (void)surfaceA;
  (void)surfaceB;
  (void)direction;
  ACTS_VERBOSE("Binary portal merging");
  throw std::logic_error{"Not implemented"};
}

std::unique_ptr<PortalLinkBase> GridPortalLink::mergeImpl(
    const PortalLinkBase& other, const RegularSurface& surfaceA,
    const RegularSurface& surfaceB, BinningValue direction,
    const Logger& logger) const {
  ACTS_VERBOSE("this: GridPortalLink<" << dim() << ">");
  ACTS_VERBOSE("Is other also GridPortalLink?");
  if (const auto* gridPortalLink =
          dynamic_cast<const GridPortalLink*>(&other)) {
    ACTS_VERBOSE("-> yes!");
    auto merged = mergeGridPortals(this, gridPortalLink, surfaceA, surfaceB,
                                   direction, logger);

    if (merged != nullptr) {
      return merged;
    }
    ACTS_VERBOSE("Grid merging failed, falling back to binary merging");
  } else {
    ACTS_VERBOSE("-> no! Falling back to binary merging");
  }
  return PortalLinkBase::mergeImpl(other, surfaceA, surfaceB, direction,
                                   logger);
}

}  // namespace Acts
