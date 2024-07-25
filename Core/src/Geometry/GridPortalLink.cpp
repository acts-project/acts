// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GridPortalLink.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <stdexcept>
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
std::unique_ptr<GridPortalLink> makeGrid(
    const std::shared_ptr<surface_t>& surface, BinningValue direction,
    const Logger& logger, std::tuple<Args...> args,
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
      if (surface->bounds().coversFullAzimuth()) {
        return makeGrid(AxisClosed);
      }
    } else if (std::is_same_v<DiscSurface, surface_t>) {
      if (dynamic_cast<const RadialBounds&>(surface->bounds())
              .coversFullAzimuth()) {
        return makeGrid(AxisClosed);
      }
    }
  }

  return makeGrid(AxisBound);
}

template <PortalSurfaceConcept surface_t, typename other_axis_t>
std::unique_ptr<GridPortalLink> mergeVariable(
    const std::shared_ptr<surface_t>& mergedSurface, const IAxis& axisA,
    const IAxis& axisB, ActsScalar /*tolerance*/, BinningValue direction,
    const Logger& logger, const other_axis_t& otherAxis) {
  ACTS_VERBOSE("Variable merge: direction is " << direction);

  std::vector<ActsScalar> binEdges;

  binEdges.reserve(axisA.getNBins() + axisB.getNBins() + 1);

  auto edgesA = axisA.getBinEdges();

  if (direction == BinningValue::binR) {
    ACTS_VERBOSE("Performing asymmetric merge");
    std::copy(edgesA.begin(), edgesA.end(), std::back_inserter(binEdges));

  } else {
    ACTS_VERBOSE("Performing symmetrized merge");
    ActsScalar halfWidth =
        (axisA.getMax() - axisA.getMin() + axisB.getMax() - axisB.getMin()) /
        2.0;
    ACTS_VERBOSE("    ~> half width: " << halfWidth);

    ActsScalar shift = axisA.getMax() - halfWidth;
    ACTS_VERBOSE("    ~> shift: " << shift);

    std::transform(edgesA.begin(), edgesA.end(), std::back_inserter(binEdges),
                   [&](ActsScalar edge) { return edge + shift; });
  }

  ActsScalar stitchPoint = binEdges.back();
  auto edgesB = axisB.getBinEdges();
  std::transform(
      std::next(edgesB.begin()), edgesB.end(), std::back_inserter(binEdges),
      [&](ActsScalar edge) { return edge - axisB.getMin() + stitchPoint; });

  return makeGrid(mergedSurface, direction, logger,
                  std::tuple{std::move(binEdges)}, otherAxis);
}

template <PortalSurfaceConcept surface_t, typename other_axis_t>
std::unique_ptr<GridPortalLink> mergeEquidistant(
    const std::shared_ptr<surface_t>& mergedSurface, const IAxis& axisA,
    const IAxis& axisB, ActsScalar tolerance, BinningValue direction,
    const Logger& logger, other_axis_t otherAxis) {
  ACTS_VERBOSE("===> potentially equidistant merge: checking bin widths");

  ActsScalar binsWidthA = (axisA.getMax() - axisA.getMin()) / axisA.getNBins();
  ActsScalar binsWidthB = (axisB.getMax() - axisB.getMin()) / axisB.getNBins();

  ACTS_VERBOSE("  ~> binWidths: " << binsWidthA << " vs " << binsWidthB);

  if (std::abs(binsWidthA - binsWidthB) < tolerance) {
    ACTS_VERBOSE("==> binWidths same: " << binsWidthA);

    ActsScalar min = std::numeric_limits<ActsScalar>::signaling_NaN();
    ActsScalar max = std::numeric_limits<ActsScalar>::signaling_NaN();

    if (direction == BinningValue::binR) {
      ACTS_VERBOSE("Performing asymmetric merge");
      min = axisA.getMin();
      max = axisB.getMax();
    } else {
      ACTS_VERBOSE("Performing symmetrized merge");

      ActsScalar halfWidth =
          (axisA.getMax() - axisA.getMin() + axisB.getMax() - axisB.getMin()) /
          2.0;

      min = -halfWidth;
      max = halfWidth;
    }

    return makeGrid(mergedSurface, direction, logger,
                    std::tuple{min, max, axisA.getNBins() + axisB.getNBins()},
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

template <PortalSurfaceConcept surface_t, typename other_axis_t>
std::unique_ptr<GridPortalLink> colinearMerge(
    const std::shared_ptr<surface_t>& mergedSurface, const IAxis& axisA,
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

template <PortalSurfaceConcept surface_t>
struct SurfaceInfo;

template <>
struct SurfaceInfo<CylinderSurface> {
  static constexpr std::string_view name = "CylinderSurface";
  static constexpr BinningValue loc0 = BinningValue::binRPhi;
  static constexpr BinningValue loc1 = BinningValue::binZ;
};

template <>
struct SurfaceInfo<DiscSurface> {
  static constexpr std::string_view name = "DiscSurface";
  static constexpr BinningValue loc0 = BinningValue::binR;
  static constexpr BinningValue loc1 = BinningValue::binPhi;
};

template <PortalSurfaceConcept surface_t>
std::unique_ptr<PortalLinkBase> mergeGridPortals(
    const GridPortalLink* a, const GridPortalLink* b, const surface_t* surfaceA,
    const surface_t* surfaceB, BinningValue direction, const Logger& logger) {
  assert(surfaceA != nullptr);
  assert(surfaceB != nullptr);
  assert(a->dim() == 2 || a->dim() == 1);
  assert(a->dim() == b->dim());

  constexpr std::string_view name = SurfaceInfo<surface_t>::name;
  constexpr BinningValue loc0 = SurfaceInfo<surface_t>::loc0;
  constexpr BinningValue loc1 = SurfaceInfo<surface_t>::loc1;

  constexpr auto tolerance = s_onSurfaceTolerance;

  auto [mergedSurface, reversed] =
      surfaceA->mergedWith(*surfaceB, direction, true, logger);
  ACTS_VERBOSE("Merged surface: " << *mergedSurface);

  // Normalize ordering of grid portals and surfaces: a is always at lower
  // range than b
  if (reversed) {
    ACTS_VERBOSE("Swapping grid portals and surfaces after merging");
    std::swap(surfaceA, surfaceB);
    std::swap(a, b);
  }

  // We do this here, because this is the only place where we've normalized the
  // ordering of grids just above: a is always "before" b, so we can iterate
  // over the bins in a unified way to fill the merged grid.
  auto fillGrid = [&](auto merged) {
    if (auto* mergedGrid = dynamic_cast<GridPortalLink*>(merged.get());
        mergedGrid != nullptr) {
      ACTS_VERBOSE("Post processing merged grid: " << mergedGrid->grid());
      GridPortalLink::fillMergedGrid(*a, *b, *mergedGrid, direction, logger);
    }
    return merged;
  };

  if (a->dim() == 1) {
    ACTS_VERBOSE("Merge two 1D GridPortalLinks on " << name << "s in "
                                                    << direction);
    if (a->direction() != b->direction()) {
      ACTS_VERBOSE("GridPortalLinks have different directions");
      ACTS_VERBOSE("=> cross merge");

      // Find the one which is already binned in the merging direction
      const auto* aligned = a->direction() == direction ? a : b;
      const auto* other = a->direction() == direction ? b : a;

      // Extend the aligned one by the other one's axis
      auto aligned2D = aligned->make2DGrid(other->grid().axes().front());
      // Edtend the other one with a single bin
      auto other2D = other->make2DGrid(nullptr);

      assert(aligned2D != nullptr);
      assert(other2D != nullptr);
      return mergeGridPortals(aligned2D.get(), other2D.get(), surfaceA,
                              surfaceB, direction, logger);
    }

    const auto& axisA = *a->grid().axes().front();
    const auto& axisB = *b->grid().axes().front();

    if (direction == loc0) {
      ACTS_VERBOSE("Grids are binned along " << a->direction());
      if (a->direction() == loc0) {
        ACTS_VERBOSE("=> colinear merge");

        return fillGrid(colinearMerge(mergedSurface, axisA, axisB, tolerance,
                                      loc0, logger, NoOtherAxis{}));

      } else {
        ACTS_VERBOSE("=> parallel merge");

        auto a2D = a->make2DGrid(nullptr);
        auto b2D = b->make2DGrid(nullptr);
        assert(a2D != nullptr);
        assert(b2D != nullptr);
        return mergeGridPortals(a2D.get(), b2D.get(), surfaceA, surfaceB,
                                direction, logger);
      }
    } else if (direction == loc1) {
      ACTS_VERBOSE("Grids are binned along " << a->direction());
      if (a->direction() == loc1) {
        ACTS_VERBOSE("=> colinear merge");

        return fillGrid(colinearMerge(mergedSurface, axisA, axisB, tolerance,
                                      loc1, logger, NoOtherAxis{}));

      } else {
        ACTS_VERBOSE("=> parallel merge");
        auto a2D = a->make2DGrid(nullptr);
        auto b2D = b->make2DGrid(nullptr);
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
    ACTS_VERBOSE("Merging two 2D GridPortalLinks on "
                 << name << "s in " << direction << " direction");

    const auto& loc0AxisA = *a->grid().axes().front();
    const auto& loc1AxisA = *a->grid().axes().back();
    const auto& loc0AxisB = *b->grid().axes().front();
    const auto& loc1AxisB = *b->grid().axes().back();

    if (direction == loc0) {
      ACTS_VERBOSE("=> colinear merge along " << loc0);
      ACTS_VERBOSE("--> Checking if " << loc1 << " axes are identical");

      if (loc1AxisA != loc1AxisB) {
        ACTS_WARNING("    ~> "
                     << loc1
                     << " axes are not identical, falling back to composite "
                        "merging");
        return nullptr;
      }
      ACTS_VERBOSE("    ~> they are!");

      return fillGrid(loc1AxisA.visit(
          [&, mergedSurface](auto axis) -> std::unique_ptr<GridPortalLink> {
            ACTS_VERBOSE("    ~> " << loc1 << " axis: " << axis);
            return colinearMerge(mergedSurface, loc0AxisA, loc0AxisB, tolerance,
                                 loc0, logger, AppendAxis{axis});
          }));

    } else if (direction == loc1) {
      ACTS_VERBOSE("=> colinear merge along " << loc1);
      ACTS_VERBOSE("--> Checking if " << loc0 << " axes are identical");

      if (loc0AxisA != loc0AxisB) {
        ACTS_WARNING("    ~> "
                     << loc0
                     << " axes are not identical, falling back to composite "
                        "merging");
        return nullptr;
      }
      ACTS_VERBOSE("    ~> they are!");

      return fillGrid(loc0AxisA.visit(
          [&, mergedSurface](auto axis) -> std::unique_ptr<GridPortalLink> {
            ACTS_VERBOSE("    ~> rPhi axis: " << axis);
            return colinearMerge(mergedSurface, loc1AxisA, loc1AxisB, tolerance,
                                 loc1, logger, PrependAxis{axis});
          }));

    } else {
      ACTS_ERROR("Invalid binning direction: " << a->direction());
      throw std::invalid_argument{"Invalid binning direction"};
    }
  }
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

  const auto* cylinder = dynamic_cast<const CylinderSurface*>(&surfaceA);
  const auto* disc = dynamic_cast<const DiscSurface*>(&surfaceA);

  if (a->dim() == b->dim()) {
    ACTS_VERBOSE("Grid both have same dimension: " << a->dim());

    if (cylinder != nullptr) {
      return mergeGridPortals(a, b, cylinder,
                              &dynamic_cast<const CylinderSurface&>(surfaceB),
                              direction, logger);
    } else if (disc != nullptr) {
      return mergeGridPortals(a, b, disc,
                              &dynamic_cast<const DiscSurface&>(surfaceB),
                              direction, logger);
    } else {
      // @TODO: Support PlaneSurface
      ACTS_VERBOSE("Surface type is not supported here, falling back");
      return nullptr;
    }
  } else {
    ACTS_VERBOSE("Grids have different dimension, extending rhs to 2D");
    const IAxis* otherAxis = nullptr;

    if (b->direction() == direction) {
      ACTS_VERBOSE("1D grid is binned in merging direction " << direction);
      ACTS_VERBOSE("~> Adding complementary axis");

      if (cylinder != nullptr) {
        otherAxis = direction == BinningValue::binRPhi
                        ? a->grid().axes().back()
                        : a->grid().axes().front();
      } else if (disc != nullptr) {
        otherAxis = direction == BinningValue::binR ? a->grid().axes().back()
                                                    : a->grid().axes().front();
      } else {
        ACTS_VERBOSE("Surface type is not supported here, falling back");
        return nullptr;
      }
    } else {
      ACTS_VERBOSE("1D grid is binned in complementary direction");
    }

    auto b2D = b->make2DGrid(otherAxis);
    ACTS_VERBOSE("-> new grid: " << b2D->grid());
    return mergeGridPortals(a, b2D.get(), surfaceA, surfaceB, direction,
                            logger);
  }
}

enum class FillDirection {
  loc0,
  loc1,
};

template <FillDirection dir, typename axis_1_t, typename axis_2_t,
          typename axis_3_t>
void fillGrid1dTo2d(const Grid<const TrackingVolume*, axis_1_t>& grid1d,
                    Grid<const TrackingVolume*, axis_2_t, axis_3_t>& grid2d) {
  const auto locSource = grid1d.numLocalBins();
  const auto locDest = grid2d.numLocalBins();

  for (std::size_t i = 0; i <= locSource[0] + 1; ++i) {
    const auto* source = grid1d.atLocalBins({i});

    if constexpr (dir == FillDirection::loc1) {
      for (std::size_t j = 0; j <= locDest[1] + 1; ++j) {
        grid2d.atLocalBins({i, j}) = source;
      }
    } else if constexpr (dir == FillDirection::loc0) {
      for (std::size_t j = 0; j <= locDest[0] + 1; ++j) {
        grid2d.atLocalBins({j, i}) = source;
      }
    }
  }
}

}  // namespace

void GridPortalLink::fillMergedGrid(const GridPortalLink& a,
                                    const GridPortalLink& b,
                                    GridPortalLink& merged,
                                    BinningValue direction,
                                    const Logger& logger) {
  ACTS_VERBOSE("Filling merged grid");
  assert(a.dim() == b.dim());
  assert(a.direction() == b.direction());
  const auto locBinsA = a.numLocalBins();
  const auto locBinsB = b.numLocalBins();

  ACTS_VERBOSE("a: " << a.grid());
  ACTS_VERBOSE("b: " << b.grid());
  ACTS_VERBOSE("merged: " << merged.grid());

  if (a.dim() == 1) {
    std::size_t nBinsA = locBinsA.at(0);
    std::size_t nBinsB = locBinsB.at(0);

    ACTS_VERBOSE("1d merge:");

    for (std::size_t i = 1; i <= nBinsA; ++i) {
      merged.atLocalBins({i}) = a.atLocalBins({i});
      ACTS_VERBOSE(" ~> a " << i << " -> m " << i << " (" << a.atLocalBins({i})
                            << ")");
    }
    for (std::size_t i = 1; i <= nBinsB; ++i) {
      merged.atLocalBins({nBinsA + i}) = b.atLocalBins({i});
      ACTS_VERBOSE(" ~> b " << i << " -> m " << (nBinsA + i) << " ("
                            << b.atLocalBins({i}) << ")");
    }
  } else {
    ACTS_VERBOSE("2d merge:");
    if (a.direction() == direction) {
      ACTS_VERBOSE("~> direction is loc0");
      std::size_t nBinsA = locBinsA.at(0);
      std::size_t nBinsB = locBinsB.at(0);
      std::size_t nBinsCommon = locBinsB.at(1);

      for (std::size_t i = 1; i <= nBinsA; ++i) {
        for (std::size_t j = 1; j <= nBinsCommon; ++j) {
          merged.atLocalBins({i, j}) = a.atLocalBins({i, j});
          ACTS_VERBOSE(" ~> a " << i << ", " << j << " -> m " << i << ", " << j
                                << " (" << a.atLocalBins({i, j}) << ")");
        }
      }

      for (std::size_t i = 1; i <= nBinsB; ++i) {
        for (std::size_t j = 1; j <= nBinsCommon; ++j) {
          std::size_t ti = i + nBinsA;
          merged.atLocalBins({ti, j}) = b.atLocalBins({i, j});
          ACTS_VERBOSE(" ~> b " << i << ", " << j << " -> m " << ti << ", " << j
                                << " (" << b.atLocalBins({i, j}) << ")");
        }
      }
    } else {
      ACTS_VERBOSE("~> direction is loc1");
      std::size_t nBinsA = locBinsA.at(1);
      std::size_t nBinsB = locBinsB.at(1);
      std::size_t nBinsCommon = locBinsB.at(0);

      for (std::size_t i = 1; i <= nBinsCommon; ++i) {
        for (std::size_t j = 1; j <= nBinsA; ++j) {
          merged.atLocalBins({i, j}) = a.atLocalBins({i, j});
          ACTS_VERBOSE(" ~> a " << i << ", " << j << " -> m " << i << ", " << j
                                << " (" << a.atLocalBins({i, j}) << ")");
        }
      }

      for (std::size_t i = 1; i <= nBinsCommon; ++i) {
        for (std::size_t j = 1; j <= nBinsB; ++j) {
          std::size_t tj = j + nBinsA;
          merged.atLocalBins({i, tj}) = b.atLocalBins({i, j});
          ACTS_VERBOSE(" ~> b " << i << ", " << j << " -> m " << tj << ", " << j
                                << " (" << b.atLocalBins({i, j}) << ")");
        }
      }
    }
  }

  //   auto swizzle = [&](std::size_t i, std::size_t j) {
  //     return a.direction() == direction ? std::pair{i, j} : std::pair{j,
  //     i};
  //   };
  //
  //   auto [ai, aj] = swizzle(0, 1);
  //
  //   std::size_t nBinsA = locBinsA.at(ai);
  //   std::size_t nBinsB = locBinsB.at(ai);
  //   std::size_t nBinsC = locBinsB.at(aj);
  //   assert(locBinsA.at(aj) == locBinsB.at(aj));
  //   std::cout << "nBinsA: " << nBinsA << ", nBinsB: " << nBinsB
  //             << ", nBinsC: " << nBinsC << std::endl;
  //   std::cout << locBinsA.at(aj) << std::endl;
  //
  //   for (std::size_t i = 1; i <= nBinsA; ++i) {
  //     for (std::size_t j = 1; j <= nBinsC; ++j) {
  //       auto [li, lj] = swizzle(i, j);
  //       merged.atLocalBins({li, lj}) = a.atLocalBins({li, lj});
  //       ACTS_VERBOSE(" ~> a " << li << ", " << lj << " -> m " << li << ", "
  //                             << lj << " (" << a.atLocalBins({li, lj}) <<
  //                             ")");
  //     }
  //   }
  //   for (std::size_t i = 1; i <= nBinsB; ++i) {
  //     for (std::size_t j = 1; j <= nBinsC; ++j) {
  //       auto [li, lj] = swizzle(i, j);
  //
  //       std::size_t ti = li;
  //       std::size_t tj = lj;
  //       // if (b.direction() == direction) {
  //       // ti += nBinsA;
  //       // } else {
  //       //   tj += nBinsA;
  //       // }
  //
  //       std::tie(ti, tj) = swizzle(ti, tj);
  //       merged.atLocalBins({ti, tj}) = b.atLocalBins({li, lj});
  //       ACTS_VERBOSE(" ~> b " << li << ", " << lj << " -> m " << ti << ", "
  //                             << tj << " (" << b.atLocalBins({li, lj}) <<
  //                             ")");
  //     }
  //   }
  // }
}

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
    throw std::logic_error{"Not implemented"};
  } else {
    throw std::invalid_argument{"Surface type is not supported"};
  }

  assert(grid != nullptr);
  grid->setVolume(volume);

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

std::unique_ptr<GridPortalLink> GridPortalLink::extendTo2D(
    const std::shared_ptr<CylinderSurface>& surface, const IAxis* other) const {
  assert(dim() == 1);
  if (direction() == BinningValue::binRPhi) {
    const auto& axisRPhi = *grid().axes().front();
    // 1D direction is binRPhi, so add a Z axis
    ActsScalar hlZ = surface->bounds().get(CylinderBounds::eHalfLengthZ);

    Axis axisZ{AxisBound, -hlZ, hlZ, 1};

    if (other == nullptr) {
      other = &axisZ;
    }

    return axisRPhi.visit([&](const auto& axis0) {
      const auto& self =
          dynamic_cast<const GridPortalLinkT<std::decay_t<decltype(axis0)>>&>(
              *this);
      return other->visit(
          [&](const auto& axis1) -> std::unique_ptr<GridPortalLink> {
            auto grid = GridPortalLink::make(surface, axis0, axis1);
            fillGrid1dTo2d<FillDirection::loc1>(self.grid(), grid->grid());
            return grid;
          });
    });

  } else {
    const auto& axisZ = *grid().axes().front();
    // 1D direction is binZ, so add an rPhi axis
    ActsScalar r = surface->bounds().get(CylinderBounds::eR);
    ActsScalar hlPhi = surface->bounds().get(CylinderBounds::eHalfPhiSector);
    ActsScalar hlRPhi = r * hlPhi;

    auto axis = [&](auto bdt) {
      Axis axisRPhi{bdt, -hlRPhi, hlRPhi, 1};

      if (other == nullptr) {
        other = &axisRPhi;
      }

      return axisZ.visit([&](const auto& axis1) {
        const auto& self =
            dynamic_cast<const GridPortalLinkT<std::decay_t<decltype(axis1)>>&>(
                *this);
        return other->visit(
            [&](const auto& axis0) -> std::unique_ptr<GridPortalLink> {
              auto grid = GridPortalLink::make(surface, axis0, axis1);
              fillGrid1dTo2d<FillDirection::loc0>(self.grid(), grid->grid());
              return grid;
            });
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
    const std::shared_ptr<DiscSurface>& surface, const IAxis* other) const {
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
      Axis axisPhi{bdt, -hlPhi, hlPhi, 1};

      if (other == nullptr) {
        other = &axisPhi;
      }

      return axisR.visit([&](const auto& axis0) {
        const auto& self =
            dynamic_cast<const GridPortalLinkT<std::decay_t<decltype(axis0)>>&>(
                *this);
        return other->visit(
            [&](const auto& axis1) -> std::unique_ptr<GridPortalLink> {
              auto grid = GridPortalLink::make(surface, axis0, axis1);
              fillGrid1dTo2d<FillDirection::loc1>(self.grid(), grid->grid());
              return grid;
            });
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

    Axis axisR{AxisBound, rMin, rMax, 1};

    if (other == nullptr) {
      other = &axisR;
    }

    return axisPhi.visit([&](const auto& axis1) {
      const auto& self =
          dynamic_cast<const GridPortalLinkT<std::decay_t<decltype(axis1)>>&>(
              *this);
      return other->visit(
          [&](const auto& axis0) -> std::unique_ptr<GridPortalLink> {
            auto grid = GridPortalLink::make(surface, axis0, axis1);
            fillGrid1dTo2d<FillDirection::loc0>(self.grid(), grid->grid());
            return grid;
          });
    });
  }
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
    ACTS_VERBOSE("Grid merging failed, falling back to composite merging");
  } else {
    ACTS_VERBOSE("-> no! Falling back to composite merging");
  }
  return PortalLinkBase::mergeImpl(other, surfaceA, surfaceB, direction,
                                   logger);
}

}  // namespace Acts
