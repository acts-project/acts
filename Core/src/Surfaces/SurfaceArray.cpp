// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ranges.hpp"

#include <ranges>
#include <type_traits>
#include <utility>

using namespace Acts::Ranges;
namespace Acts {

// implementation for pure virtual destructor of ISurfaceGridLookup
SurfaceArray::ISurfaceGridLookup::~ISurfaceGridLookup() = default;

SurfaceArray::SurfaceArray(std::unique_ptr<ISurfaceGridLookup> gridLookup,
                           std::vector<std::shared_ptr<const Surface>> surfaces,
                           const Transform3& transform)
    : p_gridLookup(std::move(gridLookup)),
      m_surfaces(std::move(surfaces)),
      m_surfacesRawPointers(unpackSmartPointers(m_surfaces)),
      m_transform(transform) {
  if (p_gridLookup != nullptr) {
    if (const auto& grid = p_gridLookup->getGridView()) {
      checkGrid(grid.value());
    }
  }
}

SurfaceArray::SurfaceArray(std::shared_ptr<const Surface> srf)
    : p_gridLookup(std::make_unique<SingleElementLookup>(srf.get())),
      m_surfaces({std::move(srf)}) {
  m_surfacesRawPointers.push_back(m_surfaces.at(0).get());
}

std::ostream& SurfaceArray::toStream(const GeometryContext& /*gctx*/,
                                     std::ostream& sl) const {
  sl << std::fixed << std::setprecision(4);
  sl << "SurfaceArray:" << std::endl;
  sl << " - no surfaces: " << m_surfaces.size() << std::endl;

  auto axes = p_gridLookup->getAxes();

  for (std::size_t j = 0; j < axes.size(); ++j) {
    AxisBoundaryType bdt = axes.at(j)->getBoundaryType();
    sl << " - axis " << (j + 1) << std::endl;
    sl << "   - boundary type: ";
    if (bdt == AxisBoundaryType::Open) {
      sl << "open";
    }
    if (bdt == AxisBoundaryType::Bound) {
      sl << "bound";
    }
    if (bdt == AxisBoundaryType::Closed) {
      sl << "closed";
    }
    sl << std::endl;
    sl << "   - type: "
       << (axes.at(j)->isEquidistant() ? "equidistant" : "variable")
       << std::endl;
    sl << "   - n bins: " << axes.at(j)->getNBins() << std::endl;
    sl << "   - bin edges: [ ";
    auto binEdges = axes.at(j)->getBinEdges();
    for (std::size_t i = 0; i < binEdges.size(); ++i) {
      if (i > 0) {
        sl << ", ";
      }
      auto binEdge = binEdges.at(i);
      // Do not display negative zeroes
      sl << ((std::abs(binEdge) >= 5e-4) ? binEdge : 0.0);
    }
    sl << " ]" << std::endl;
  }
  return sl;
}

void SurfaceArray::checkGrid(AnyGridConstView<SurfaceVector> grid) {
  std::set allSurfaces =
      m_surfaces |
      std::views::transform([](const auto& sp) { return sp.get(); }) |
      to<std::set>;
  std::set<const Surface*> seenSurface;
  auto bins = grid.numLocalBins();
  for (std::size_t i = 0; i <= bins.at(0); ++i) {
    for (std::size_t j = 0; j <= bins.at(1); ++j) {
      const auto& surfaces = grid.atLocalBins({i, j});
      for (const auto& srf : surfaces) {
        seenSurface.insert(srf);
      }
    }
  }

  if (allSurfaces != seenSurface) {
    std::set<const Surface*> diff;
    std::ranges::set_difference(allSurfaces, seenSurface,
                                std::inserter(diff, diff.begin()));

    throw std::logic_error(
        std::format("SurfaceArray grid does not contain all surfaces provided! "
                    "{} surfaces not seen",
                    diff.size()));
  }
}

namespace {

/// @brief Lookup helper which encapsulates a @c Grid
/// @tparam Axes The axes used for the grid
template <class Axis1, class Axis2>
struct SurfaceGridLookup2 : SurfaceArray::ISurfaceGridLookup {
  /// Grid type storing surface vectors with two axes
  using Grid_t = Grid<SurfaceVector, Axis1, Axis2>;

  /// Construct a surface grid lookup
  /// @param representative The surface which is used as representative
  /// @param tolerance The tolerance used for intersection checks
  /// @param axes The axes used for the grid
  /// @param bValues Optional vector of axis directions for binning
  SurfaceGridLookup2(std::shared_ptr<RegularSurface> representative,
                     double tolerance, std::tuple<Axis1, Axis2> axes,
                     std::vector<AxisDirection> bValues = {})
      : m_representative(std::move(representative)),
        m_tolerance(tolerance),
        m_grid(std::move(axes)),
        m_binValues(std::move(bValues)) {
    m_neighborMap.resize(m_grid.size());
  }

  /// @brief Fill provided surfaces into the contained @c Grid.
  ///
  /// This is done by iterating, accessing the referencePosition, lookup
  /// and append.
  /// Also populates the neighbor map by combining the filled bins of
  /// all bins around a given one.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces Input surface pointers
  void fill(const GeometryContext& gctx,
            const SurfaceVector& surfaces) override {
    for (const Surface* surface : surfaces) {
      const std::size_t globalBin = fillSurfaceToBinMapping(gctx, *surface);
      if (globalBin == 0) {
        continue;
      }

      fillBinToSurfaceMapping(gctx, *surface, globalBin);
    }

    populateNeighborCache();
  }

  const SurfaceVector& lookup(const Vector3& position,
                              const Vector3& direction) const override {
    return m_grid.at(findGlobalBin(position, direction,
                                   std::numeric_limits<double>::infinity()));
  }

  /// @brief Performs lookup at global bin and returns bin content as
  /// reference
  /// @param bin Global lookup bin
  /// @return @c SurfaceVector at given bin
  SurfaceVector& lookup(std::size_t bin) override { return m_grid.at(bin); }

  /// @brief Performs lookup at global bin and returns bin content as const
  /// reference
  /// @param bin Global lookup bin
  /// @return @c SurfaceVector at given bin
  const SurfaceVector& lookup(std::size_t bin) const override {
    return m_grid.at(bin);
  }

  /// @brief Performs a lookup at @c pos, but returns neighbors as well
  ///
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return @c SurfaceVector at given bin. Copy of all bins selected
  const SurfaceVector& neighbors(const Vector3& position,
                                 const Vector3& direction) const override {
    return m_neighborMap.at(findGlobalBin(
        position, direction, std::numeric_limits<double>::infinity()));
  }

  /// @brief Returns the total size of the grid (including under/overflow
  /// bins)
  /// @return Size of the grid data structure
  std::size_t size() const override { return m_grid.size(); }

  /// @brief The binning values described by this surface grid lookup
  /// They are in order of the axes
  /// @return Vector of axis directions for binning
  std::vector<AxisDirection> binningValues() const override {
    return m_binValues;
  }

  /// @brief Gets the center position of bin @c bin in global coordinates
  /// @param bin the global bin index
  /// @return The bin center
  Vector3 getBinCenter(std::size_t bin) const override {
    auto gctx = GeometryContext::dangerouslyDefaultConstruct();
    return getBinCenterImpl(gctx, bin);
  }

  /// @brief Returns copies of the axes used in the grid as @c AnyAxis
  /// @return The axes
  /// @note This returns copies. Use for introspection and querying.
  std::vector<const IAxis*> getAxes() const override {
    auto arr = m_grid.axes();
    return std::vector<const IAxis*>(arr.begin(), arr.end());
  }

  std::optional<AnyGridConstView<SurfaceVector>> getGridView() const override {
    return AnyGridConstView<SurfaceVector>{m_grid};
  }

  const Surface* surfaceRepresentation() const override {
    return m_representative.get();
  }

  /// @brief Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  bool isValidBin(std::size_t bin) const override {
    std::array<std::size_t, 2> indices = m_grid.localBinsFromGlobalBin(bin);
    std::array<std::size_t, 2> nBins = m_grid.numLocalBins();
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t idx = indices.at(i);
      if (idx <= 0 || idx >= nBins.at(i) + 1) {
        return false;
      }
    }
    return true;
  }

 private:
  /// map surface center to grid
  std::size_t fillSurfaceToBinMapping(const GeometryContext& gctx,
                                      const Surface& surface) {
    const Vector3 pos = surface.referencePosition(gctx, AxisDirection::AxisR);
    const Vector3 normal = m_representative->normal(gctx, pos);
    const std::size_t globalBin = findGlobalBin(pos, normal, m_tolerance);
    if (globalBin != 0) {
      m_grid.at(globalBin).push_back(&surface);
    }
    return globalBin;
  };

  /// flood fill neighboring bins given a starting bin
  void fillBinToSurfaceMapping(const GeometryContext& gctx,
                               const Surface& surface, std::size_t startBin) {
    const std::array<std::size_t, 2> startIndices =
        m_grid.localBinsFromGlobalBin(startBin);
    const auto startNeighborIndices =
        m_grid.neighborHoodIndices(startIndices, 1u);

    std::set<std::size_t> visited({startBin});
    std::vector<std::size_t> queue(startNeighborIndices.begin(),
                                   startNeighborIndices.end());

    while (!queue.empty()) {
      const std::size_t current = queue.back();
      queue.pop_back();
      if (visited.contains(current)) {
        continue;
      }

      const std::array<std::size_t, 2> currentIndices =
          m_grid.localBinsFromGlobalBin(current);
      visited.insert(current);

      const std::array<double, 2> gridLocal = m_grid.binCenter(currentIndices);
      const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
      const Vector3 normal = m_representative->normal(gctx, surfaceLocal);
      const Vector3 global =
          m_representative->localToGlobal(gctx, surfaceLocal, normal);

      const Intersection3D intersection =
          surface.intersect(gctx, global, normal, BoundaryTolerance::None())
              .closest();
      if (!intersection.isValid() ||
          std::abs(intersection.pathLength()) > m_tolerance) {
        continue;
      }
      m_grid.at(current).push_back(&surface);

      const auto neighborIndices =
          m_grid.neighborHoodIndices(currentIndices, 1u);
      queue.insert(queue.end(), neighborIndices.begin(), neighborIndices.end());
    }
  };

  void populateNeighborCache() {
    // calculate neighbors for every bin and store in map
    for (std::size_t i = 0; i < m_grid.size(); i++) {
      if (!isValidBin(i)) {
        continue;
      }
      const std::array<std::size_t, 2> indices =
          m_grid.localBinsFromGlobalBin(i);
      std::vector<const Surface*>& neighbors = m_neighborMap.at(i);
      neighbors.clear();

      for (std::size_t idx : m_grid.neighborHoodIndices(indices, 1u)) {
        const std::vector<const Surface*>& binContent = m_grid.at(idx);
        std::copy(binContent.begin(), binContent.end(),
                  std::back_inserter(neighbors));
      }

      std::ranges::sort(neighbors);
      auto last = std::ranges::unique(neighbors);
      neighbors.erase(last.begin(), last.end());
      neighbors.shrink_to_fit();
    }
  }

  Vector3 getBinCenterImpl(const GeometryContext& gctx, std::size_t bin) const {
    const std::array<double, 2> gridLocal =
        m_grid.binCenter(m_grid.localBinsFromGlobalBin(bin));
    const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
    return m_representative->localToGlobal(gctx, surfaceLocal);
  }

  const CylinderBounds* getCylinderBounds() const {
    return dynamic_cast<const CylinderBounds*>(&m_representative->bounds());
  }

  Vector2 gridToSurfaceLocal(std::array<double, 2> gridLocal) const {
    Vector2 surfaceLocal = Eigen::Map<Vector2>(gridLocal.data());
    if (const CylinderBounds* bounds = getCylinderBounds(); bounds != nullptr) {
      surfaceLocal[0] *= bounds->get(CylinderBounds::eR);
    }
    return surfaceLocal;
  }
  std::array<double, 2> surfaceToGridLocal(Vector2 local) const {
    std::array<double, 2> gridLocal = {local[0], local[1]};
    if (const CylinderBounds* bounds = getCylinderBounds(); bounds != nullptr) {
      gridLocal[0] /= bounds->get(CylinderBounds::eR);
    }
    return gridLocal;
  }

  std::size_t findGlobalBin(const Vector3& position, const Vector3& direction,
                            double tolerance) const {
    auto gctx = GeometryContext::dangerouslyDefaultConstruct();

    const Intersection3D intersection =
        m_representative
            ->intersect(gctx, position, direction,
                        BoundaryTolerance::Infinite())
            .closest();
    if (!intersection.isValid() ||
        std::abs(intersection.pathLength()) > tolerance) {
      return 0;  // overflow bin
    }
    const Vector2 surfaceLocal =
        m_representative
            ->globalToLocal(gctx, intersection.position(), direction)
            .value();
    const std::array<double, 2> gridLocal = surfaceToGridLocal(surfaceLocal);
    return m_grid.globalBinFromPosition(gridLocal);
  }

  std::shared_ptr<RegularSurface> m_representative;
  double m_tolerance{};
  Grid_t m_grid;
  std::vector<AxisDirection> m_binValues;
  std::vector<SurfaceVector> m_neighborMap;
};

}  // namespace

std::unique_ptr<SurfaceArray::ISurfaceGridLookup>
SurfaceArray::makeSurfaceGridLookup(
    std::shared_ptr<RegularSurface> representative, double tolerance,
    std::tuple<const IAxis*, const IAxis*> axes,
    std::vector<AxisDirection> bValues) {
  const auto& [iAxisA, iAxisB] = axes;

  return iAxisA->visit([&]<typename axis_a_t>(const axis_a_t& axisA) {
    return iAxisB->visit(
        [&]<typename axis_b_t>(const axis_b_t& axisB)
            -> std::unique_ptr<SurfaceArray::ISurfaceGridLookup> {
          return std::make_unique<SurfaceGridLookup2<axis_a_t, axis_b_t>>(
              representative, tolerance,
              std::tuple<axis_a_t, axis_b_t>{axisA, axisB}, bValues);
        });
  });
}

}  // namespace Acts
