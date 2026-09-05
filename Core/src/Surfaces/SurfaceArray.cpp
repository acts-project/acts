// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Ranges.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/detail/MultiAxisHelper.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <ranges>
#include <set>
#include <stdexcept>
#include <utility>

namespace Acts {

/// Base interface for all surface lookups.
struct SurfaceArray::ISurfaceGridLookup {
  virtual ~ISurfaceGridLookup() = default;

  /// Fill provided surfaces into the contained @c Grid.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param surfaces Input surface pointers
  virtual void fill(const GeometryContext& gctx,
                    std::span<const Surface* const> surfaces) = 0;

  /// Get all surfaces in bin given by the global bin index
  /// @param bin the global bin index
  /// @return span of surface pointers of the bin at that position
  virtual std::span<const Surface* const> at(std::size_t bin) const = 0;

  /// Performs lookup at @c pos and returns bin content as const reference
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return A span of surface pointers
  virtual std::span<const Surface* const> at(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const = 0;

  /// Get all surfaces in bin given by local grid indices and neighbor
  /// distance.
  /// @param gridIndices the local grid indices
  /// @param neighborDistance the neighbor distance to include in the lookup,
  ///        per axis
  /// @return span of surface pointers of the bin at that position and its neighbors
  virtual std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::array<std::uint8_t, 2> neighborDistance) const = 0;

  /// Performs a lookup at @c pos, but returns neighbors as well
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Lookup position
  /// @param direction Lookup direction
  /// @return A span of surface pointers
  virtual std::span<const Surface* const> neighbors(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const = 0;

  /// Returns the total size of the grid (including under/overflow bins)
  /// @return Size of the grid data structure
  virtual std::size_t size() const = 0;

  /// Gets the center position of bin @c bin in global coordinates
  /// @param bin the global bin index
  /// @return The bin center
  virtual Vector3 getBinCenter(std::size_t bin) const = 0;

  /// Returns copies of the axes used in the grid as @c AnyAxis
  /// @return The axes
  /// @note This returns copies. Use for introspection and querying.
  virtual std::vector<const IAxis*> getAxes() const = 0;

  /// Get the representative surface used for this lookup
  /// @return Surface pointer
  virtual const Surface* surfaceRepresentation() const = 0;

  /// Checks if global bin is valid
  /// @param bin the global bin index
  /// @return bool if the bin is valid
  /// @note Valid means that the index points to a bin which is not a under
  ///       or overflow bin or out of range in any axis.
  virtual bool isValidBin(std::size_t bin) const = 0;

  /// The binning values described by this surface grid lookup. They are in
  /// order of the axes (optional) and empty for eingle lookups
  /// @return Vector of axis directions for binning
  virtual std::vector<AxisDirection> binningValues() const { return {}; }

  /// Get the number of local bins in each dimension. This is used to
  /// determine the size of the grid for neighbor lookups.
  /// @return Array of number of local bins in each dimension
  virtual std::array<std::size_t, 2> numLocalBins() const = 0;

  /// Get the bounds on the neighbor window this lookup serves. The window
  /// itself is derived per axis from the crossing angle and clamped to them.
  /// @return Neighbor window bounds per axis
  virtual SurfaceArray::NeighborWindow neighborWindow() const = 0;
};

namespace {

struct SingleElementLookupImpl final : SurfaceArray::ISurfaceGridLookup {
  explicit SingleElementLookupImpl(const Surface* element)
      : m_element({element}) {}

  std::span<const Surface* const> at(std::size_t bin) const override {
    if (bin != 0) {
      throw std::out_of_range(
          "SingleElementLookupImpl only contains one bin with index 0");
    }
    return m_element;
  }

  std::span<const Surface* const> at(
      const GeometryContext& /*gctx*/, const Vector3& /*position*/,
      const Vector3& /*direction*/) const override {
    return m_element;
  }

  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::array<std::uint8_t, 2> neighborDistance) const override {
    if (gridIndices != std::array<std::size_t, 2>{0, 0} ||
        neighborDistance != std::array<std::uint8_t, 2>{0, 0}) {
      throw std::out_of_range(
          "SingleElementLookupImpl only contains one bin with zero neighbor "
          "distance");
    }
    return m_element;
  }

  std::span<const Surface* const> neighbors(
      const GeometryContext& /*gctx*/, const Vector3& /*position*/,
      const Vector3& /*direction*/) const override {
    return m_element;
  }

  std::size_t size() const override { return 1; }

  Vector3 getBinCenter(std::size_t /*bin*/) const override {
    return Vector3(0, 0, 0);
  }

  std::vector<const IAxis*> getAxes() const override { return {}; }

  const Surface* surfaceRepresentation() const override { return nullptr; }

  void fill(const GeometryContext& /*gctx*/,
            std::span<const Surface* const> /*surfaces*/) override {
    // no-op: the single element is already fixed at construction time
  }

  bool isValidBin(std::size_t bin) const override { return bin == 0; }

  std::array<std::size_t, 2> numLocalBins() const override { return {1, 1}; }

  SurfaceArray::NeighborWindow neighborWindow() const override {
    return {.min = {0, 0}, .max = {0, 0}};
  }

 private:
  std::vector<const Surface*> m_element;
};

template <class Axis1, class Axis2>
struct SurfaceGridLookupImpl final : SurfaceArray::ISurfaceGridLookup {
  SurfaceGridLookupImpl(std::shared_ptr<RegularSurface> representative,
                        double tolerance, std::tuple<Axis1, Axis2> axes,
                        std::vector<AxisDirection> binValues = {},
                        SurfaceArray::NeighborWindow neighborWindow = {})
      : m_representative(std::move(representative)),
        m_tolerance(tolerance),
        m_axes(std::move(axes)),
        m_binValues(std::move(binValues)),
        m_neighborWindow(neighborWindow) {
    if (m_neighborWindow.min[0] > m_neighborWindow.max[0] ||
        m_neighborWindow.min[1] > m_neighborWindow.max[1]) {
      throw std::invalid_argument("neighbor window floor exceeds its bound");
    }
    m_fillingGrid.resize(size());
  }

  void fill(const GeometryContext& gctx,
            std::span<const Surface* const> surfaces) override {
    for (const Surface* surface : surfaces) {
      fillSurfaceFootprint(gctx, *surface);
    }

    for (std::vector<const Surface*>& binSurfaces : m_fillingGrid) {
      std::ranges::sort(binSurfaces);
      const auto last = std::ranges::unique(binSurfaces);
      binSurfaces.erase(last.begin(), last.end());
      binSurfaces.shrink_to_fit();
    }

    checkGrid(surfaces);

    populateNeighborCache();
  }

  std::span<const Surface* const> at(std::size_t globalBin) const override {
    return m_fillingGrid.at(globalBin);
  }

  std::span<const Surface* const> at(const GeometryContext& gctx,
                                     const Vector3& position,
                                     const Vector3& direction) const override {
    const std::optional<GridIndex> localBins = findLocalBin2D(
        gctx, position, direction, std::numeric_limits<double>::infinity());
    if (!localBins.has_value()) {
      return {};
    }
    return surfacePack(neighborPackIndex(*localBins, {0, 0}));
  }

  std::span<const Surface* const> neighbors(
      std::array<std::size_t, 2> gridIndices,
      std::array<std::uint8_t, 2> neighborDistance) const override {
    // beyond the bound there is no cached pack, and the index would wrap into
    // the next bin rather than run off the end
    if (neighborDistance[0] > m_neighborWindow.max[0] ||
        neighborDistance[1] > m_neighborWindow.max[1]) {
      throw std::out_of_range(
          "neighbor distance exceeds the neighbor window bound");
    }
    // local bins run from the underflow bin 0 to the overflow bin nBins + 1
    const GridIndex nBins = numLocalBins2D();
    if (gridIndices[0] > nBins[0] + 1 || gridIndices[1] > nBins[1] + 1) {
      throw std::out_of_range("local bin is outside the grid");
    }
    return surfacePack(neighborPackIndex(gridIndices, neighborDistance));
  }

  std::span<const Surface* const> neighbors(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const override {
    const std::optional<Vector2> surfaceLocal = findSurfaceLocal(
        gctx, position, direction, std::numeric_limits<double>::infinity());
    if (!surfaceLocal.has_value()) {
      return {};
    }

    const GridPoint gridLocal = surfaceToGridLocal(*surfaceLocal);
    const GridIndex localBins = localBinsFromPosition2D(gridLocal);

    const GridDistance neighborDistance =
        crossingNeighborDistance(gctx, *surfaceLocal, direction, localBins);

    return surfacePack(neighborPackIndex(localBins, neighborDistance));
  }

  std::size_t size() const override {
    const GridIndex nBins = numLocalBins2D();
    return (nBins[0] + 2) * (nBins[1] + 2);
  }

  std::vector<AxisDirection> binningValues() const override {
    return m_binValues;
  }

  Vector3 getBinCenter(std::size_t bin) const override {
    const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
    const GridPoint gridLocal = binCenter(localBinsFromGlobalBin2D(bin));
    const Vector2 surfaceLocal = gridToSurfaceLocal(gridLocal);
    return m_representative->localToGlobal(gctx, surfaceLocal);
  }

  std::vector<const IAxis*> getAxes() const override {
    return {&std::get<0>(m_axes), &std::get<1>(m_axes)};
  }

  const Surface* surfaceRepresentation() const override {
    return m_representative.get();
  }

  bool isValidBin(std::size_t globalBin) const override {
    const GridIndex indices = localBinsFromGlobalBin2D(globalBin);
    return isValidBin(indices);
  }

  std::array<std::size_t, 2> numLocalBins() const override {
    return numLocalBins2D();
  }

  SurfaceArray::NeighborWindow neighborWindow() const override {
    return m_neighborWindow;
  }

 private:
  using GridIndex = std::array<std::size_t, 2>;
  using GridPoint = std::array<double, 2>;
  /// Neighbor window half width in bins, per axis
  using GridDistance = std::array<std::uint8_t, 2>;
  /// Offset into @c m_surfacePacks and number of surfaces
  using SurfacePackRange = std::pair<std::uint32_t, std::uint32_t>;

  std::shared_ptr<RegularSurface> m_representative;
  double m_tolerance{};
  // needs to be a tuple for the grid_helper functions
  std::tuple<Axis1, Axis2> m_axes;
  std::vector<AxisDirection> m_binValues;
  SurfaceArray::NeighborWindow m_neighborWindow{};

  // legacy grid for filling and for deprecated lookup methods.
  // TODO: remove this once deprecated lookup methods are removed and filling is
  // done directly into the neighbor cache
  std::vector<std::vector<const Surface*>> m_fillingGrid;

  // containers to store the surfaces in the custom grid. the packs are indexed
  // per (bin, distance along axis 0, distance along axis 1), so the index array
  // grows with the square of the maximum distance - it holds ranges rather than
  // spans to keep that affordable.
  std::vector<const Surface*> m_surfacePacks;
  std::vector<SurfacePackRange> m_neighborSurfacePacks;

  std::span<const Surface* const> surfacePack(std::size_t packIndex) const {
    const SurfacePackRange& range = m_neighborSurfacePacks.at(packIndex);
    return {m_surfacePacks.data() + range.first, range.second};
  }

  bool isValidBin(const GridIndex& indices) const {
    const GridIndex nBins = numLocalBins2D();
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const std::size_t idx = indices.at(i);
      if (idx <= 0 || idx >= nBins.at(i) + 1) {
        return false;
      }
    }
    return true;
  }

  GridIndex numLocalBins2D() const {
    return {std::get<0>(m_axes).getNBins(), std::get<1>(m_axes).getNBins()};
  }

  GridIndex localBinsFromPosition2D(const GridPoint& point) const {
    return detail::MultiAxisHelper::getLocalBinsFromPoint(point, m_axes);
  }

  GridIndex localBinsFromGlobalBin2D(std::size_t globalBin) const {
    return detail::MultiAxisHelper::getLocalBinsFromGlobalBin(globalBin,
                                                              m_axes);
  }

  std::size_t globalBinFromLocalBins2D(const GridIndex& localBins) const {
    return detail::MultiAxisHelper::getGlobalBinFromLocalBins(localBins,
                                                              m_axes);
  }

  std::size_t neighborDistanceStride(std::size_t axis) const {
    return static_cast<std::size_t>(m_neighborWindow.max.at(axis)) + 1;
  }

  std::size_t neighborPackIndex(const GridIndex& localBins,
                                const GridDistance& neighborDistance) const {
    const std::size_t globalGridBin =
        detail::MultiAxisHelper::getGlobalBinFromLocalBins(localBins, m_axes);
    return (globalGridBin * neighborDistanceStride(0) + neighborDistance[0]) *
               neighborDistanceStride(1) +
           neighborDistance[1];
  }

  GridPoint binCenter(const GridIndex& localBins) const {
    return detail::MultiAxisHelper::getBinCenter(localBins, m_axes);
  }

  /// Orthogonal projection of a global point onto the representative surface,
  /// in grid coordinates.
  std::optional<GridPoint> projectToGrid(const GeometryContext& gctx,
                                         const Vector3& position) const {
    const Vector3 normal = m_representative->normal(gctx, position);
    const std::optional<Vector2> surfaceLocal = findSurfaceLocal(
        gctx, position, normal, std::numeric_limits<double>::infinity());
    if (!surfaceLocal.has_value()) {
      return std::nullopt;
    }
    return surfaceToGridLocal(*surfaceLocal);
  }

  /// Whether the axis wraps, so bin distances have to be taken the short way
  /// around.
  bool axisIsClosed(std::size_t axis) const {
    return axis == 0 ? std::get<0>(m_axes).getBoundaryType() ==
                           AxisBoundaryType::Closed
                     : std::get<1>(m_axes).getBoundaryType() ==
                           AxisBoundaryType::Closed;
  }

  /// Distance in bins between two local bins along one axis. Under- and
  /// overflow are clamped to the edge bin, which is as far as a window can
  /// usefully reach.
  std::uint8_t axisBinDistance(std::size_t axis, std::size_t from,
                               std::size_t to) const {
    const std::size_t nBins = numLocalBins2D().at(axis);
    const std::size_t a = std::clamp<std::size_t>(from, 1, nBins);
    const std::size_t b = std::clamp<std::size_t>(to, 1, nBins);
    std::size_t distance = a > b ? a - b : b - a;
    if (axisIsClosed(axis)) {
      distance = std::min(distance, nBins - distance);
    }
    return clampValue<std::uint8_t>(distance);
  }

  /// Below this the track is treated as running along the layer and the window
  /// is opened all the way. Any smaller incidence needs more bins than the
  /// cache can hold anyway.
  static constexpr double s_minIncidence = 1e-4;

  /// How many bins the track can move along each axis while it is inside the
  /// layer.
  ///
  /// The lookup happens where the track crosses the representative surface, but
  /// a surface is registered where it projects onto that surface. In between,
  /// the track traverses the layer thickness and slides along the layer by the
  /// crossing angle - on a barrel almost entirely in z. The window has to cover
  /// the bins that displacement spans, and per axis, because widening the other
  /// one only multiplies the candidate count. The result is clamped to the
  /// configured window bounds.
  GridDistance crossingNeighborDistance(const GeometryContext& gctx,
                                        const Vector2& surfaceLocal,
                                        const Vector3& direction,
                                        const GridIndex& localBins) const {
    const GridDistance maximum = m_neighborWindow.max;

    const Vector3 normal = m_representative->normal(gctx, surfaceLocal);
    const double incidence = std::abs(normal.dot(direction));
    if (incidence < s_minIncidence) {
      return maximum;
    }

    const Vector3 crossing =
        m_representative->localToGlobal(gctx, surfaceLocal);
    const double halfPath = m_tolerance / incidence;

    GridDistance neighborDistance{};
    for (const double side : {-1., 1.}) {
      const std::optional<GridPoint> edge =
          projectToGrid(gctx, crossing + side * halfPath * direction);
      if (!edge.has_value()) {
        return maximum;
      }
      const GridIndex edgeBins = localBinsFromPosition2D(*edge);
      for (std::size_t axis = 0; axis < neighborDistance.size(); ++axis) {
        neighborDistance.at(axis) = std::max(
            neighborDistance.at(axis),
            axisBinDistance(axis, localBins.at(axis), edgeBins.at(axis)));
      }
    }

    return {std::clamp(neighborDistance[0], m_neighborWindow.min[0],
                       m_neighborWindow.max[0]),
            std::clamp(neighborDistance[1], m_neighborWindow.min[1],
                       m_neighborWindow.max[1])};
  }

  /// Register a surface in every bin its projection onto the representative
  /// surface overlaps. The projected outline gives the bins it passes through,
  /// the interior is filled per column, which is exact for a convex outline.
  void fillSurfaceFootprint(const GeometryContext& gctx,
                            const Surface& surface) {
    // the surface has to sit within the layer this grid represents
    const Vector3 reference =
        surface.referencePosition(gctx, AxisDirection::AxisR);
    const Vector3 referenceNormal = m_representative->normal(gctx, reference);
    if (!findSurfaceLocal(gctx, reference, referenceNormal, m_tolerance)
             .has_value()) {
      return;
    }

    // resolution of the outline: segments per quarter circle for curved
    // bounds, and samples per segment, since a straight edge does not project
    // to a straight line in grid coordinates
    constexpr unsigned int nSamples = 32;

    const Polyhedron polyhedron =
        surface.polyhedronRepresentation(gctx, nSamples);

    // columns are keyed by the closed axis, if any, so that the span fill runs
    // along an axis whose bins do not wrap
    const bool firstAxisClosed =
        std::get<0>(m_axes).getBoundaryType() == AxisBoundaryType::Closed;
    const std::size_t iColumn = firstAxisClosed ? 0 : 1;
    const std::size_t iSpan = firstAxisClosed ? 1 : 0;

    // column bin -> [min, max] bin along the other axis
    std::map<std::size_t, std::pair<std::size_t, std::size_t>> columns;
    const auto addSample = [&](const Vector3& point) {
      const std::optional<GridPoint> gridLocal = projectToGrid(gctx, point);
      if (!gridLocal.has_value()) {
        return;
      }
      const GridIndex indices = localBinsFromPosition2D(*gridLocal);
      if (!isValidBin(indices)) {
        return;
      }
      const auto [it, inserted] =
          columns.try_emplace(indices[iColumn], indices[iSpan], indices[iSpan]);
      if (!inserted) {
        it->second.first = std::min(it->second.first, indices[iSpan]);
        it->second.second = std::max(it->second.second, indices[iSpan]);
      }
    };

    for (const Polyhedron::FaceType& face : polyhedron.faces) {
      for (std::size_t i = 0; i < face.size(); ++i) {
        const Vector3& from = polyhedron.vertices.at(face.at(i));
        const Vector3& to =
            polyhedron.vertices.at(face.at((i + 1) % face.size()));
        for (std::size_t k = 0; k < nSamples; ++k) {
          const double t = static_cast<double>(k) / nSamples;
          addSample(from + t * (to - from));
        }
      }
    }

    for (const auto& [column, span] : columns) {
      for (std::size_t i = span.first; i <= span.second; ++i) {
        GridIndex indices{};
        indices[iColumn] = column;
        indices[iSpan] = i;
        m_fillingGrid.at(globalBinFromLocalBins2D(indices)).push_back(&surface);
      }
    }
  }

  /// calculate neighbors for every bin and window size and store in map
  void populateNeighborCache() {
    m_surfacePacks.clear();
    m_neighborSurfacePacks.clear();

    m_neighborSurfacePacks.assign(
        size() * neighborDistanceStride(0) * neighborDistanceStride(1), {0, 0});

    std::vector<const Surface*> surfacePack;
    std::map<std::vector<const Surface*>, SurfacePackRange> surfacesToPackRange;
    for (std::size_t inputGlobalBin = 0; inputGlobalBin < m_fillingGrid.size();
         ++inputGlobalBin) {
      const GridIndex indices = localBinsFromGlobalBin2D(inputGlobalBin);

      if (!isValidBin(indices)) {
        continue;
      }

      for (std::uint8_t distance0 = 0; distance0 <= m_neighborWindow.max[0];
           ++distance0) {
        for (std::uint8_t distance1 = 0; distance1 <= m_neighborWindow.max[1];
             ++distance1) {
          surfacePack.clear();

          const auto span0 = std::get<0>(m_axes).neighborHoodIndices(
              indices[0], std::pair<int, int>{-distance0, distance0});
          const auto span1 = std::get<1>(m_axes).neighborHoodIndices(
              indices[1], std::pair<int, int>{-distance1, distance1});
          for (const std::size_t bin0 : span0) {
            for (const std::size_t bin1 : span1) {
              const std::vector<const Surface*>& binContent =
                  m_fillingGrid.at(globalBinFromLocalBins2D({bin0, bin1}));
              std::copy(binContent.begin(), binContent.end(),
                        std::back_inserter(surfacePack));
            }
          }

          std::ranges::sort(surfacePack);
          const auto last = std::ranges::unique(surfacePack);
          surfacePack.erase(last.begin(), last.end());

          const std::size_t packIndex =
              neighborPackIndex(indices, {distance0, distance1});

          if (const auto it = surfacesToPackRange.find(surfacePack);
              it != surfacesToPackRange.end()) {
            m_neighborSurfacePacks[packIndex] = it->second;
          } else {
            throw_assert(m_surfacePacks.size() + surfacePack.size() <=
                             std::numeric_limits<std::uint32_t>::max(),
                         "surface pack storage exceeds the 32 bit index range");
            const SurfacePackRange surfacePackRange = {
                static_cast<std::uint32_t>(m_surfacePacks.size()),
                static_cast<std::uint32_t>(surfacePack.size())};
            m_surfacePacks.insert(m_surfacePacks.end(), surfacePack.begin(),
                                  surfacePack.end());
            surfacesToPackRange[surfacePack] = surfacePackRange;
            m_neighborSurfacePacks[packIndex] = surfacePackRange;
          }
        }
      }
    }

    m_surfacePacks.shrink_to_fit();
  }

  void checkGrid(std::span<const Surface* const> surfaces) {
    const std::set<const Surface*> allSurfaces(surfaces.begin(),
                                               surfaces.end());

    std::set<const Surface*> seenSurface;
    for (std::size_t globalBin = 0; globalBin < m_fillingGrid.size();
         ++globalBin) {
      for (const Surface* surface : m_fillingGrid.at(globalBin)) {
        seenSurface.insert(surface);
      }
    }

    if (allSurfaces != seenSurface) {
      std::set<const Surface*> diff;
      std::ranges::set_difference(allSurfaces, seenSurface,
                                  std::inserter(diff, diff.begin()));

      throw std::logic_error(std::format(
          "SurfaceArray grid does not contain all surfaces provided! "
          "{} surfaces not seen",
          diff.size()));
    }
  }

  const CylinderBounds* getCylinderBounds() const {
    return dynamic_cast<const CylinderBounds*>(&m_representative->bounds());
  }

  Vector2 gridToSurfaceLocal(const GridPoint& gridLocal) const {
    Vector2 surfaceLocal = {gridLocal[0], gridLocal[1]};
    if (const CylinderBounds* bounds = getCylinderBounds(); bounds != nullptr) {
      surfaceLocal[0] *= bounds->get(CylinderBounds::eR);
    }
    return surfaceLocal;
  }

  GridPoint surfaceToGridLocal(const Vector2& local) const {
    GridPoint gridLocal = {local[0], local[1]};
    if (const CylinderBounds* bounds = getCylinderBounds(); bounds != nullptr) {
      gridLocal[0] /= bounds->get(CylinderBounds::eR);
    }
    return gridLocal;
  }

  std::optional<Vector2> findSurfaceLocal(const GeometryContext& gctx,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          double tolerance) const {
    const Intersection3D intersection =
        m_representative
            ->intersect(gctx, position, direction,
                        BoundaryTolerance::Infinite())
            .closest();
    if (!intersection.isValid() ||
        std::abs(intersection.pathLength()) > tolerance) {
      return std::nullopt;
    }
    const Vector2 surfaceLocal =
        m_representative
            ->globalToLocal(gctx, intersection.position(), direction)
            .value();
    return surfaceLocal;
  }

  std::optional<GridIndex> findLocalBin2D(const GeometryContext& gctx,
                                          const Vector3& position,
                                          const Vector3& direction,
                                          double tolerance) const {
    const std::optional<Vector2> surfaceLocal =
        findSurfaceLocal(gctx, position, direction, tolerance);
    if (!surfaceLocal.has_value()) {
      return std::nullopt;
    }
    const GridPoint gridLocal = surfaceToGridLocal(*surfaceLocal);
    return localBinsFromPosition2D(gridLocal);
  }
};

std::unique_ptr<SurfaceArray::ISurfaceGridLookup> makeSurfaceGridLookup(
    std::shared_ptr<RegularSurface> representative, double tolerance,
    std::tuple<const IAxis&, const IAxis&> axes,
    SurfaceArray::NeighborWindow neighborWindow) {
  const auto& [iAxisA, iAxisB] = axes;

  return iAxisA.visit([&]<typename axis_a_t>(const axis_a_t& axisA) {
    return iAxisB.visit(
        [&]<typename axis_b_t>(const axis_b_t& axisB)
            -> std::unique_ptr<SurfaceArray::ISurfaceGridLookup> {
          return std::make_unique<SurfaceGridLookupImpl<axis_a_t, axis_b_t>>(
              std::move(representative), tolerance,
              std::tuple<axis_a_t, axis_b_t>{axisA, axisB},
              std::vector<AxisDirection>(), neighborWindow);
        });
  });
}

}  // namespace

SurfaceArray::SurfaceArray(std::shared_ptr<const Surface> srf)
    : m_gridLookup(std::make_unique<SingleElementLookupImpl>(srf.get())),
      m_surfaces({std::move(srf)}) {
  m_surfacesRawPointers.push_back(m_surfaces.at(0).get());
}

SurfaceArray::SurfaceArray(const GeometryContext& gctx,
                           std::vector<std::shared_ptr<const Surface>> surfaces,
                           std::shared_ptr<RegularSurface> representative,
                           double tolerance,
                           std::tuple<const IAxis&, const IAxis&> axes,
                           NeighborWindow neighborWindow) {
  m_gridLookup = makeSurfaceGridLookup(std::move(representative), tolerance,
                                       axes, neighborWindow);
  m_surfaces = std::move(surfaces);
  m_surfacesRawPointers =
      m_surfaces |
      std::views::transform(
          [](const std::shared_ptr<const Surface>& sp) { return sp.get(); }) |
      Ranges::to<std::vector>;
  m_gridLookup->fill(gctx, m_surfacesRawPointers);
}

SurfaceArray::SurfaceArray(SurfaceArray&& other) noexcept = default;

SurfaceArray& SurfaceArray::operator=(SurfaceArray&& other) noexcept = default;

SurfaceArray::~SurfaceArray() = default;

std::span<const Surface* const> SurfaceArray::at(std::size_t bin) const {
  return m_gridLookup->at(bin);
}

std::span<const Surface* const> SurfaceArray::at(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  return m_gridLookup->at(gctx, position, direction);
}

std::span<const Surface* const> SurfaceArray::neighbors(
    std::array<std::size_t, 2> gridIndices,
    std::uint8_t neighborDistance) const {
  return m_gridLookup->neighbors(gridIndices,
                                 {neighborDistance, neighborDistance});
}

std::span<const Surface* const> SurfaceArray::neighbors(
    std::array<std::size_t, 2> gridIndices,
    std::array<std::uint8_t, 2> neighborDistance) const {
  return m_gridLookup->neighbors(gridIndices, neighborDistance);
}

std::span<const Surface* const> SurfaceArray::neighbors(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction) const {
  return m_gridLookup->neighbors(gctx, position, direction);
}

std::size_t SurfaceArray::size() const {
  return m_gridLookup->size();
}

Vector3 SurfaceArray::getBinCenter(std::size_t bin) const {
  return m_gridLookup->getBinCenter(bin);
}

std::vector<const IAxis*> SurfaceArray::getAxes() const {
  return m_gridLookup->getAxes();
}

bool SurfaceArray::isValidBin(std::size_t bin) const {
  return m_gridLookup->isValidBin(bin);
}

std::vector<AxisDirection> SurfaceArray::binningValues() const {
  return m_gridLookup->binningValues();
}

std::ostream& SurfaceArray::toStream(const GeometryContext& /*gctx*/,
                                     std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(4);
  sl << "SurfaceArray:" << std::endl;
  sl << " - no surfaces: " << m_surfaces.size() << std::endl;

  const std::vector<const IAxis*> axes = m_gridLookup->getAxes();

  for (const auto [j, axis] : enumerate(axes)) {
    const AxisBoundaryType bdt = axis->getBoundaryType();
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
    sl << "   - type: " << (axis->isEquidistant() ? "equidistant" : "variable")
       << std::endl;
    sl << "   - n bins: " << axis->getNBins() << std::endl;
    sl << "   - bin edges: [ ";
    const std::vector<double> binEdges = axis->getBinEdges();
    for (const auto [i, binEdge] : enumerate(binEdges)) {
      if (i > 0) {
        sl << ", ";
      }
      // Do not display negative zeroes
      sl << ((std::abs(binEdge) >= 5e-4) ? binEdge : 0.0);
    }
    sl << " ]" << std::endl;
  }
  return sl;
}

const Surface* SurfaceArray::surfaceRepresentation() const {
  return m_gridLookup->surfaceRepresentation();
}

std::array<std::size_t, 2> SurfaceArray::numLocalBins() const {
  return m_gridLookup->numLocalBins();
}

SurfaceArray::NeighborWindow SurfaceArray::neighborWindow() const {
  return m_gridLookup->neighborWindow();
}

}  // namespace Acts
