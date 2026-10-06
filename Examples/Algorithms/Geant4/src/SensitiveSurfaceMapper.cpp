// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/HashCombine.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TransformHelpers.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsExamples/Geant4/AlgebraConverters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>
#include <utility>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4Polyhedron.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <boost/container/flat_set.hpp>
#include <boost/geometry.hpp>

// Add some type traits for boost::geometry, so we can use the machinery
// directly with Acts::Vector2 / Eigen::Matrix
namespace boost::geometry::traits {

template <typename T, int D>
struct tag<Eigen::Matrix<T, D, 1>> {
  using type = point_tag;
};
template <typename T, int D>
struct dimension<Eigen::Matrix<T, D, 1>> : std::integral_constant<int, D> {};
template <typename T, int D>
struct coordinate_type<Eigen::Matrix<T, D, 1>> {
  using type = T;
};
template <typename T, int D>
struct coordinate_system<Eigen::Matrix<T, D, 1>> {
  using type = boost::geometry::cs::cartesian;
};

template <typename T, int D, std::size_t Index>
struct access<Eigen::Matrix<T, D, 1>, Index> {
  static_assert(Index < D, "Out of range");
  using Point = Eigen::Matrix<T, D, 1>;
  using CoordinateType = typename coordinate_type<Point>::type;
  static inline CoordinateType get(Point const& p) { return p[Index]; }
  static inline void set(Point& p, CoordinateType const& value) {
    p[Index] = value;
  }
};

}  // namespace boost::geometry::traits

namespace {

void writeG4Polyhedron(
    Acts::IVisualization3D& visualizer, const G4Polyhedron& polyhedron,
    const Acts::Transform3& trafo = Acts::Transform3::Identity(),
    Acts::Color color = {0, 0, 0}) {
  constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;

  for (int i = 1; i <= polyhedron.GetNoFacets(); ++i) {
    // This is a bit ugly but I didn't find an easy way to compute this
    // beforehand.
    constexpr std::size_t maxPoints = 1000;
    G4Point3D points[maxPoints];
    int nPoints = 0;
    polyhedron.GetFacet(i, nPoints, points);
    assert(static_cast<std::size_t>(nPoints) < maxPoints);

    std::vector<Acts::Vector3> faces;
    for (int j = 0; j < nPoints; ++j) {
      faces.emplace_back(points[j][0] * convertLength,
                         points[j][1] * convertLength,
                         points[j][2] * convertLength);
      faces.back() = trafo * faces.back();
    }

    visualizer.face(faces, color);
  }
}

}  // namespace

namespace ActsExamples::Geant4 {

SensitiveCandidates::SensitiveCandidates(
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    std::unique_ptr<const Acts::Logger> _logger)
    : m_trackingGeo{trackingGeometry}, m_logger{std::move(_logger)} {}
std::vector<const Acts::Surface*> SensitiveCandidates::queryPosition(
    const Acts::GeometryContext& gctx, const Acts::Vector3& position) const {
  std::vector<const Acts::Surface*> surfaces{};
  ACTS_VERBOSE("Try to fetch the surfaces close to " << position.transpose());

  switch (m_trackingGeo->geometryVersion()) {
    using enum Acts::TrackingGeometry::GeometryVersion;
    case Gen1: {
      // In case we do not find a layer at this position for whatever reason
      const auto layer = m_trackingGeo->associatedLayer(gctx, position);
      if (layer == nullptr) {
        return surfaces;
      }

      const auto surfaceArray = layer->surfaceArray();
      if (surfaceArray == nullptr) {
        return surfaces;
      }

      for (const auto& surface : surfaceArray->surfaces()) {
        if (surface->isSensitive()) {
          surfaces.push_back(surface);
        }
      }
      break;
    }
    case Gen3: {
      const auto* refVolume =
          m_trackingGeo->resolveLowestTrackingVolume(gctx, position).value();
      if (refVolume != nullptr) {
        constexpr bool restrictToSensitives = true;
        refVolume->visitSurfaces(
            [&](const Acts::Surface* surface) { surfaces.push_back(surface); },
            restrictToSensitives);
      }
      break;
    }
  }
  return surfaces;
}
std::optional<const Acts::GeometryObject*> SensitiveCandidatesBase::queryRegion(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::Vector3& /*position*/) const {
  return std::nullopt;
}

std::optional<const Acts::GeometryObject*> SensitiveCandidates::queryRegion(
    const Acts::GeometryContext& gctx, const Acts::Vector3& position) const {
  // Must resolve the same object that `queryPosition` takes the surfaces from
  switch (m_trackingGeo->geometryVersion()) {
    using enum Acts::TrackingGeometry::GeometryVersion;
    case Gen1:
      return m_trackingGeo->associatedLayer(gctx, position);
    case Gen3:
      return m_trackingGeo->resolveLowestTrackingVolume(gctx, position).value();
  }
  return std::nullopt;
}

std::vector<const Acts::Surface*> SensitiveCandidates::queryAll() const {
  std::vector<const Acts::Surface*> surfaces;

  constexpr bool restrictToSensitives = true;
  m_trackingGeo->visitSurfaces(
      [&](auto surface) { surfaces.push_back(surface); }, restrictToSensitives);

  return surfaces;
}

/// Candidate surfaces of one region together with their centers, and a grid
/// over the centers. It finds the first candidate centered at a position
/// without computing and comparing the centers of all candidates.
class SensitiveSurfaceMapper::CandidateIndex {
 public:
  /// @param surfaces the candidate surfaces, in the order they are matched
  /// @param gctx the geometry context to compute the centers
  /// @param buildGrid whether to build the grid; without it the center
  ///        lookup is a linear scan, which is cheaper for a single lookup
  CandidateIndex(std::vector<const Acts::Surface*> surfaces,
                 const Acts::GeometryContext& gctx, bool buildGrid)
      : m_surfaces{std::move(surfaces)} {
    m_centers.reserve(m_surfaces.size());
    for (const auto [i, surface] : Acts::enumerate(m_surfaces)) {
      m_centers.push_back(surface->center(gctx));
      if (surface->bounds().type() == Acts::SurfaceBounds::eAnnulus) {
        m_annulusIndices.push_back(i);
      }
    }
    // Two centers that compare equal must be at most one cell apart
    m_useGrid = buildGrid && m_comparator.compare<3>(
                                 Acts::Vector3::Zero(),
                                 Acts::Vector3{0.5 * s_cellSize, 0., 0.}) != 0;
    if (!m_useGrid) {
      return;
    }
    for (const auto [i, center] : Acts::enumerate(m_centers)) {
      if (const auto cell = cellOf(center); cell.has_value()) {
        m_grid[*cell].push_back(i);
      } else {
        m_unbinned.push_back(i);
      }
    }
  }

  bool empty() const { return m_surfaces.empty(); }
  std::size_t size() const { return m_surfaces.size(); }
  const Acts::Surface& surface(std::size_t i) const { return *m_surfaces[i]; }
  const Acts::Vector3& center(std::size_t i) const { return m_centers[i]; }

  /// Indices of the candidates with annulus bounds, in ascending order
  const std::vector<std::size_t>& annulusIndices() const {
    return m_annulusIndices;
  }

  /// Index of the first candidate whose center compares equal to
  /// @p position, i.e. the one a linear scan in candidate order would find
  std::optional<std::size_t> firstCenterMatch(
      const Acts::Vector3& position) const {
    const auto matches = [&](std::size_t i) {
      return m_comparator.compare<3>(m_centers[i], position) == 0;
    };
    const auto cell = m_useGrid ? cellOf(position) : std::nullopt;
    if (!cell.has_value()) {
      for (std::size_t i = 0; i < m_centers.size(); ++i) {
        if (matches(i)) {
          return i;
        }
      }
      return std::nullopt;
    }
    std::size_t first = m_centers.size();
    const auto check = [&](const std::vector<std::size_t>& indices) {
      for (const std::size_t i : indices) {
        if (i < first && matches(i)) {
          first = i;
        }
      }
    };
    check(m_unbinned);
    for (std::int64_t dx = -1; dx <= 1; ++dx) {
      for (std::int64_t dy = -1; dy <= 1; ++dy) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
          const CellKey key{(*cell)[0] + dx, (*cell)[1] + dy, (*cell)[2] + dz};
          if (const auto it = m_grid.find(key); it != m_grid.end()) {
            check(it->second);
          }
        }
      }
    }
    if (first == m_centers.size()) {
      return std::nullopt;
    }
    return first;
  }

 private:
  using CellKey = std::array<std::int64_t, 3>;
  struct CellHash {
    std::size_t operator()(const CellKey& key) const {
      return Acts::hashMixAndCombine(key[0], key[1], key[2]);
    }
  };

  /// Cell size of the center grid, much larger than the comparator tolerance
  static constexpr double s_cellSize = 1. * Acts::UnitConstants::mm;

  /// Grid cell of a position, std::nullopt if it is not finite or too large
  /// to be binned. Such centers are always compared, such positions use the
  /// linear scan.
  static std::optional<CellKey> cellOf(const Acts::Vector3& position) {
    constexpr double maxCell = 1e15;
    CellKey key{};
    for (std::size_t i = 0; i < 3; ++i) {
      const double cell = std::floor(position[i] / s_cellSize);
      // Also false for NaN
      if (!(std::abs(cell) < maxCell)) {
        return std::nullopt;
      }
      key[i] = static_cast<std::int64_t>(cell);
    }
    return key;
  }

  std::vector<const Acts::Surface*> m_surfaces;
  std::vector<Acts::Vector3> m_centers;
  std::vector<std::size_t> m_annulusIndices;
  Acts::detail::TransformComparator m_comparator{};
  bool m_useGrid{false};
  std::unordered_map<CellKey, std::vector<std::size_t>, CellHash> m_grid;
  std::vector<std::size_t> m_unbinned;
};

/// Candidate lookups shared over one traversal of the Geant4 tree
struct SensitiveSurfaceMapper::Cache {
  /// Per region returned by `SensitiveCandidatesBase::queryRegion`
  std::unordered_map<const Acts::GeometryObject*, CandidateIndex> regions;
  /// All sensitive surfaces, the fallback if no region has candidates
  std::optional<CandidateIndex> all;
};

SensitiveSurfaceMapper::SensitiveSurfaceMapper(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

void SensitiveSurfaceMapper::remapSensitiveNames(
    State& state, const Acts::GeometryContext& gctx,
    G4VPhysicalVolume* g4PhysicalVolume,
    const Acts::Transform3& motherTransform) const {
  Cache cache;
  remapSensitiveNames(state, cache, gctx, g4PhysicalVolume, motherTransform);
}

void SensitiveSurfaceMapper::remapSensitiveNames(
    State& state, Cache& cache, const Acts::GeometryContext& gctx,
    G4VPhysicalVolume* g4PhysicalVolume,
    const Acts::Transform3& motherTransform) const {
  // Make sure the unit conversion is correct

  auto g4LogicalVolume = g4PhysicalVolume->GetLogicalVolume();
  auto g4SensitiveDetector = g4LogicalVolume->GetSensitiveDetector();

  // Get the transform of the G4 object
  Acts::Transform3 localG4ToGlobal{Acts::Transform3::Identity()};
  {
    auto g4Translation = g4PhysicalVolume->GetTranslation();
    auto g4Rotation = g4PhysicalVolume->GetRotation();
    Acts::Vector3 g4RelPosition = convertPosition(g4Translation);
    Acts::Translation3 translation(g4RelPosition);
    if (g4Rotation == nullptr) {
      localG4ToGlobal = motherTransform * translation;
    } else {
      Acts::RotationMatrix3 rotation;
      rotation << g4Rotation->xx(), g4Rotation->yx(), g4Rotation->zx(),
          g4Rotation->xy(), g4Rotation->yy(), g4Rotation->zy(),
          g4Rotation->xz(), g4Rotation->yz(), g4Rotation->zz();
      localG4ToGlobal =
          motherTransform * Acts::makeTransform3(rotation, g4RelPosition);
    }
  }

  const Acts::Vector3 g4AbsPosition = localG4ToGlobal.translation();

  if (G4int nDaughters = g4LogicalVolume->GetNoDaughters(); nDaughters > 0) {
    // Step down to all daughters
    for (G4int id = 0; id < nDaughters; ++id) {
      remapSensitiveNames(state, cache, gctx, g4LogicalVolume->GetDaughter(id),
                          localG4ToGlobal);
    }
  }

  const std::string& volumeName{g4LogicalVolume->GetName()};
  const std::string& volumeMaterialName{
      g4LogicalVolume->GetMaterial()->GetName()};

  const bool isSensitive = g4SensitiveDetector != nullptr;
  const bool isMappedMaterial =
      Acts::rangeContainsValue(m_cfg.materialMappings, volumeMaterialName);
  const bool isMappedVolume =
      Acts::rangeContainsSubstring(m_cfg.volumeMappings, volumeName);

  if (!(isSensitive || isMappedMaterial || isMappedVolume)) {
    ACTS_VERBOSE("Did not try mapping '"
                 << g4PhysicalVolume->GetName() << "' at "
                 << g4AbsPosition.transpose()
                 << " because g4SensitiveDetector (=" << g4SensitiveDetector
                 << ") is null and volume name (=" << volumeName
                 << ") and material name (=" << volumeMaterialName
                 << ") were not found");
    return;
  }
  ACTS_VERBOSE("Attempt to map " << g4PhysicalVolume->GetName() << "' at "
                                 << g4AbsPosition.transpose()
                                 << " to the tracking geometry");

  // Query the candidates at the first polyhedron vertex that has any. They
  // are cached per region, so the lookup is built only once per region.
  const CandidateIndex* candidates = nullptr;
  std::optional<CandidateIndex> uncachedCandidates;
  const auto g4Polyhedron = g4LogicalVolume->GetSolid()->GetPolyhedron();
  for (int i = 1; i < g4Polyhedron->GetNoVertices(); ++i) {
    auto vtx = convertPosition(g4Polyhedron->GetVertex(i));
    auto vtxGlobal = localG4ToGlobal * vtx;

    if (const auto region =
            m_cfg.candidateSurfaces->queryRegion(gctx, vtxGlobal);
        region.has_value()) {
      auto it = cache.regions.find(*region);
      if (it == cache.regions.end()) {
        it = cache.regions
                 .try_emplace(
                     *region,
                     m_cfg.candidateSurfaces->queryPosition(gctx, vtxGlobal),
                     gctx, true)
                 .first;
      }
      candidates = &it->second;
    } else {
      uncachedCandidates.emplace(
          m_cfg.candidateSurfaces->queryPosition(gctx, vtxGlobal), gctx, false);
      candidates = &*uncachedCandidates;
    }

    if (!candidates->empty()) {
      break;
    }
  }

  // Fall back to query all surfaces
  if (candidates == nullptr || candidates->empty()) {
    ACTS_DEBUG("No candidate surfaces for volume '" << volumeName << "' at "
                                                    << g4AbsPosition.transpose()
                                                    << ", query all surfaces");
    if (!cache.all.has_value()) {
      cache.all.emplace(m_cfg.candidateSurfaces->queryAll(), gctx, true);
    }
    candidates = &*cache.all;
  }

  ACTS_VERBOSE("Found " << candidates->size() << " candidate surfaces for "
                        << volumeName);

  // The match is the first candidate that either has its center at the G4
  // position or, for annulus bounds, has its bounds centroid inside the G4
  // solid. So only annulus candidates before the first center match need the
  // centroid check.
  const Acts::Surface* mappedSurface = nullptr;
  const auto centerMatch = candidates->firstCenterMatch(g4AbsPosition);
  const std::size_t nBeforeCenterMatch =
      centerMatch.value_or(candidates->size());
  for (const std::size_t i : candidates->annulusIndices()) {
    if (i >= nBeforeCenterMatch) {
      break;
    }
    const Acts::Surface& candidateSurface = candidates->surface(i);
    const auto& bounds =
        *static_cast<const Acts::AnnulusBounds*>(&candidateSurface.bounds());

    const auto vertices = bounds.vertices(0);

    constexpr bool clockwise = false;
    constexpr bool closed = false;
    using Polygon =
        boost::geometry::model::polygon<Acts::Vector2, clockwise, closed>;

    Polygon poly;
    boost::geometry::assign_points(poly, vertices);

    Acts::Vector2 boundsCentroidSurfaceFrame = Acts::Vector2::Zero();
    boost::geometry::centroid(poly, boundsCentroidSurfaceFrame);

    Acts::Vector3 boundsCentroidGlobal{boundsCentroidSurfaceFrame[0],
                                       boundsCentroidSurfaceFrame[1], 0.0};
    boundsCentroidGlobal =
        candidateSurface.localToGlobalTransform(gctx) * boundsCentroidGlobal;

    const auto boundsCentroidG4Frame =
        localG4ToGlobal.inverse() * boundsCentroidGlobal;

    if (g4LogicalVolume->GetSolid()->Inside(
            convertPosition(boundsCentroidG4Frame)) != EInside::kOutside) {
      ACTS_VERBOSE("Successful match with centroid matching");
      mappedSurface = &candidateSurface;
      break;
    }
  }
  if (mappedSurface == nullptr && centerMatch.has_value()) {
    ACTS_DEBUG("Successful match with center: "
               << candidates->center(*centerMatch).transpose()
               << ", G4-position: " << g4AbsPosition.transpose());
    mappedSurface = &candidates->surface(*centerMatch);
  }

  Acts::detail::TransformComparator trfSorter{};

  if (mappedSurface == nullptr) {
    ACTS_DEBUG("No mapping found for '"
               << volumeName << "' with material '" << volumeMaterialName
               << "' at position " << g4AbsPosition.transpose());
    state.missingVolumes.emplace_back(g4PhysicalVolume, localG4ToGlobal);
    return;
  }

  // A mapped surface was found, a new name will be set that G4PhysVolume
  ACTS_DEBUG("Matched " << volumeName << " to " << mappedSurface->geometryId()
                        << " at position " << g4AbsPosition.transpose());
  // Check if the prefix is not yet assigned
  if (volumeName.find(mappingPrefix) == std::string::npos) {
    // Set the new name
    std::string mappedName = std::string(mappingPrefix) + volumeName;
    g4PhysicalVolume->SetName(mappedName);
  }
  if (state.g4VolumeToSurfaces.find(g4PhysicalVolume) ==
      state.g4VolumeToSurfaces.end()) {
    state.g4VolumeToSurfaces.insert(
        std::make_pair(g4PhysicalVolume, SurfacePosMap_t{trfSorter}));
  }
  // Insert into the multi-map
  if (!state.g4VolumeToSurfaces[g4PhysicalVolume]
           .insert(std::make_pair(g4AbsPosition, mappedSurface))
           .second) {
    ACTS_WARNING("Duplicate surface found for " << volumeName << " @ "
                                                << g4AbsPosition.transpose());
  }
}

bool SensitiveSurfaceMapper::checkMapping(
    const State& state, const Acts::GeometryContext& gctx,
    bool writeMissingG4VolsAsObj, bool writeMissingSurfacesAsObj) const {
  auto allSurfaces = m_cfg.candidateSurfaces->queryAll();
  std::ranges::sort(allSurfaces);

  std::vector<const Acts::Surface*> found;
  for (const auto& [_, surfaceMap] : state.g4VolumeToSurfaces) {
    for (const auto& [__, surfacePtr] : surfaceMap) {
      found.push_back(surfacePtr);
    }
  }
  std::ranges::sort(found);
  auto newEnd = std::unique(found.begin(), found.end());
  found.erase(newEnd, found.end());

  std::vector<const Acts::Surface*> missing;
  std::set_difference(allSurfaces.begin(), allSurfaces.end(), found.begin(),
                      found.end(), std::back_inserter(missing));

  ACTS_INFO("Number of overall sensitive surfaces: " << allSurfaces.size());
  ACTS_INFO("Number of mapped volume->surface mappings: " << found.size());
  ACTS_INFO(
      "Number of sensitive surfaces that are not mapped: " << missing.size());
  ACTS_INFO("Number of G4 volumes without a matching Surface: "
            << state.missingVolumes.size());

  if (writeMissingG4VolsAsObj) {
    Acts::ObjVisualization3D visualizer;
    for (const auto& [g4vol, trafo] : state.missingVolumes) {
      auto polyhedron = g4vol->GetLogicalVolume()->GetSolid()->GetPolyhedron();
      writeG4Polyhedron(visualizer, *polyhedron, trafo);
    }

    std::ofstream os("missing_g4_volumes.obj");
    visualizer.write(os);
  }

  if (writeMissingSurfacesAsObj) {
    Acts::ObjVisualization3D visualizer;
    Acts::ViewConfig vcfg;
    vcfg.quarterSegments = 720;
    for (auto srf : missing) {
      Acts::GeometryView3D::drawSurface(visualizer, *srf, gctx,
                                        Acts::Transform3::Identity(), vcfg);
    }

    std::ofstream os("missing_acts_surfaces.obj");
    visualizer.write(os);
  }

  return missing.empty();
}

}  // namespace ActsExamples::Geant4
