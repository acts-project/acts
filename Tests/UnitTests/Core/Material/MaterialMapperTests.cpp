// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IAssignmentFinder.hpp"
#include "Acts/Material/ISurfaceMaterialAccumulator.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <map>
#include <numbers>
#include <random>
#include <span>
#include <stdexcept>

using namespace Acts;

namespace ActsTests {

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

/// @brief Interface for the material mapping that seeks the possible
/// assignment candidates for the material interactiosn
class IntersectSurfacesFinder : public IAssignmentFinder {
 public:
  std::vector<const Surface*> surfaces;

  /// @brief Interface method for generating assignment candidates for the
  /// material interaction assignment to surfaces or volumes
  ///
  /// @param gctx is the geometry context
  /// @param mctx is the magnetic field context
  /// @param position is the position of the initial ray
  /// @param direction is the direction of initial ray
  ///
  /// @return a vector of Surface Assignments and Volume Assignments
  std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
            std::vector<IAssignmentFinder::VolumeAssignment>>
  assignmentCandidates(const GeometryContext& gctx,
                       const MagneticFieldContext& /*ignored*/,
                       const Vector3& position,
                       const Vector3& direction) const override {
    std::vector<IAssignmentFinder::SurfaceAssignment> surfaceAssignments;
    std::vector<IAssignmentFinder::VolumeAssignment> volumeAssignments;
    // Intersect the surfaces
    for (auto& surface : surfaces) {
      // Get the intersection
      MultiIntersection3D multiIntersection = surface->intersect(
          gctx, position, direction, BoundaryTolerance::None());
      // One solution, take it
      if (multiIntersection.size() == 1u &&
          multiIntersection.at(0).status() >= IntersectionStatus::reachable &&
          multiIntersection.at(0).pathLength() >= 0.0) {
        surfaceAssignments.push_back(
            {surface, multiIntersection.at(0).position(), direction});
        continue;
      }
      if (multiIntersection.size() > 1u) {
        // Multiple intersections, take the closest
        Intersection3D closestForward = multiIntersection.closestForward();
        if (closestForward.status() >= IntersectionStatus::reachable &&
            closestForward.pathLength() > 0.0) {
          surfaceAssignments.push_back(
              {surface, closestForward.position(), direction});
          continue;
        }
      }
    }
    return {surfaceAssignments, volumeAssignments};
  }
};

/// @brief Interface for the material mapping, this is the accumulation step
class MaterialBlender : public ISurfaceMaterialAccumulator {
 public:
  explicit MaterialBlender(
      const std::vector<std::shared_ptr<Surface>>& surfaces = {})
      : m_surfaces(surfaces) {}

  /// The state of the material accumulator, this is used
  /// to cache information across tracks/events
  class State final : public ISurfaceMaterialAccumulator::State {
   public:
    std::map<const Surface*, AccumulatedMaterialSlab> accumulatedMaterial;
  };

  /// Factory for creating the state
  std::unique_ptr<ISurfaceMaterialAccumulator::State> createState(
      const GeometryContext& /*gctx*/) const override {
    auto state = std::make_unique<State>();
    for (auto& surface : m_surfaces) {
      state->accumulatedMaterial[surface.get()] = AccumulatedMaterialSlab();
    }
    return state;
  };

  /// @brief Accumulate the material interaction on the surface
  ///
  /// @param state is the state of the accumulator
  /// @param interactions is the material interactions, with assigned surfaces
  /// @param surfacesWithoutAssignment are the surfaces without assignment
  ///
  /// @note this the track average over the binned material
  void accumulateTracks(ISurfaceMaterialAccumulator::State& state,
                        const GeometryContext& /*gctx*/,
                        std::span<const TrackAccumulation> tracks,
                        detail::WorkerPool* /*pool*/) const override {
    for (const auto& track : tracks) {
      accumulateTrack(state, *track.interactions);
    }
  }

  void accumulateTrack(
      ISurfaceMaterialAccumulator::State& state,
      const std::vector<MaterialInteraction>& interactions) const {
    auto cState = static_cast<State*>(&state);
    for (const auto& mi : interactions) {
      // Get the surface
      const Surface* surface = mi.surface;
      // Get the accumulated material
      auto accMaterial = cState->accumulatedMaterial.find(surface);
      if (accMaterial == cState->accumulatedMaterial.end()) {
        throw std::invalid_argument(
            "Surface material is not found, inconsistent configuration.");
      }
      // Accumulate the material
      accMaterial->second.accumulate(mi.materialSlab);
    }
    // Average over the track
    for (auto& [surface, accumulatedMaterial] : cState->accumulatedMaterial) {
      accumulatedMaterial.trackAverage();
    }
  }

  /// Finalize the surface material maps
  ///
  /// @param state the state of the accumulator
  ///
  /// @note this does the run average over the (binned) material
  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(ISurfaceMaterialAccumulator::State& state,
                   const GeometryContext& /*gctx*/) const override {
    auto cState = static_cast<State*>(&state);

    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
        materialMaps;
    for (auto& [surface, accumulatedMaterial] : cState->accumulatedMaterial) {
      materialMaps[surface->geometryId()] =
          std::make_shared<HomogeneousSurfaceMaterial>(
              accumulatedMaterial.totalAverage().first);
    }
    return materialMaps;
  }

 private:
  std::vector<std::shared_ptr<Surface>> m_surfaces;
};

BOOST_AUTO_TEST_SUITE(MaterialSuite)

/// @brief This test checks the data flow of the material mapper, it is not
/// a test of the single components, which are tested individually
///
/// @note Currently only surface mapping is implemented, volume mapping/assigning
/// is being added in future PRs
BOOST_AUTO_TEST_CASE(MaterialMapperFlowTest) {
  // Create a vector of surfaces
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 50.0,
                                           100.0)};

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // The assigner
  auto assigner = std::make_shared<IntersectSurfacesFinder>();
  assigner->surfaces = {surfaces[0].get(), surfaces[1].get(),
                        surfaces[2].get()};

  // The accumulator - which blends all the material
  auto accumulator = std::make_shared<MaterialBlender>(surfaces);

  // Create the mapper
  MaterialMapper::Config mmConfig;
  mmConfig.assignmentFinder = assigner;
  mmConfig.surfaceMaterialAccumulator = accumulator;

  MaterialMapper mapper(mmConfig);

  auto state = mapper.createState(tContext);
  BOOST_CHECK(state.get() != nullptr);
  BOOST_CHECK(state->surfaceMaterialAccumulatorState.get() != nullptr);

  std::vector<RecordedMaterialTrack> mappedTracks;
  std::vector<RecordedMaterialTrack> unmappedTracks;

  // Track loop
  Vector3 position(0., 0., 0.);
  for (unsigned int it = 0; it < 11; ++it) {
    Vector3 direction =
        Vector3(0.9 + it * 0.02, 1.1 - it * 0.02, 0.).normalized();
    RecordedMaterialTrack mTrack{{position, direction}, {}};
    for (unsigned int im = 0; im < 60; ++im) {
      MaterialInteraction mi;
      mi.materialSlab = MaterialSlab(
          Material::fromMassDensity(it + 1, it + 1, it + 1, it + 1, it + 1),
          0.1);
      mi.position = position + (im + 1) * direction;
      mi.direction = direction;
      mTrack.second.materialInteractions.push_back(mi);
    }
    auto [mapped, unmapped] = mapper.mapMaterial(*state, tContext, {}, mTrack);
    mappedTracks.push_back(mapped);
    unmappedTracks.push_back(unmapped);
  }

  // Get the maps
  auto maps = mapper.finalizeMaps(*state, tContext);

  BOOST_CHECK(maps.surfaceMaterials.size() == 3);
  BOOST_CHECK(maps.volumeMaterials.empty());
}

BOOST_AUTO_TEST_CASE(MaterialMapperBatchMatchesSerial) {
  std::vector<std::shared_ptr<Surface>> surfaces;
  for (double r : {20., 30., 50., 80.}) {
    surfaces.push_back(
        Surface::makeShared<CylinderSurface>(Transform3::Identity(), r, 100.));
  }
  BinUtility binning(6, -std::numbers::pi, std::numbers::pi, closed,
                     AxisDirection::AxisPhi);
  binning += BinUtility(5, -100., 100., open, AxisDirection::AxisZ);
  std::vector<const Surface*> surfacePtrs;
  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
    surface->assignSurfaceMaterial(
        std::make_shared<ProtoSurfaceMaterial>(binning));
    surfacePtrs.push_back(surface.get());
  }

  IntersectionMaterialAssigner::Config assignerConfig;
  assignerConfig.surfaces = surfacePtrs;
  BinnedSurfaceMaterialAccumulator::Config accumulatorConfig;
  accumulatorConfig.materialSurfaces = surfacePtrs;

  MaterialMapper::Config mmConfig;
  mmConfig.assignmentFinder =
      std::make_shared<IntersectionMaterialAssigner>(assignerConfig);
  mmConfig.surfaceMaterialAccumulator =
      std::make_shared<BinnedSurfaceMaterialAccumulator>(accumulatorConfig);
  MaterialMapper serialMapper(mmConfig);
  mmConfig.numThreads = 4;
  MaterialMapper parallelMapper(mmConfig);

  // Straight tracks from the origin with interactions along the way, some of
  // them leaving through the cylinder ends
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> uniform(0., 1.);
  std::vector<RecordedMaterialTrack> tracks(300);
  for (auto& [startDir, recorded] : tracks) {
    const double phi = std::numbers::pi * (2 * uniform(rng) - 1);
    const double cotTheta = 2 * uniform(rng) - 1;
    startDir = {Vector3::Zero(),
                Vector3(std::cos(phi), std::sin(phi), cotTheta).normalized()};
    for (double s = 1.; s < 120.; s += 0.5 + 4 * uniform(rng)) {
      MaterialInteraction mi;
      mi.position = s * startDir.second;
      mi.direction = startDir.second;
      mi.materialSlab =
          MaterialSlab(Material::fromMassDensity(1 + 100 * uniform(rng),
                                                 1 + 300 * uniform(rng), 10.,
                                                 5., 1 + uniform(rng)),
                       0.1 + uniform(rng));
      recorded.materialInteractions.push_back(mi);
    }
  }

  auto serialState = serialMapper.createState(tContext);
  std::vector<std::pair<RecordedMaterialTrack, RecordedMaterialTrack>>
      serialResults;
  for (const auto& track : tracks) {
    serialResults.push_back(
        serialMapper.mapMaterial(*serialState, tContext, {}, track));
  }

  auto parallelState = parallelMapper.createState(tContext);
  std::vector<std::pair<RecordedMaterialTrack, RecordedMaterialTrack>>
      parallelResults;
  for (std::size_t begin = 0; begin < tracks.size(); begin += 64) {
    std::vector<const RecordedMaterialTrack*> batch;
    for (std::size_t t = begin; t < std::min(begin + 64, tracks.size()); ++t) {
      batch.push_back(&tracks[t]);
    }
    for (auto& result :
         parallelMapper.mapMaterial(*parallelState, tContext, {}, batch)) {
      parallelResults.push_back(std::move(result));
    }
  }

  // Identical mapped and unmapped tracks
  BOOST_REQUIRE_EQUAL(serialResults.size(), parallelResults.size());
  std::size_t nMapped = 0;
  for (std::size_t t = 0; t < serialResults.size(); ++t) {
    for (auto [serial, parallel] :
         {std::pair{&serialResults[t].first, &parallelResults[t].first},
          std::pair{&serialResults[t].second, &parallelResults[t].second}}) {
      const auto& sInteractions = serial->second.materialInteractions;
      const auto& pInteractions = parallel->second.materialInteractions;
      BOOST_REQUIRE_EQUAL(sInteractions.size(), pInteractions.size());
      for (std::size_t i = 0; i < sInteractions.size(); ++i) {
        BOOST_CHECK_EQUAL(sInteractions[i].surface, pInteractions[i].surface);
        BOOST_CHECK(sInteractions[i].intersection ==
                    pInteractions[i].intersection);
      }
      BOOST_CHECK_EQUAL(serial->second.materialInX0,
                        parallel->second.materialInX0);
    }
    nMapped += serialResults[t].first.second.materialInteractions.size();
  }
  BOOST_CHECK_GT(nMapped, 0u);

  // Identical maps
  auto serialMaps = serialMapper.finalizeMaps(*serialState, tContext);
  auto parallelMaps = parallelMapper.finalizeMaps(*parallelState, tContext);
  BOOST_REQUIRE_EQUAL(serialMaps.surfaceMaterials.size(), surfaces.size());
  BOOST_REQUIRE_EQUAL(parallelMaps.surfaceMaterials.size(), surfaces.size());
  for (const auto& [key, material] : serialMaps.surfaceMaterials) {
    auto* sm = dynamic_cast<const BinnedSurfaceMaterial*>(material.get());
    auto* pm = dynamic_cast<const BinnedSurfaceMaterial*>(
        parallelMaps.surfaceMaterials.at(key).get());
    BOOST_REQUIRE(sm != nullptr && pm != nullptr);
    BOOST_CHECK(sm->fullMaterial() == pm->fullMaterial());
  }
}

/// Counts the calls of both accumulation methods, overriding only one
template <bool overrideBatch, bool overrideSingle>
class CallCounter : public ISurfaceMaterialAccumulator {
 public:
  mutable std::size_t batchCalls = 0;
  mutable std::size_t singleCalls = 0;

  std::unique_ptr<ISurfaceMaterialAccumulator::State> createState(
      const GeometryContext& /*gctx*/) const override {
    return std::make_unique<ISurfaceMaterialAccumulator::State>();
  }

  void accumulateTracks(ISurfaceMaterialAccumulator::State& state,
                        const GeometryContext& gctx,
                        std::span<const TrackAccumulation> tracks,
                        detail::WorkerPool* pool) const override {
    if constexpr (overrideBatch) {
      batchCalls += tracks.size();
    } else {
      ISurfaceMaterialAccumulator::accumulateTracks(state, gctx, tracks, pool);
    }
  }

  ACTS_PUSH_IGNORE_DEPRECATED()
  void accumulate(ISurfaceMaterialAccumulator::State& state,
                  const GeometryContext& gctx,
                  const std::vector<MaterialInteraction>& interactions,
                  const std::vector<IAssignmentFinder::SurfaceAssignment>&
                      surfacesWithoutAssignment) const override {
    if constexpr (overrideSingle) {
      ++singleCalls;
    } else {
      ISurfaceMaterialAccumulator::accumulate(state, gctx, interactions,
                                              surfacesWithoutAssignment);
    }
  }
  ACTS_POP_IGNORE_DEPRECATED()

  std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>
  finalizeMaterial(ISurfaceMaterialAccumulator::State& /*state*/,
                   const GeometryContext& /*gctx*/) const override {
    return {};
  }
};

BOOST_AUTO_TEST_CASE(DeprecatedSingleTrackAccumulation) {
  std::vector<MaterialInteraction> interactions;
  std::vector<IAssignmentFinder::SurfaceAssignment> empty;
  std::vector<ISurfaceMaterialAccumulator::TrackAccumulation> tracks(
      3, {&interactions, &empty});

  // An implementation of the deprecated method serves batch calls
  CallCounter<false, true> legacy;
  auto legacyState = legacy.createState(tContext);
  legacy.accumulateTracks(*legacyState, tContext, tracks, nullptr);
  BOOST_CHECK_EQUAL(legacy.singleCalls, 3u);

  // The deprecated method does not forward to the batch method
  CallCounter<true, false> modern;
  auto modernState = modern.createState(tContext);
  ACTS_PUSH_IGNORE_DEPRECATED()
  BOOST_CHECK_THROW(
      modern.accumulate(*modernState, tContext, interactions, empty),
      std::logic_error);
  ACTS_POP_IGNORE_DEPRECATED()
  modern.accumulateTracks(*modernState, tContext, tracks, nullptr);
  BOOST_CHECK_EQUAL(modern.batchCalls, 3u);

  // Overriding neither is an error
  CallCounter<false, false> broken;
  auto brokenState = broken.createState(tContext);
  BOOST_CHECK_THROW(
      broken.accumulateTracks(*brokenState, tContext, tracks, nullptr),
      std::logic_error);
}

BOOST_AUTO_TEST_CASE(MaterialMapperInvalidTest) {
  // The assigner
  auto assigner = std::make_shared<IntersectSurfacesFinder>();

  // The accumulator - which blends all the material
  auto accumulator = std::make_shared<MaterialBlender>();

  // Create the mapper
  MaterialMapper::Config mmConfigInvalid;
  mmConfigInvalid.assignmentFinder = nullptr;
  mmConfigInvalid.surfaceMaterialAccumulator = accumulator;

  BOOST_CHECK_THROW(auto mapperIA = MaterialMapper(mmConfigInvalid),
                    std::invalid_argument);

  // BOOST_CHECK_THROW(MaterialMapper(mmConfigInvalid), std::invalid_argument);

  mmConfigInvalid.assignmentFinder = assigner;
  mmConfigInvalid.surfaceMaterialAccumulator = nullptr;

  BOOST_CHECK_THROW(auto mapperIS = MaterialMapper(mmConfigInvalid),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
