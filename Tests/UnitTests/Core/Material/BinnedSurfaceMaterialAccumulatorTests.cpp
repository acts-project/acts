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
#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "Acts/Utilities/detail/WorkerPool.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

using namespace Acts;

namespace ActsTests {

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

namespace {

/// Accumulate a single track, without worker pool
void accumulateTrack(
    const ISurfaceMaterialAccumulator& accumulator,
    ISurfaceMaterialAccumulator::State& state,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>& empty) {
  ISurfaceMaterialAccumulator::TrackAccumulation track{&interactions, &empty};
  accumulator.accumulateTracks(state, tContext, std::span(&track, 1), nullptr);
}

/// The original per-track accumulation, kept as the reference: average the
/// first touched bin of every surface, then apply the empty bin correction
void referenceAccumulate(
    BinnedSurfaceMaterialAccumulator::State& state, bool emptyBinCorrection,
    const std::vector<MaterialInteraction>& interactions,
    const std::vector<IAssignmentFinder::SurfaceAssignment>& empty) {
  std::map<AccumulatedSurfaceMaterial*, std::array<std::size_t, 3>> touched;
  for (const auto& mi : interactions) {
    auto& accMaterial = state.accumulatedMaterial.at(mi.surface->geometryId());
    auto bin = accMaterial.accumulate(mi.intersection, mi.materialSlab,
                                      mi.pathCorrection);
    touched.emplace(&accMaterial, bin);
  }
  for (const auto& [accMaterial, bin] : touched) {
    accMaterial->trackAverage(std::vector<std::array<std::size_t, 3>>{bin},
                              true);
  }
  if (emptyBinCorrection) {
    for (const auto& missed : empty) {
      state.accumulatedMaterial.at(missed.surface->geometryId())
          .trackAverage(missed.position, true);
    }
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(MaterialSuite)

BOOST_AUTO_TEST_CASE(InvalidSetupTest) {
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
  };

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // First is homogeneous
  MaterialSlab mp = MaterialSlab::Nothing();
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(mp, 1.));

  // Second is empty - invalid

  BinnedSurfaceMaterialAccumulator::Config bsmaConfig;
  bsmaConfig.materialSurfaces = {surfaces[0].get(), surfaces[1].get()};
  bsmaConfig.emptyBinCorrection = true;

  BinnedSurfaceMaterialAccumulator bsma(
      bsmaConfig,
      getDefaultLogger("BinnedSurfaceMaterialAccumulator", Logging::VERBOSE));

  // Generate the state - this throws because the second surface has no
  // material assigned.
  BOOST_CHECK_THROW(bsma.createState(tContext), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(AccumulationTest) {
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 50.0,
                                           100.0)};

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // Accepted materials are:
  // - homogeneous
  // - Prot0
  // - Binned (remapping)

  // First is homogeneous
  MaterialSlab mp = MaterialSlab::Nothing();
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(mp, 1.));

  // Second surface is binned Phi / Z
  BinUtility sb1(4, -std::numbers::pi, std::numbers::pi, closed,
                 AxisDirection::AxisPhi);
  sb1 += BinUtility(2, -100., 100., open, AxisDirection::AxisZ);
  surfaces[1u]->assignSurfaceMaterial(
      std::make_shared<ProtoSurfaceMaterial>(sb1));

  // Third is binned
  std::vector<MaterialSlab> mps = {mp, mp, mp};
  BinUtility sb2(3, -100., 100., open, AxisDirection::AxisZ);
  surfaces[2u]->assignSurfaceMaterial(
      std::make_shared<BinnedSurfaceMaterial>(sb2, mps));

  BinnedSurfaceMaterialAccumulator::Config bsmaConfig;
  bsmaConfig.materialSurfaces = {surfaces[0].get(), surfaces[1].get(),
                                 surfaces[2].get()};
  bsmaConfig.emptyBinCorrection = true;

  BinnedSurfaceMaterialAccumulator bsma(
      bsmaConfig,
      getDefaultLogger("BinnedSurfaceMaterialAccumulator", Logging::VERBOSE));

  // Generate the state
  auto state = bsma.createState(tContext);

  auto cState =
      static_cast<const BinnedSurfaceMaterialAccumulator::State*>(state.get());

  BOOST_CHECK_EQUAL(cState->accumulatedMaterial.size(), 3u);

  // Intersections
  // Track 0:
  // - Surface 0 hit once
  // - Surface 1 hit twice
  // - Surface 2 empty hit
  Vector3 p(0, 0, 0);
  Vector3 d0 = Vector3(50, 1, 0).normalized();

  MaterialInteraction m00;
  m00.surface = surfaces[0u].get();
  m00.position = 20 * d0;
  m00.direction = d0;

  MaterialInteraction m01;
  m01.surface = surfaces[1u].get();
  m01.position = 30 * d0;
  m01.direction = d0;

  MaterialInteraction m02;
  m02.surface = surfaces[1u].get();
  m02.position = 30 * d0;
  m02.direction = d0;

  std::vector<MaterialInteraction> mInteractions = {m01, m01, m02};
  std::vector<IAssignmentFinder::SurfaceAssignment> emptyHits = {
      {surfaces[2u].get(), 50 * d0, d0}};

  accumulateTrack(bsma, *state, mInteractions, emptyHits);

  // Track 1:
  // - Surface 0 empty hit
  // - Surface 1 hit once
  // - Surface 2 hit once
  Vector3 d1 = Vector3(10, 10, 0).normalized();
  MaterialInteraction m11;
  m11.surface = surfaces[1u].get();
  m11.position = 30 * d1;
  m11.direction = d1;

  MaterialInteraction m12;
  m12.surface = surfaces[2u].get();
  m12.position = 50 * d1;
  m12.direction = d1;

  mInteractions = {m11, m12};
  emptyHits = {{surfaces[0u].get(), 50 * d1, d1}};
  accumulateTrack(bsma, *state, mInteractions, emptyHits);

  // Get the maps
  auto maps = bsma.finalizeMaterial(*state, tContext);

  BOOST_CHECK_EQUAL(maps.size(), 3u);

  auto m0Itr = maps.find(GeometryIdentifier().withSensitive(1));
  BOOST_CHECK(m0Itr != maps.end());
  BOOST_CHECK(m0Itr->second != nullptr);

  auto m1Itr = maps.find(GeometryIdentifier().withSensitive(2));
  BOOST_CHECK(m1Itr != maps.end());
  BOOST_CHECK(m1Itr->second != nullptr);

  auto m2Itr = maps.find(GeometryIdentifier().withSensitive(3));
  BOOST_CHECK(m2Itr != maps.end());
  BOOST_CHECK(m2Itr->second != nullptr);

  // Check the material
  // map0 should be homogeneous
  auto m0 = m0Itr->second;
  const HomogeneousSurfaceMaterial* hm0 =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(m0.get());
  BOOST_CHECK(hm0 != nullptr);

  // map1 should be binned
  auto m1 = m1Itr->second;
  const BinnedSurfaceMaterial* bm1 =
      dynamic_cast<const BinnedSurfaceMaterial*>(m1.get());
  BOOST_CHECK(bm1 != nullptr);

  // map2 should be binned
  auto m2 = m2Itr->second;
  const BinnedSurfaceMaterial* bm2 =
      dynamic_cast<const BinnedSurfaceMaterial*>(m2.get());
  BOOST_CHECK(bm2 != nullptr);

  // Check failures
  auto invalidSurface =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 40.0, 100.0);
  invalidSurface->assignGeometryId(GeometryIdentifier().withSensitive(4));

  // Invalid surface amongst material
  MaterialInteraction mXX;
  mXX.surface = invalidSurface.get();
  mXX.position = 50 * d1;
  mXX.direction = d1;
  BOOST_CHECK_THROW(accumulateTrack(bsma, *state, {mXX}, {}),
                    std::invalid_argument);

  // Invalid surface amongst empty hits
  BOOST_CHECK_THROW(
      accumulateTrack(bsma, *state, {}, {{invalidSurface.get(), 50 * d1, d1}}),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ProtoGridResolutionTest) {
  using enum AxisDirection;

  // A full cylinder and a full disc carrying deferred binning specs
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30.0, 80.0)};

  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }

  // Cylinder: deferred equidistant binning in (rphi, z), given in swapped
  // order to exercise the direction based re-ordering
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<ProtoGridSurfaceMaterial>(
          MultiAxisSpec2D({AxisSpec::DeferredEquidistant(10, AxisZ),
                           AxisSpec::DeferredEquidistant(4, AxisRPhi)})));

  // Disc: deferred variable binning in r, deferred equidistant binning in phi
  surfaces[1u]->assignSurfaceMaterial(
      std::make_shared<ProtoGridSurfaceMaterial>(MultiAxisSpec2D(
          {AxisSpec::DeferredVariable({0., 0.2, 1.}, std::nullopt, AxisR),
           AxisSpec::DeferredEquidistant(8, AxisPhi)})));

  BinnedSurfaceMaterialAccumulator::Config bsmaConfig;
  bsmaConfig.materialSurfaces = {surfaces[0].get(), surfaces[1].get()};

  BinnedSurfaceMaterialAccumulator bsma(
      bsmaConfig,
      getDefaultLogger("BinnedSurfaceMaterialAccumulator", Logging::VERBOSE));

  auto state = bsma.createState(tContext);
  auto cState =
      static_cast<const BinnedSurfaceMaterialAccumulator::State*>(state.get());
  BOOST_CHECK_EQUAL(cState->accumulatedMaterial.size(), 2u);

  // Cylinder binning: rphi first (re-ordered), closed on the full azimuth
  const auto& cylBinning =
      cState->accumulatedMaterial.at(GeometryIdentifier().withSensitive(1))
          .binUtility()
          .binningData();
  BOOST_CHECK_EQUAL(cylBinning.size(), 2u);
  BOOST_CHECK_EQUAL(cylBinning[0].binvalue, AxisRPhi);
  BOOST_CHECK_EQUAL(cylBinning[0].bins(), 4u);
  BOOST_CHECK_EQUAL(cylBinning[0].option, closed);
  BOOST_CHECK_CLOSE(cylBinning[0].max, 20. * std::numbers::pi, 1e-4);
  BOOST_CHECK_EQUAL(cylBinning[1].binvalue, AxisZ);
  BOOST_CHECK_EQUAL(cylBinning[1].bins(), 10u);
  BOOST_CHECK_EQUAL(cylBinning[1].option, open);
  BOOST_CHECK_EQUAL(cylBinning[1].min, -100.);
  BOOST_CHECK_EQUAL(cylBinning[1].max, 100.);

  // Disc binning: variable r edges scaled onto [30, 80], closed phi
  const auto& discBinning =
      cState->accumulatedMaterial.at(GeometryIdentifier().withSensitive(2))
          .binUtility()
          .binningData();
  BOOST_CHECK_EQUAL(discBinning.size(), 2u);
  BOOST_CHECK_EQUAL(discBinning[0].binvalue, AxisR);
  BOOST_CHECK_EQUAL(discBinning[0].bins(), 2u);
  BOOST_CHECK_EQUAL(discBinning[0].min, 30.);
  BOOST_CHECK_EQUAL(discBinning[0].max, 80.);
  BOOST_CHECK_CLOSE(discBinning[0].boundaries()[1], 40., 1e-4);
  BOOST_CHECK_EQUAL(discBinning[1].binvalue, AxisPhi);
  BOOST_CHECK_EQUAL(discBinning[1].bins(), 8u);
  BOOST_CHECK_EQUAL(discBinning[1].option, closed);
}

BOOST_AUTO_TEST_CASE(BatchAccumulationMatchesSerial) {
  // Homogeneous, proto-binned and binned surfaces
  std::vector<std::shared_ptr<Surface>> surfaces = {
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 20.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30.0, 100.0),
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 50.0,
                                           100.0)};
  for (auto [is, surface] : enumerate(surfaces)) {
    surface->assignGeometryId(GeometryIdentifier().withSensitive(is + 1));
  }
  surfaces[0u]->assignSurfaceMaterial(
      std::make_shared<HomogeneousSurfaceMaterial>(MaterialSlab::Nothing()));
  BinUtility sb1(8, -std::numbers::pi, std::numbers::pi, closed,
                 AxisDirection::AxisPhi);
  sb1 += BinUtility(4, -100., 100., open, AxisDirection::AxisZ);
  surfaces[1u]->assignSurfaceMaterial(
      std::make_shared<ProtoSurfaceMaterial>(sb1));
  BinUtility sb2(5, -100., 100., open, AxisDirection::AxisZ);
  surfaces[2u]->assignSurfaceMaterial(std::make_shared<BinnedSurfaceMaterial>(
      sb2, std::vector<MaterialSlab>(5, MaterialSlab::Nothing())));

  // Random tracks: each surface is hit (possibly several times, in several
  // bins), missed with an empty hit, or not touched at all
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uniform(0., 1.);
  struct Track {
    std::vector<MaterialInteraction> interactions;
    std::vector<IAssignmentFinder::SurfaceAssignment> empty;
  };
  std::vector<Track> tracks(500);
  for (auto& track : tracks) {
    for (const auto& surface : surfaces) {
      const double r = surface->bounds().values()[0];
      const double phi = std::numbers::pi * (2 * uniform(rng) - 1);
      const double z = 190 * uniform(rng) - 95;
      const Vector3 position(r * std::cos(phi), r * std::sin(phi), z);
      const double roll = uniform(rng);
      if (roll < 0.6) {
        const int nHits = 1 + static_cast<int>(3 * uniform(rng));
        for (int h = 0; h < nHits; ++h) {
          MaterialInteraction mi;
          mi.surface = surface.get();
          // Later hits may land in a different bin of the same surface
          mi.intersection =
              h == 0 ? position
                     : Vector3(r * std::cos(phi + h), r * std::sin(phi + h), z);
          mi.pathCorrection = 1. + uniform(rng);
          mi.materialSlab =
              MaterialSlab(Material::fromMassDensity(
                               1 + 100 * uniform(rng), 1 + 300 * uniform(rng),
                               1 + 50 * uniform(rng), 1 + 20 * uniform(rng),
                               0.1 + 10 * uniform(rng)),
                           0.1 + uniform(rng));
          track.interactions.push_back(mi);
        }
      } else if (roll < 0.85) {
        track.empty.push_back({surface.get(), position, position.normalized()});
      }
    }
    std::ranges::shuffle(track.interactions, rng);
  }

  auto checkEqual = [](const auto& lhs, const auto& rhs) {
    BOOST_REQUIRE_EQUAL(lhs.size(), rhs.size());
    for (const auto& [geoID, material] : lhs) {
      const auto& other = rhs.at(geoID);
      if (auto* hl =
              dynamic_cast<const HomogeneousSurfaceMaterial*>(material.get())) {
        auto* hr = dynamic_cast<const HomogeneousSurfaceMaterial*>(other.get());
        BOOST_REQUIRE(hr != nullptr);
        BOOST_CHECK(*hl == *hr);
      } else {
        auto* bl = dynamic_cast<const BinnedSurfaceMaterial*>(material.get());
        auto* br = dynamic_cast<const BinnedSurfaceMaterial*>(other.get());
        BOOST_REQUIRE(bl != nullptr && br != nullptr);
        BOOST_CHECK(bl->fullMaterial() == br->fullMaterial());
      }
    }
  };

  for (bool emptyBinCorrection : {true, false}) {
    BinnedSurfaceMaterialAccumulator::Config cfg;
    cfg.materialSurfaces = {surfaces[0].get(), surfaces[1].get(),
                            surfaces[2].get()};
    cfg.emptyBinCorrection = emptyBinCorrection;
    BinnedSurfaceMaterialAccumulator bsma(cfg);

    // Reference: the original per-track accumulation
    auto serial = bsma.createState(tContext);
    for (const auto& track : tracks) {
      referenceAccumulate(
          static_cast<BinnedSurfaceMaterialAccumulator::State&>(*serial),
          emptyBinCorrection, track.interactions, track.empty);
    }
    auto reference = bsma.finalizeMaterial(*serial, tContext);

    // One track at a time
    auto single = bsma.createState(tContext);
    for (const auto& track : tracks) {
      accumulateTrack(bsma, *single, track.interactions, track.empty);
    }
    checkEqual(reference, bsma.finalizeMaterial(*single, tContext));

    // Batches, without pool and with pools of different size
    for (std::size_t nThreads : {0u, 1u, 4u}) {
      std::optional<detail::WorkerPool> pool;
      if (nThreads > 0) {
        pool.emplace(nThreads);
      }
      auto batched = bsma.createState(tContext);
      for (std::size_t begin = 0; begin < tracks.size(); begin += 128) {
        std::vector<ISurfaceMaterialAccumulator::TrackAccumulation> batch;
        for (std::size_t t = begin; t < std::min(begin + 128, tracks.size());
             ++t) {
          batch.push_back({&tracks[t].interactions, &tracks[t].empty});
        }
        bsma.accumulateTracks(*batched, tContext, batch,
                              pool.has_value() ? &*pool : nullptr);
      }
      checkEqual(reference, bsma.finalizeMaterial(*batched, tContext));
    }
  }

  // An unknown surface is reported from the worker threads
  auto unknown =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 40.0, 100.0);
  unknown->assignGeometryId(GeometryIdentifier().withSensitive(4));
  MaterialInteraction mi;
  mi.surface = unknown.get();
  std::vector<MaterialInteraction> interactions = {mi};
  std::vector<IAssignmentFinder::SurfaceAssignment> none;
  BinnedSurfaceMaterialAccumulator::Config cfg;
  cfg.materialSurfaces = {surfaces[0].get()};
  BinnedSurfaceMaterialAccumulator bsma(cfg);
  auto state = bsma.createState(tContext);
  detail::WorkerPool pool(4);
  std::vector<ISurfaceMaterialAccumulator::TrackAccumulation> batch(
      8, {&interactions, &none});
  BOOST_CHECK_THROW(bsma.accumulateTracks(*state, tContext, batch, &pool),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
