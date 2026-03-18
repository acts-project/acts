// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/detail/fwd.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/FreeTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/TrackParamsLookupAccumulator.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace Acts;

namespace ActsTests {

auto gctx = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(TrackFindingSuite)

using Axis = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
using AxisGen = GridAxisGenerators::EqOpenEqOpen;

using CellBound = std::pair<std::shared_ptr<BoundTrackParameters>,
                            std::shared_ptr<BoundTrackParameters>>;

using GridBound = Grid<CellBound, Axis, Axis>;
using AccBound = TrackParamsLookupAccumulator<GridBound>;

using CellFree = std::pair<std::shared_ptr<FreeTrackParameters>,
                           std::shared_ptr<FreeTrackParameters>>;

using GridFree = Grid<CellFree, Axis, Axis>;
using AccFree = TrackParamsLookupAccumulator<GridFree>;

AxisGen axisGen{{-1, 1}, 2, {-1, 1}, 2};

BOOST_AUTO_TEST_CASE(Exceptions) {
  // Instantiate grid
  GridBound grid(axisGen());
  AccBound acc(grid);

  // Create a reference surface for bound parameters
  auto transform = Transform3::Identity();
  auto bounds1 = std::make_shared<RectangleBounds>(1, 1);
  auto bounds2 = std::make_shared<RectangleBounds>(2, 2);

  auto surf1 = Surface::makeShared<PlaneSurface>(transform, bounds1);

  auto surf2 = Surface::makeShared<PlaneSurface>(transform, bounds2);

  // Create parameters to accumulate
  Vector4 pos{1, 2, 0, 4};
  Vector3 dir{1, 0, 0};
  double P = 1.;

  auto hypothesis1 = ParticleHypothesis::electron();
  auto hypothesis2 = ParticleHypothesis::muon();

  auto pars1 = BoundTrackParameters::create(gctx, surf1, pos, dir, 1. / P,
                                            std::nullopt, hypothesis1)
                   .value();

  auto pars2 = BoundTrackParameters::create(gctx, surf2, pos, dir, 1. / P,
                                            std::nullopt, hypothesis1)
                   .value();

  auto pars3 = BoundTrackParameters::create(gctx, surf1, pos, dir, 1. / P,
                                            std::nullopt, hypothesis2)
                   .value();

  auto pars4 = BoundTrackParameters::create(gctx, surf1, pos, dir, -1. / P,
                                            std::nullopt, hypothesis2)
                   .value();

  // Get the point of the grid
  auto bin = grid.localBinsFromGlobalBin(2);
  auto center = grid.binCenter(bin);
  Vector2 loc{center.at(0), center.at(1)};

  // Fill in grid
  acc.addTrack(pars1, pars1, loc);

  // Different reference surfaces
  BOOST_CHECK_THROW(acc.addTrack(pars2, pars2, loc), std::invalid_argument);

  // Different particle hypotheses
  BOOST_CHECK_THROW(acc.addTrack(pars3, pars3, loc), std::invalid_argument);

  // Different charges
  BOOST_CHECK_THROW(acc.addTrack(pars4, pars4, loc), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Accumulation) {
  // Instantiate grids
  GridBound gridBound(axisGen());
  AccBound accBound(gridBound);

  GridFree gridFree(axisGen());
  AccFree accFree(gridFree);

  // Create a reference surface for bound parameters
  auto transform = Transform3::Identity();
  auto bounds = std::make_shared<RectangleBounds>(1, 1);
  auto surf = Surface::makeShared<PlaneSurface>(transform, bounds);

  auto hypothesis = ParticleHypothesis::electron();

  std::vector<Vector4> avgPoss;
  std::vector<Vector3> avgMoms;
  Vector4 pos{1, 2, 0, 4};
  for (std::size_t i = 0; i < gridBound.size(); i++) {
    // Create parameters to accumulate
    std::array<Vector4, 4> fourPositions = {pos * (i + 1), pos * (i + 2),
                                            pos * (i + 3), pos * (i + 4)};

    std::array<double, 4> thetas = {
        std::numbers::pi / (i + 1), std::numbers::pi / (i + 2),
        std::numbers::pi / (i + 3), std::numbers::pi / (i + 4)};

    std::array<double, 4> phis = {
        2 * std::numbers::pi / (i + 1), 2 * std::numbers::pi / (i + 2),
        2 * std::numbers::pi / (i + 3), 2 * std::numbers::pi / (i + 4)};

    double P = 1.5 * (i + 1);

    // Get the point of the grid
    auto bin = gridBound.localBinsFromGlobalBin(i);
    auto center = gridBound.binCenter(bin);
    Vector2 loc{center.at(0), center.at(1)};

    // Accumulate
    Vector4 avgPos{0, 0, 0, 0};
    Vector3 avgMom{0, 0, 0};
    for (std::size_t j = 0; j < 4; j++) {
      Vector3 direction{std::sin(thetas.at(j)) * std::cos(phis.at(j)),
                        std::sin(thetas.at(j)) * std::sin(phis.at(j)),
                        std::cos(thetas.at(j))};

      avgPos += fourPositions.at(j);
      avgMom += P * direction;

      // Fill in each grid
      auto parsBound = BoundTrackParameters::create(
                           gctx, surf, fourPositions.at(j), direction, 1. / P,
                           std::nullopt, hypothesis)
                           .value();

      auto parsFree = FreeTrackParameters(fourPositions.at(j), direction,
                                          1. / P, std::nullopt, hypothesis);

      accBound.addTrack(parsBound, parsBound, loc);
      accFree.addTrack(parsFree, parsFree, loc);
    }
    avgPoss.push_back(avgPos / fourPositions.size());
    avgMoms.push_back(avgMom / fourPositions.size());
  }

  // Finalize and compare
  GridBound avgGridBound = accBound.finalizeLookup();
  GridFree avgGridFree = accFree.finalizeLookup();
  for (std::size_t i = 0; i < avgGridBound.size(); i++) {
    auto [ipBound, refBound] = avgGridBound.at(i);
    auto [ipFree, refFree] = avgGridFree.at(i);

    Vector4 avgPos = avgPoss.at(i);

    Vector3 avgMom = avgMoms.at(i);
    Vector3 avgDir = avgMom.normalized();
    double avgP = avgMom.norm();

    CHECK_CLOSE_ABS(ipBound->fourPosition(gctx), avgPos, 1e-3);
    CHECK_CLOSE_ABS(ipBound->direction(), avgDir, 1e-3);
    CHECK_CLOSE_ABS(ipBound->absoluteMomentum(), avgP, 1e-3);

    CHECK_CLOSE_ABS(ipFree->fourPosition(), avgPos, 1e-3);
    CHECK_CLOSE_ABS(ipFree->direction(), avgDir, 1e-3);
    CHECK_CLOSE_ABS(ipFree->absoluteMomentum(), avgP, 1e-3);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
