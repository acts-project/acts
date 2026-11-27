// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

namespace {

// Track direction from the origin, slightly off the stereo strips
constexpr std::array<float, 3> kDirection{1.f, 4e-4f, 0.5f};
constexpr float kStripHalfLength = 12_mm;
constexpr float kStereoAngle = 0.02f;
constexpr float kModuleGap = 2_mm;

// Adds a strip space point of a barrel stereo module at radius `r`, crossed
// by `kDirection`, at its nominal position shifted by `offset`.
void addStripSpacePoint(SpacePointContainer& container, float r,
                        const std::array<float, 3>& offset) {
  const std::array<float, 3> innerHalf{
      0, -kStripHalfLength * std::sin(kStereoAngle),
      kStripHalfLength * std::cos(kStereoAngle)};
  const std::array<float, 3> outerHalf{
      0, kStripHalfLength * std::sin(kStereoAngle),
      kStripHalfLength * std::cos(kStereoAngle)};

  // strip center such that the track hit at `x` lies on the strip
  const auto center = [&](float x, const std::array<float, 3>& half) {
    const float s = kDirection[1] * x / half[1];
    return std::array<float, 3>{x, kDirection[1] * x - s * half[1],
                                kDirection[2] * x - s * half[2]};
  };
  const std::array<float, 3> innerCenter = center(r, innerHalf);
  const std::array<float, 3> outerCenter = center(r + kModuleGap, outerHalf);

  auto sp = container.createSpacePoint();
  const float x = r + offset[0];
  const float y = kDirection[1] * r + offset[1];
  sp.xy() = {x, y};
  sp.zr() = {kDirection[2] * r + offset[2], std::hypot(x, y)};
  sp.varianceZ() = 1_mm;
  sp.varianceR() = 0;
  sp.outerStripCalibrationDetails() = {
      outerCenter,
      {outerCenter[0] - innerCenter[0], outerCenter[1] - innerCenter[1],
       outerCenter[2] - innerCenter[2]},
      outerHalf,
      innerHalf};
}

// Finds the top doublets of the first space point among the others
DoubletsForMiddleSp findDoublets(const SpacePointContainer& container,
                                 std::optional<float> tolerance) {
  DoubletSeedFinder::Config config;
  config.candidateDirection = Direction::Forward();
  config.deltaRMin = 10_mm;
  config.deltaRMax = 200_mm;
  config.collisionRegionMin = -50_mm;
  config.collisionRegionMax = 50_mm;
  config.stripCalibrationTolerance = tolerance;
  const auto finder =
      DoubletSeedFinder::create(DoubletSeedFinder::DerivedConfig(config, 2_T));

  const ConstSpacePointProxy middle = container[0];
  DoubletsForMiddleSp doublets;
  finder->createDoublets(middle, DoubletSeedFinder::computeMiddleSpInfo(middle),
                         container.range({1, container.size()}), doublets);
  return doublets;
}

SpacePointContainer makeContainer() {
  return SpacePointContainer(
      SpacePointColumns::PackedXY | SpacePointColumns::PackedZR |
      SpacePointColumns::VarianceZ | SpacePointColumns::VarianceR |
      SpacePointColumns::StripCalibrationDetails);
}

// Adds a pixel space point at (x, y), with z on a line of cotTheta 0.5 from
// the origin
void addPixelSpacePoint(SpacePointContainer& container, float x, float y) {
  auto sp = container.createSpacePoint();
  const float r = std::hypot(x, y);
  sp.xy() = {x, y};
  sp.zr() = {0.5f * r, r};
  sp.varianceZ() = 0.01_mm;
  sp.varianceR() = 0.01_mm;
}

// Finds the doublets of the first space point among those from `first` to
// before `last`, with the interaction point cut
DoubletsForMiddleSp findPixelDoublets(const SpacePointContainer& container,
                                      Direction direction,
                                      SpacePointIndex first,
                                      SpacePointIndex last) {
  DoubletSeedFinder::Config config;
  config.candidateDirection = direction;
  config.interactionPointCut = true;
  const auto finder =
      DoubletSeedFinder::create(DoubletSeedFinder::DerivedConfig(config, 2_T));

  const ConstSpacePointProxy middle = container[0];
  DoubletsForMiddleSp doublets;
  finder->createDoublets(middle, DoubletSeedFinder::computeMiddleSpInfo(middle),
                         container.range({first, last}), doublets);
  return doublets;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(DoubletStripCalibrationUsesCalibratedPositions) {
  // nominal positions off along z, as from a wrong formation direction
  SpacePointContainer container = makeContainer();
  addStripSpacePoint(container, 300_mm, {0, 0, 3_mm});
  addStripSpacePoint(container, 400_mm, {0, 0, -3_mm});

  const DoubletsForMiddleSp nominal = findDoublets(container, std::nullopt);
  BOOST_REQUIRE_EQUAL(nominal.size(), 1u);
  BOOST_CHECK_CLOSE(nominal[0].cotTheta(), kDirection[2] - 0.06f, 0.1);

  const DoubletsForMiddleSp calibrated = findDoublets(container, 1.1f);
  BOOST_REQUIRE_EQUAL(calibrated.size(), 1u);
  BOOST_CHECK_SMALL(calibrated[0].cotTheta() - kDirection[2], 5e-3f);
}

BOOST_AUTO_TEST_CASE(DoubletStripCalibrationRejectsMissedStrips) {
  // the chord points far off in phi and misses the strips
  SpacePointContainer container = makeContainer();
  addStripSpacePoint(container, 300_mm, {0, 0, 0});
  addStripSpacePoint(container, 400_mm, {0, 20_mm, 0});

  BOOST_CHECK_EQUAL(findDoublets(container, std::nullopt).size(), 1u);
  BOOST_CHECK_EQUAL(findDoublets(container, 1.1f).size(), 0u);
}

BOOST_AUTO_TEST_CASE(DoubletInteractionPointCutsOverBlocks) {
  // more candidates than two blocks of the finder, every third far from the
  // interaction point and too curved, the others on lines through it: one
  // call over all of them gives the same doublets in the same order as one
  // call per candidate
  constexpr std::size_t nCandidates = 600;
  for (const bool isBottom : {false, true}) {
    const float rM = isBottom ? 400_mm : 100_mm;
    SpacePointContainer container(
        SpacePointColumns::PackedXY | SpacePointColumns::PackedZR |
        SpacePointColumns::VarianceZ | SpacePointColumns::VarianceR);
    addPixelSpacePoint(container, rM, 0);
    for (std::size_t i = 0; i < nCandidates; ++i) {
      const float r = isBottom ? 60_mm + 0.46_mm * i : 150_mm + 0.3_mm * i;
      if (i % 3 == 0) {
        const float slope = isBottom ? -0.3f : 0.6f;
        addPixelSpacePoint(container, r, slope * (r - rM));
      } else {
        const float phi = 1e-4f * static_cast<float>(i % 5);
        addPixelSpacePoint(container, r * std::cos(phi), r * std::sin(phi));
      }
    }
    const Direction direction =
        isBottom ? Direction::Backward() : Direction::Forward();

    const auto last = static_cast<SpacePointIndex>(container.size());
    const DoubletsForMiddleSp together =
        findPixelDoublets(container, direction, 1, last);
    DoubletsForMiddleSp alone;
    for (SpacePointIndex i = 1; i < last; ++i) {
      const DoubletsForMiddleSp one =
          findPixelDoublets(container, direction, i, i + 1);
      for (const auto doublet : one) {
        alone.emplace_back(doublet.spacePointIndex(), doublet.cotTheta(),
                           doublet.iDeltaR(), doublet.er(), doublet.u(),
                           doublet.v(), doublet.x(), doublet.y());
      }
    }

    BOOST_CHECK_GT(together.size(), 0u);
    BOOST_CHECK_LT(together.size(), 2 * nCandidates / 3 + 1);
    BOOST_REQUIRE_EQUAL(together.size(), alone.size());
    // exact: the same expressions on the same inputs, so the same bits
    for (DoubletsForMiddleSp::Index k = 0; k < together.size(); ++k) {
      BOOST_CHECK_EQUAL(together[k].spacePointIndex(),
                        alone[k].spacePointIndex());
      BOOST_CHECK_EQUAL(together[k].cotTheta(), alone[k].cotTheta());
      BOOST_CHECK_EQUAL(together[k].iDeltaR(), alone[k].iDeltaR());
      BOOST_CHECK_EQUAL(together[k].er(), alone[k].er());
      BOOST_CHECK_EQUAL(together[k].u(), alone[k].u());
      BOOST_CHECK_EQUAL(together[k].v(), alone[k].v());
      BOOST_CHECK_EQUAL(together[k].x(), alone[k].x());
      BOOST_CHECK_EQUAL(together[k].y(), alone[k].y());
    }
  }
}

BOOST_AUTO_TEST_CASE(DoubletsOfSeveralGroups) {
  // groups of candidates, an empty one and some longer than a block among
  // them: one call over all of them gives the same doublets in the same order
  // as one call per group
  constexpr std::size_t nCandidates = 600;
  const std::vector<std::array<SpacePointIndex, 2>> groupBounds = {
      {1, 41}, {41, 41}, {41, 341}, {341, 342}, {342, 601}};
  for (const bool isBottom : {false, true}) {
    const float rM = isBottom ? 400_mm : 100_mm;
    SpacePointContainer container(
        SpacePointColumns::PackedXY | SpacePointColumns::PackedZR |
        SpacePointColumns::VarianceZ | SpacePointColumns::VarianceR);
    addPixelSpacePoint(container, rM, 0);
    for (std::size_t i = 0; i < nCandidates; ++i) {
      const float r = isBottom ? 60_mm + 0.46_mm * i : 150_mm + 0.3_mm * i;
      const float slope = (i % 3 == 0) ? (isBottom ? -0.3f : 0.6f) : 0.f;
      addPixelSpacePoint(container, r, slope * (r - rM));
    }
    const Direction direction =
        isBottom ? Direction::Backward() : Direction::Forward();

    DoubletSeedFinder::Config config;
    config.candidateDirection = direction;
    config.interactionPointCut = true;
    const auto finder = DoubletSeedFinder::create(
        DoubletSeedFinder::DerivedConfig(config, 2_T));
    const SpacePointContainer& spacePoints = container;
    const ConstSpacePointProxy middle = spacePoints[0];
    std::vector<SpacePointContainer::ConstRange> groups;
    for (const auto& [first, last] : groupBounds) {
      groups.push_back(spacePoints.range({first, last}));
    }
    DoubletsForMiddleSp together;
    finder->createDoublets(
        middle, DoubletSeedFinder::computeMiddleSpInfo(middle),
        std::span<const SpacePointContainer::ConstRange>(groups), together);

    DoubletsForMiddleSp alone;
    for (const auto& [first, last] : groupBounds) {
      const DoubletsForMiddleSp one =
          findPixelDoublets(container, direction, first, last);
      for (const auto doublet : one) {
        alone.emplace_back(doublet.spacePointIndex(), doublet.cotTheta(),
                           doublet.iDeltaR(), doublet.er(), doublet.u(),
                           doublet.v(), doublet.x(), doublet.y());
      }
    }

    BOOST_CHECK_GT(together.size(), 0u);
    BOOST_REQUIRE_EQUAL(together.size(), alone.size());
    for (DoubletsForMiddleSp::Index k = 0; k < together.size(); ++k) {
      BOOST_CHECK_EQUAL(together[k].spacePointIndex(),
                        alone[k].spacePointIndex());
      BOOST_CHECK_EQUAL(together[k].cotTheta(), alone[k].cotTheta());
      BOOST_CHECK_EQUAL(together[k].iDeltaR(), alone[k].iDeltaR());
      BOOST_CHECK_EQUAL(together[k].er(), alone[k].er());
      BOOST_CHECK_EQUAL(together[k].u(), alone[k].u());
      BOOST_CHECK_EQUAL(together[k].v(), alone[k].v());
      BOOST_CHECK_EQUAL(together[k].x(), alone[k].x());
      BOOST_CHECK_EQUAL(together[k].y(), alone[k].y());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
