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

}  // namespace

BOOST_AUTO_TEST_SUITE(SeedingSuite)

BOOST_AUTO_TEST_CASE(DoubletStripCalibrationUsesCalibratedPositions) {
  // nominal positions off along z, as from a wrong formation direction
  SpacePointContainer container = makeContainer();
  addStripSpacePoint(container, 300_mm, {0, 0, 3_mm});
  addStripSpacePoint(container, 400_mm, {0, 0, -3_mm});

  const DoubletsForMiddleSp nominal = findDoublets(container, std::nullopt);
  BOOST_REQUIRE_EQUAL(nominal.size(), 1u);
  BOOST_CHECK_CLOSE(nominal[0].cotTheta, kDirection[2] - 0.06f, 0.1);

  const DoubletsForMiddleSp calibrated = findDoublets(container, 1.1f);
  BOOST_REQUIRE_EQUAL(calibrated.size(), 1u);
  BOOST_CHECK_SMALL(calibrated[0].cotTheta - kDirection[2], 5e-3f);
}

BOOST_AUTO_TEST_CASE(DoubletStripCalibrationRejectsMissedStrips) {
  // the chord points far off in phi and misses the strips
  SpacePointContainer container = makeContainer();
  addStripSpacePoint(container, 300_mm, {0, 0, 0});
  addStripSpacePoint(container, 400_mm, {0, 20_mm, 0});

  BOOST_CHECK_EQUAL(findDoublets(container, std::nullopt).size(), 1u);
  BOOST_CHECK_EQUAL(findDoublets(container, 1.1f).size(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
