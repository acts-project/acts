// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <random>

using namespace Acts;
using namespace Acts::Experimental::detail;
using namespace Acts::UnitLiterals;
using RandomEngine = std::mt19937;

namespace Acts::Test {

class StrawTestPoint {
 public:
  using ResidualIdx = FastStrawLineFitter::ResidualIdx;
  StrawTestPoint(const Vector3& pos, const double driftR,
                 const double driftRUncert)
      : m_pos{pos}, m_driftR{driftR} {
    m_cov[toUnderlying(ResidualIdx::bending)] = Acts::pow(driftRUncert, 2);
  }
  /// @brief Straw tube's direction
  const Vector3& localPosition() const { return m_pos; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_wireDir; }
  /// @brief Dummy return not used in test
  const Vector3& toNextSensor() const { return m_wireDir; }
  /// @brief Dummy return not used in test
  const Vector3& planeNormal() const { return m_wireDir; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR; }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_cov; }
  /// @brief Dummy return not used in test
  double time() const { return 0.; }

  bool isStraw() const { return true; }
  bool hasTime() const { return false; }
  bool measuresLoc0() const { return false; }
  bool measuresLoc1() const { return false; }

 private:
  Vector3 m_pos{Vector3::Zero()};
  Vector3 m_wireDir{Vector3::UnitX()};
  double m_driftR{0.};
  std::array<double, 3> m_cov{Acts::filledArray<double, 3>(0.)};
};
static_assert(Acts::Experimental::CompositeSpacePoint<StrawTestPoint>);

using TestStrawCont_t = std::vector<std::unique_ptr<StrawTestPoint>>;

using Line_t = CompSpacePointAuxiliaries::Line_t;

Line_t generateLine(RandomEngine& engine) {
  using ParIndex = Line_t::ParIndex;
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(ParIndex::x0)] = 0.;
  linePars[toUnderlying(ParIndex::phi)] = 90._degree;
  linePars[toUnderlying(ParIndex::y0)] = (engine() % 10000 - 5000.) / 10.;
  linePars[toUnderlying(ParIndex::theta)] = (engine() % 1800) * 0.1_degree;
  std::cout << "Generated parameters theta: "
            << (linePars[toUnderlying(ParIndex::theta)] / 1._degree)
            << ", y0: " << linePars[toUnderlying(ParIndex::y0)] << std::endl;

  Line_t line{};
  line.updateParameters(linePars);
  return line;
}

double calcDriftUncert(const double driftR) {
  return 0.15_mm * Acts::pow(1._mm + Acts::abs(driftR), -1);
}

TestStrawCont_t generateStrawCircles(const Line_t& trajLine,
                                     RandomEngine& engine) {
  const Vector3 posStaggering{0., std::cos(60._degree), std::sin(60._degree)};
  const Vector3 negStaggering{0., -std::cos(60._degree), std::sin(60._degree)};
  /// Number of tube layers per multilayer
  constexpr std::uint32_t nLayersPerMl = 4;
  /// Number of overall tubelayers
  constexpr std::uint32_t nTubeLayers = nLayersPerMl * 2;
  constexpr double tubeRadius = 15._mm;
  constexpr double tubeLayerDist = 60._cm;

  std::array<Vector3, nTubeLayers> tubePositions{
      filledArray<Vector3, nTubeLayers>(Vector3{0., tubeRadius, tubeRadius})};

  for (std::uint32_t l = 1; l < nTubeLayers; ++l) {
    const Vector3& layStag{l % 2 == 1 ? posStaggering : negStaggering};
    tubePositions[l] = tubePositions[l - 1] + 2. * tubeRadius * layStag;

    if (l == nLayersPerMl) {
      tubePositions[l] += tubeLayerDist * Vector3::UnitZ();
    }
  }
  /// Print the staggering
  std::cout << "##############################################" << std::endl;
  for (std::uint32_t l = 0; l < nTubeLayers; ++l) {
    std::cout << " *** " << (l + 1) << " - " << toString(tubePositions[l])
              << std::endl;
  }
  std::cout << "##############################################" << std::endl;
  TestStrawCont_t circles{};
  for (const auto& stag : tubePositions) {
    auto planeExtpLow = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() - tubeRadius);
    auto planeExtpHigh = Acts::PlanarHelper::intersectPlane(
        trajLine.position(), trajLine.direction(), Vector3::UnitZ(),
        stag.z() + tubeRadius);

    std::cout << "extrapolated to plane " << toString(planeExtpLow.position())
              << " " << toString(planeExtpHigh.position()) << std::endl;
    const auto dToFirstLow = static_cast<std::int32_t>(std::ceil(
        (planeExtpLow.position().y() - stag.y()) / (2. * tubeRadius)));
    const auto dToFirstHigh = static_cast<std::int32_t>(std::ceil(
        (planeExtpHigh.position().y() - stag.y()) / (2. * tubeRadius)));

    const std::int32_t dT = dToFirstHigh > dToFirstLow ? 1 : -1;
    for (std::int32_t tN = dToFirstLow; tN != dToFirstHigh; tN += dT) {
      const Vector3 tube = stag + 2. * tN * tubeRadius * Vector3::UnitY();
      const double rad = Acts::detail::LineHelper::signedDistance(
          tube, Vector3::UnitX(), trajLine.position(), trajLine.direction());
      if (std::abs(rad) > tubeRadius) {
        continue;
      }
      std::normal_distribution<> dist{rad, calcDriftUncert(rad)};
      const double smearedR = std::abs(dist(engine));
      if (smearedR > tubeRadius) {
        continue;
      }
      circles.emplace_back(std::make_unique<StrawTestPoint>(
          tube, smearedR, calcDriftUncert(smearedR)));
    }
  }
  std::cout << " Track hit in total " << circles.size() << " tubes "
            << std::endl;
  return circles;
}

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(StrawDriftTimeCase) {
  constexpr std::uint32_t nTrials = 1000;
  RandomEngine engine{1419};
  for (std::uint32_t n = 0; n < nTrials; ++n) {
    auto track = generateLine(engine);
    auto strawPoints = generateStrawCircles(track, engine);
  }
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
