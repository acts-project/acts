// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include <random>

using namespace Acts;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;

using RandomEngine = std::mt19937;
using uniform = std::uniform_real_distribution<double>;
using Line_t = CompSpacePointAuxiliaries::Line_t;
constexpr auto logLvl = Acts::Logging::Level::INFO;

namespace Acts::Test {

class FitTestSpacePoint {
 public:
  const Vector3& localPosition() const { return m_position; }
  /// @brief Wire direction
  const Vector3& sensorDirection() const { return m_sensorDir; }
  /// @brief To next sensor in the plane
  const Vector3& toNextSensor() const { return m_toNextSen; }
  /// @brief To next straw layer
  const Vector3& planeNormal() const { return m_planeNorm; }
  /// @brief Measurement's radius
  double driftRadius() const { return m_driftR.value_or(0.); }
  /// @brief Measurement's covariance
  const std::array<double, 3>& covariance() const { return m_covariance; }
  /// @brief Time of record
  double time() const { return m_time.value_or(0.); }
  /// @brief All measurements are straws
  bool isStraw() const { return m_driftR.has_value(); }
  /// @brief Check whether the space point has a time value
  bool hasTime() const { return m_time.has_value(); }
  /// @brief Check whether the space point measures the non-bending direction
  bool measuresLoc0() const { return m_measLoc0; }
  /// @brief Check whether the space point measures the bending direction
  bool measuresLoc1() const { return m_measLoc1 || isStraw(); }

 private:
  Vector3 m_position{Vector3::Zero()};
  Vector3 m_sensorDir{Vector3::UnitX()};
  Vector3 m_toNextSen{Vector3::UnitY()};
  Vector3 m_planeNorm{m_sensorDir.cross(m_toNextSen).normalized()};
  std::optional<double> m_driftR{std::nullopt};
  std::optional<double> m_time{std::nullopt};
  std::array<double, 3> m_covariance{filledArray<double, 3>(0)};
  bool m_measLoc0{false};
  bool m_measLoc1{false};
};

using Container_t = std::vector<std::shared_ptr<FitTestSpacePoint>>;

static_assert(CompositeSpacePoint<FitTestSpacePoint>);
static_assert(CompositeSpacePointContainer<Container_t>);
class SpCalibrator {
 public:
  Container_t calibrate(const Acts::CalibrationContext& /*ctx*/,
                        const Vector3& trackPos, const Vector3& trackDir,
                        const double timeOffSet,
                        const Container_t& uncalibCont) const;

  void updateSigns(const Vector3& trackPos, const Vector3& trackDir,
                   Container_t& measurements) const;
};
static_assert(
    CompositeSpacePointCalibrator<SpCalibrator, Container_t, Container_t>);

ACTS_LOCAL_LOGGER(getDefaultLogger("StrawLineFitterTest", logLvl));

BOOST_AUTO_TEST_SUITE(FastStrawLineFitTests)

BOOST_AUTO_TEST_CASE(SimpleLineFit) {
  RandomEngine engine{1503};

  using FitOpts_t =
      CompositeSpacePointLineFitter::FitOptions<Container_t, SpCalibrator>;
  using FitResult_t = CompositeSpacePointLineFitter::FitResult<Container_t>;

  CompositeSpacePointLineFitter::Config cfg{};

  // CompositeSpacePointLineFitter fitter{cfg};

  Container_t measurements{};
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
