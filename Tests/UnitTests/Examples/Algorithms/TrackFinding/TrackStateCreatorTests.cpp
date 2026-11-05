// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/TrackFinding/TrackStateCreator.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {

namespace {

/// The calibrated form of a measurement of dimension `dim` with unit variance.
/// Only the first coordinate differs from the prediction.
template <std::size_t dim>
struct CalibratedMeasurement {
  Vector<dim> params;
  SquareMatrix<dim> cov{SquareMatrix<dim>::Identity()};
  std::array<std::uint8_t, dim> indices = [] {
    std::array<std::uint8_t, dim> result{};
    for (std::size_t i = 0; i < dim; ++i) {
      result[i] = static_cast<std::uint8_t>(i);
    }
    return result;
  }();

  static constexpr std::size_t size() { return dim; }
  const std::array<std::uint8_t, dim>& subspaceIndices() const {
    return indices;
  }
  const Vector<dim>& parameters() const { return params; }
  const SquareMatrix<dim>& covariance() const { return cov; }
};

/// Creator over a flat list of loc0 values, one measurement of dimension `dim`
/// per entry
template <std::size_t dim>
class TestCreator final
    : public ActsExamples::TrackStateCreatorBase<TestCreator<dim>> {
 public:
  using State =
      typename ActsExamples::TrackStateCreatorBase<TestCreator<dim>>::State;

  std::vector<double> values;
  std::vector<IndexSourceLink> sourceLinks;

  void setValues(std::vector<double> v) {
    values = std::move(v);
    sourceLinks.clear();
    for (Index i = 0; i < values.size(); ++i) {
      sourceLinks.emplace_back(GeometryIdentifier(), i);
    }
  }

  auto measurementRange(const State& /*state*/) const {
    return std::ranges::subrange(sourceLinks.begin(), sourceLinks.end());
  }

  CalibratedMeasurement<dim> calibrate(
      const State& /*state*/, const IndexSourceLink& measurement) const {
    CalibratedMeasurement<dim> result;
    result.params = Vector<dim>::Zero();
    result.params[0] = values.at(measurement.index());
    return result;
  }
};

template <std::size_t dim = 1>
struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  CalibrationContext cctx;
  std::unique_ptr<const Logger> logger =
      getDefaultLogger("TrackStateCreatorTests", Logging::INFO);

  std::shared_ptr<const Surface> surface =
      CurvilinearSurface(Vector3::Zero(), Vector3::UnitZ()).planeSurface();

  VectorMultiTrajectory trajectory;
  TestCreator<dim> creator;
  TrackStateSelectionCuts cuts;

  /// Predicted parameters at zero with unit variance, so the chi2 of a
  /// measurement at `value` is `value^2 / 2`
  CkfTypes::BoundState boundState{
      BoundTrackParameters(surface, BoundVector::Zero(),
                           BoundMatrix(BoundMatrix::Identity()),
                           ParticleHypothesis::pion()),
      BoundMatrix::Identity(), 1.5};

  auto run() {
    creator.setCuts(cuts);
    return creator.createTrackStates(gctx, cctx, *surface, boundState,
                                     kTrackIndexInvalid, trajectory, *logger);
  }
};

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

BOOST_AUTO_TEST_SUITE(ExamplesTrackStateCreatorSuite)

/// A single candidate with a non-finite chi2 must not be selected
BOOST_AUTO_TEST_CASE(NonFiniteChi2IsRejected) {
  Fixture<> f;
  f.creator.setValues({kNaN});
  f.cuts.chi2CutOffOutlier = std::numeric_limits<double>::infinity();

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result->empty());
}

/// A finite candidate must be selected next to a non-finite one, and a
/// candidate above the measurement cut must stay an outlier
BOOST_AUTO_TEST_CASE(FiniteChi2SelectedAlongsideNonFinite) {
  Fixture<> f;
  f.cuts.chi2CutOffMeasurement = 15;
  f.cuts.chi2CutOffOutlier = 25;

  // chi2 = 2
  f.creator.setValues({kNaN, 2.});
  auto result = f.run();
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  auto trackState = f.trajectory.getTrackState(result->front());
  BOOST_CHECK(!trackState.typeFlags().isOutlier());
  BOOST_CHECK_CLOSE(trackState.chi2(), 2., 1e-9);

  // chi2 = 18
  f.creator.setValues({kNaN, 6.});
  result = f.run();
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  trackState = f.trajectory.getTrackState(result->front());
  BOOST_CHECK(trackState.typeFlags().isOutlier());
  BOOST_CHECK_CLOSE(trackState.chi2(), 18., 1e-9);
}

/// No branch allowed on a surface is a configuration error
BOOST_AUTO_TEST_CASE(ZeroBranchesFails) {
  Fixture<> f;
  f.creator.setValues({0.});
  f.cuts.maxBranchesPerSurface = 0;

  auto result = f.run();

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    CombinatorialKalmanFilterError::MeasurementSelectionFailed);
}

/// Selection of one measurement with the given chi2: (selected, isOutlier)
template <std::size_t dim>
std::pair<bool, bool> selectWithChi2(double chi2, bool perDimension) {
  Fixture<dim> f;
  f.cuts.chi2CutOffMeasurement = 15;
  f.cuts.chi2CutOffOutlier = 25;
  f.cuts.chi2CutPerDimension = perDimension;
  // chi2 = value^2 / 2
  f.creator.setValues({std::sqrt(2 * chi2)});

  auto result = f.run();
  BOOST_REQUIRE(result.ok());
  if (result->empty()) {
    return {false, false};
  }
  const auto trackState = f.trajectory.getTrackState(result->front());
  return {true, trackState.typeFlags().isOutlier()};
}

/// Without the flag every measurement dimension uses the same cuts
BOOST_AUTO_TEST_CASE(Chi2CutSameForAllDimensions) {
  BOOST_CHECK(selectWithChi2<1>(14.9, false) == std::pair(true, false));
  BOOST_CHECK(selectWithChi2<3>(15.1, false) == std::pair(true, true));
  BOOST_CHECK(selectWithChi2<3>(25.1, false) == std::pair(false, false));
}

/// With the flag the cuts are the 1D cuts for 1D measurements and the chi2
/// values with the same p-value for higher dimensions: 18.28 / 28.74 for 2D
/// and 20.96 / 31.81 for 3D
BOOST_AUTO_TEST_CASE(Chi2CutPerDimension) {
  BOOST_CHECK(selectWithChi2<1>(14.9, true) == std::pair(true, false));
  BOOST_CHECK(selectWithChi2<1>(15.1, true) == std::pair(true, true));

  BOOST_CHECK(selectWithChi2<2>(18.2, true) == std::pair(true, false));
  BOOST_CHECK(selectWithChi2<2>(18.4, true) == std::pair(true, true));

  BOOST_CHECK(selectWithChi2<3>(20.9, true) == std::pair(true, false));
  BOOST_CHECK(selectWithChi2<3>(21.0, true) == std::pair(true, true));
  BOOST_CHECK(selectWithChi2<3>(31.7, true) == std::pair(true, true));
  BOOST_CHECK(selectWithChi2<3>(31.9, true) == std::pair(false, false));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
