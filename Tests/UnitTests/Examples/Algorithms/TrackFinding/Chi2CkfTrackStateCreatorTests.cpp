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
#include "ActsExamples/TrackFinding/Chi2CkfTrackStateCreatorBase.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <ranges>
#include <vector>

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {

namespace {

/// The calibrated form of a 1D measurement on `eBoundLoc0` with unit variance
struct CalibratedMeasurement {
  Vector<1> params;
  SquareMatrix<1> cov{SquareMatrix<1>::Identity()};
  std::array<std::uint8_t, 1> indices{eBoundLoc0};

  static constexpr std::size_t size() { return 1; }
  const std::array<std::uint8_t, 1>& subspaceIndices() const { return indices; }
  const Vector<1>& parameters() const { return params; }
  const SquareMatrix<1>& covariance() const { return cov; }
};

/// Creator over a flat list of loc0 values, one measurement per entry
class CustomTrackStateCreator final
    : public ActsExamples::Chi2CkfTrackStateCreatorBase<
          CustomTrackStateCreator> {
 public:
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

  CalibratedMeasurement calibrate(const State& /*state*/,
                                  const IndexSourceLink& measurement) const {
    return CalibratedMeasurement{Vector<1>{values.at(measurement.index())}};
  }
};

struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  CalibrationContext cctx;
  std::unique_ptr<const Logger> logger =
      getDefaultLogger("Chi2CkfTrackStateCreatorTests", Logging::INFO);

  std::shared_ptr<const Surface> surface =
      CurvilinearSurface(Vector3::Zero(), Vector3::UnitZ()).planeSurface();

  VectorMultiTrajectory trajectory;
  CustomTrackStateCreator creator;

  /// Predicted loc0 = 0 with unit variance, so the chi2 of a measurement at
  /// `value` is `value^2 / 2`
  CkfTypes::BoundState boundState{
      BoundTrackParameters(surface, BoundVector::Zero(),
                           BoundMatrix(BoundMatrix::Identity()),
                           ParticleHypothesis::pion()),
      BoundMatrix::Identity(), 1.5};

  auto run() {
    return creator.createTrackStates({.geoContext = &gctx,
                                      .calibrationContext = &cctx,
                                      .surface = surface.get(),
                                      .boundState = &boundState,
                                      .prevTip = kTrackIndexInvalid,
                                      .trajectory = &trajectory,
                                      .logger = logger.get()});
  }
};

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

BOOST_AUTO_TEST_SUITE(ExamplesChi2CkfTrackStateCreatorSuite)

/// A single candidate with a non-finite chi2 must not be selected
BOOST_AUTO_TEST_CASE(NonFiniteChi2IsRejected) {
  Fixture f;
  f.creator.setValues({kNaN});
  f.creator.cuts.chi2CutOffOutlier = std::numeric_limits<double>::infinity();

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result->empty());
}

/// A finite candidate must be selected next to a non-finite one, and a
/// candidate above the measurement cut must stay an outlier
BOOST_AUTO_TEST_CASE(FiniteChi2SelectedAlongsideNonFinite) {
  Fixture f;
  f.creator.cuts.chi2CutOffMeasurement = 15;
  f.creator.cuts.chi2CutOffOutlier = 25;

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
  Fixture f;
  f.creator.setValues({0.});
  f.creator.cuts.maxBranchesPerSurface = 0;

  auto result = f.run();

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    CombinatorialKalmanFilterError::MeasurementSelectionFailed);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
