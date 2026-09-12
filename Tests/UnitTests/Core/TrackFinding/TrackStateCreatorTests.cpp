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
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilterError.hpp"
#include "Acts/TrackFinding/TrackStateCreatorBase.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <ranges>
#include <vector>

using namespace Acts;

namespace ActsTests {

using TrackContainer =
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                   detail::ValueHolder>;

namespace {

/// A 1d measurement on `eBoundLoc0` with an identity-scaled variance.
struct TestMeasurement {
  double value{};
  double variance{1.};
  std::size_t id{};
};

/// The calibrated form of a @ref TestMeasurement, i.e. about the smallest
/// thing satisfying @ref Acts::MeasurementConcept. Its dimension is known at
/// compile time, so the creator never dispatches on it at runtime.
struct CalibratedTestMeasurement {
  Vector<1> params;
  SquareMatrix<1> cov;
  std::array<std::uint8_t, 1> indices{eBoundLoc0};

  static constexpr std::size_t size() { return 1; }
  const std::array<std::uint8_t, 1>& subspaceIndices() const { return indices; }
  const Vector<1>& parameters() const { return params; }
  const SquareMatrix<1>& covariance() const { return cov; }
};

static_assert(StaticMeasurementConcept<CalibratedTestMeasurement>,
              "The test measurement should be usable without a runtime "
              "dispatch on its dimension");

/// Minimal creator over a flat list of `TestMeasurement`.
struct TestCreator final
    : public TrackStateCreatorBase<TestCreator, TrackContainer> {
  using Base = TrackStateCreatorBase<TestCreator, TrackContainer>;
  using State = typename Base::State;

  std::vector<TestMeasurement> measurements;
  /// Only measurements with these ids pass the preselection
  std::vector<std::size_t> preselected;

  /// Cuts the selection below applies
  std::size_t numMeasurements = 1;
  double chi2Measurement = std::numeric_limits<double>::max();
  double chi2Outlier = std::numeric_limits<double>::max();
  /// If set, the selection fails with this error
  std::optional<CombinatorialKalmanFilterError> selectionError;

  mutable std::size_t nSelectMeasurements{0};
  /// Ids of the measurements which were calibrated, in order. Note that a
  /// selected measurement is calibrated twice by default, once to get its
  /// chi2 and once to fill the track state.
  mutable std::vector<std::size_t> calibrated;

  auto measurementRange(const State& /*state*/) const {
    return std::ranges::subrange(measurements.begin(), measurements.end());
  }

  bool hasMeasurementPreselection(const State& /*state*/) const {
    return !preselected.empty();
  }

  bool preselectMeasurement(const State& /*state*/,
                            const TestMeasurement& measurement) const {
    return std::ranges::find(preselected, measurement.id) != preselected.end();
  }

  SourceLink sourceLink(const State& /*state*/,
                        const TestMeasurement& measurement) const {
    return SourceLink{measurement.id};
  }

  CalibratedTestMeasurement calibrate(
      const State& /*state*/, const TestMeasurement& measurement) const {
    calibrated.push_back(measurement.id);
    return CalibratedTestMeasurement{
        Vector<1>{measurement.value},
        SquareMatrix<1>{SquareMatrix<1>::Constant(measurement.variance)},
        {static_cast<std::uint8_t>(eBoundLoc0)}};
  }

  /// The chi2 policy the shipped creators use: nothing below the measurement
  /// cut means a single outlier if it is below the outlier cut, otherwise
  /// nothing at all; measurements and outliers are never mixed.
  template <typename range_t, typename candidates_t>
  Result<void> selectMeasurements(const State& state, range_t&& range,
                                  candidates_t& candidates) const {
    ++nSelectMeasurements;

    if (selectionError.has_value()) {
      return *selectionError;
    }
    if (numMeasurements == 0) {
      return CombinatorialKalmanFilterError::MeasurementSelectionFailed;
    }

    if (hasMeasurementPreselection(state)) {
      for (const TestMeasurement& measurement : range) {
        if (preselectMeasurement(state, measurement)) {
          candidates.push_back({measurement, 0., false});
        }
      }
    }
    if (candidates.empty()) {
      for (const TestMeasurement& measurement : range) {
        candidates.push_back({measurement, 0., false});
      }
    }

    for (auto& candidate : candidates) {
      candidate.chi2 = calculatePredictedChi2(
          state, calibrate(state, candidate.measurement));
    }

    const auto best = std::ranges::min_element(
        candidates, {}, [](const auto& candidate) { return candidate.chi2; });
    const bool isOutlier = best->chi2 >= chi2Measurement;
    if (isOutlier) {
      const auto selected = *best;
      candidates.clear();
      if (selected.chi2 < chi2Outlier) {
        candidates.push_back({selected.measurement, selected.chi2, true});
      }
      return Result<void>::success();
    }

    candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                    [this](const auto& candidate) {
                                      return candidate.chi2 >= chi2Measurement;
                                    }),
                     candidates.end());
    std::ranges::sort(candidates, {},
                      [](const auto& candidate) { return candidate.chi2; });
    candidates.erase(
        candidates.begin() + std::min(numMeasurements, candidates.size()),
        candidates.end());
    return Result<void>::success();
  }
};

struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  CalibrationContext cctx;
  std::unique_ptr<const Logger> logger =
      getDefaultLogger("TrackStateCreatorTests", Logging::INFO);

  std::shared_ptr<const Surface> surface =
      CurvilinearSurface(Vector3::Zero(), Vector3::UnitZ()).planeSurface();

  VectorMultiTrajectory trajectory;
  TestCreator creator;

  /// Predicted parameters at loc0 = 0 with unit variance, so the chi2 of a
  /// measurement at `value` with variance `v` is `value^2 / (1 + v)`.
  CkfTypes::BoundState boundState{
      BoundTrackParameters(surface, BoundVector::Zero(),
                           BoundMatrix(BoundMatrix::Identity()),
                           ParticleHypothesis::pion()),
      BoundMatrix::Identity(), 1.5};

  void setCuts(std::size_t numMeasurements, double chi2Measurement,
               double chi2Outlier) {
    creator.numMeasurements = numMeasurements;
    creator.chi2Measurement = chi2Measurement;
    creator.chi2Outlier = chi2Outlier;
  }

  Result<CkfTypes::BranchVector<TrackIndexType>> run() {
    return creator.createTrackStates(gctx, cctx, *surface, boundState,
                                     kTrackIndexInvalid, trajectory, *logger);
  }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFindingSuite)

/// A surface without measurements must not even resolve the cuts. Otherwise
/// surfaces which are not covered by the selector configuration would turn
/// into hard errors instead of holes.
BOOST_AUTO_TEST_CASE(NoMeasurementsCreatesNothing) {
  Fixture f;
  f.setCuts(1u, 10., 100.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result->empty());
  BOOST_CHECK_EQUAL(f.trajectory.size(), 0u);
  BOOST_CHECK_EQUAL(f.creator.nSelectMeasurements, 0u);
}

/// A selection which cannot make up its mind must abort the track finding
/// rather than silently turn the surface into a hole.
BOOST_AUTO_TEST_CASE(SelectionErrorIsPropagated) {
  Fixture f;
  f.creator.measurements = {{0., 1., 0u}};
  f.creator.selectionError =
      CombinatorialKalmanFilterError::MeasurementSelectionFailed;

  auto result = f.run();

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    CombinatorialKalmanFilterError::MeasurementSelectionFailed);
}

BOOST_AUTO_TEST_CASE(ZeroNumMeasurementsFails) {
  Fixture f;
  f.creator.measurements = {{0., 1., 0u}};
  f.setCuts(0u, 10., 100.);

  auto result = f.run();

  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    CombinatorialKalmanFilterError::MeasurementSelectionFailed);
}

/// Nothing passes the measurement cut but the best candidate is within the
/// outlier cut: exactly one outlier state, never a mix of the two.
BOOST_AUTO_TEST_CASE(BestCandidateBecomesSingleOutlier) {
  Fixture f;
  // chi2 = value^2 / 2 -> 8, 18, 32
  f.creator.measurements = {{4., 1., 0u}, {6., 1., 1u}, {8., 1., 2u}};
  f.setCuts(3u, 5., 20.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);

  const auto trackState = f.trajectory.getTrackState(result->front());
  BOOST_CHECK(trackState.typeFlags().isOutlier());
  BOOST_CHECK(!trackState.typeFlags().isHole());
  BOOST_CHECK_CLOSE(trackState.chi2(), 8., 1e-9);
  BOOST_CHECK_EQUAL(trackState.getUncalibratedSourceLink().get<std::size_t>(),
                    0u);
  BOOST_CHECK_EQUAL(f.trajectory.size(), 1u);
}

/// Nothing passes either cut: no track state at all, the caller creates the
/// hole. Creating one here would make the CKF produce two states.
BOOST_AUTO_TEST_CASE(NothingCompatibleCreatesNoTrackState) {
  Fixture f;
  f.creator.measurements = {{4., 1., 0u}, {6., 1., 1u}};
  f.setCuts(3u, 5., 6.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(result->empty());
  BOOST_CHECK_EQUAL(f.trajectory.size(), 0u);
}

BOOST_AUTO_TEST_CASE(SinglePassingCandidateIsPicked) {
  Fixture f;
  // only the last one passes, so it must be found even though it is not first
  f.creator.measurements = {{4., 1., 0u}, {6., 1., 1u}, {1., 1., 2u}};
  f.setCuts(3u, 5., 100.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);

  const auto trackState = f.trajectory.getTrackState(result->front());
  BOOST_CHECK(!trackState.typeFlags().isOutlier());
  BOOST_CHECK(trackState.typeFlags().isMeasurement());
  BOOST_CHECK_EQUAL(trackState.getUncalibratedSourceLink().get<std::size_t>(),
                    2u);
}

/// More candidates pass than are allowed: the best ones survive and they come
/// out ordered by ascending chi2, which is what the CKF branching relies on.
BOOST_AUTO_TEST_CASE(SelectedAreTruncatedAndOrderedByChi2) {
  Fixture f;
  // chi2 = value^2 / 2 -> 4.5, 0.5, 2, 8
  f.creator.measurements = {
      {3., 1., 0u}, {1., 1., 1u}, {2., 1., 2u}, {4., 1., 3u}};
  f.setCuts(3u, 10., 100.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 3u);

  std::vector<std::size_t> ids;
  double previousChi2 = -1.;
  for (const TrackIndexType index : *result) {
    const auto trackState = f.trajectory.getTrackState(index);
    BOOST_CHECK_GT(trackState.chi2(), previousChi2);
    previousChi2 = trackState.chi2();
    ids.push_back(trackState.getUncalibratedSourceLink().get<std::size_t>());
  }
  const std::vector<std::size_t> expected = {1u, 2u, 0u};
  BOOST_CHECK_EQUAL_COLLECTIONS(ids.begin(), ids.end(), expected.begin(),
                                expected.end());
}

/// The first track state on a surface owns the predicted parameters and the
/// jacobian, all further ones share them.
BOOST_AUTO_TEST_CASE(TrackStatesShareTheirPrediction) {
  Fixture f;
  f.creator.measurements = {{1., 1., 0u}, {2., 1., 1u}};
  f.setCuts(2u, 10., 100.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 2u);

  const auto first = f.trajectory.getTrackState(result->at(0));
  const auto second = f.trajectory.getTrackState(result->at(1));

  BOOST_CHECK_EQUAL(first.predicted().data(), second.predicted().data());
  BOOST_CHECK_EQUAL(first.jacobian().data(), second.jacobian().data());
  BOOST_CHECK_NE(first.effectiveCalibrated().data(),
                 second.effectiveCalibrated().data());
}

/// Everything the combinatorial Kalman filter reads off a new track state.
BOOST_AUTO_TEST_CASE(TrackStateIsFullyPopulated) {
  Fixture f;
  f.creator.measurements = {{1., 2., 7u}};
  f.setCuts(1u, 10., 100.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);

  const auto trackState = f.trajectory.getTrackState(result->front());

  BOOST_CHECK(trackState.typeFlags().hasParameters());
  BOOST_CHECK(trackState.typeFlags().isMeasurement());
  BOOST_CHECK(!trackState.typeFlags().isOutlier());
  BOOST_CHECK(!trackState.typeFlags().isHole());

  BOOST_CHECK_EQUAL(trackState.pathLength(), 1.5);
  // chi2 = 1^2 / (1 + 2), stored as float
  BOOST_CHECK_EQUAL(trackState.chi2(), static_cast<float>(1. / 3.));
  BOOST_CHECK_EQUAL(&trackState.referenceSurface(), f.surface.get());
  BOOST_CHECK_EQUAL(trackState.predicted(), BoundVector::Zero());
  BOOST_CHECK_EQUAL(trackState.predictedCovariance(), BoundMatrix::Identity());
  BOOST_CHECK_EQUAL(trackState.jacobian(), BoundMatrix::Identity());

  // the filtered parameters must not be allocated, the caller allocates them
  // right before the Kalman update writes them
  BOOST_CHECK(!trackState.hasFiltered());
  BOOST_CHECK(!trackState.hasSmoothed());

  BOOST_REQUIRE_EQUAL(trackState.calibratedSize(), 1u);
  BOOST_CHECK_EQUAL(trackState.effectiveCalibrated()[0], 1.);
  BOOST_CHECK_EQUAL(trackState.effectiveCalibratedCovariance()(0, 0), 2.);
  // the projector is serialized including its unused slots, which the
  // calibrators leave at zero
  const BoundSubspaceIndices expectedIndices = {eBoundLoc0, 0, 0, 0, 0, 0};
  BOOST_CHECK(trackState.projectorSubspaceIndices() == expectedIndices);

  BOOST_CHECK_EQUAL(trackState.getUncalibratedSourceLink().get<std::size_t>(),
                    7u);
}

BOOST_AUTO_TEST_CASE(PreselectionSkipsMeasurements) {
  Fixture f;
  f.creator.measurements = {{1., 1., 0u}, {2., 1., 1u}, {3., 1., 2u}};
  f.creator.preselected = {2u};
  f.setCuts(3u, 100., 200.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  BOOST_CHECK_EQUAL(f.trajectory.getTrackState(result->front())
                        .getUncalibratedSourceLink()
                        .get<std::size_t>(),
                    2u);
  // the rejected measurements are never calibrated at all
  BOOST_CHECK(std::ranges::all_of(f.creator.calibrated,
                                  [](std::size_t id) { return id == 2u; }));
}

/// If the preselection rejects everything it is ignored, matching the "stay on
/// seed" behaviour of only restricting surfaces the seed actually touches.
BOOST_AUTO_TEST_CASE(PreselectionRejectingAllFallsBackToAll) {
  Fixture f;
  f.creator.measurements = {{1., 1., 0u}, {2., 1., 1u}};
  f.creator.preselected = {42u};
  f.setCuts(3u, 100., 200.);

  auto result = f.run();

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->size(), 2u);
  // both measurements were considered despite matching no preselection
  BOOST_CHECK(std::ranges::find(f.creator.calibrated, 0u) !=
              f.creator.calibrated.end());
  BOOST_CHECK(std::ranges::find(f.creator.calibrated, 1u) !=
              f.creator.calibrated.end());
}

/// The chi2 the creator selects on has to be the same quantity that is
/// computed from a filled track state, otherwise the selection silently
/// drifts away from the rest of the track EDM.
BOOST_AUTO_TEST_CASE(Chi2MatchesTrackStateBasedCalculation) {
  Fixture f;
  // deliberately a value whose chi2 is not exactly representable
  f.creator.measurements = {{1.3, 2.7, 0u}};
  f.setCuts(1u, 100., 200.);

  auto result = f.run();
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->size(), 1u);
  const auto trackState = f.trajectory.getTrackState(result->front());

  // build the equivalent track state by hand and run the independent track
  // state based implementation on it
  VectorMultiTrajectory reference;
  auto referenceState = reference.makeTrackState(
      TrackStatePropMask::Predicted | TrackStatePropMask::Calibrated);
  referenceState.predicted() = BoundVector::Zero();
  referenceState.predictedCovariance() = BoundMatrix::Identity();
  referenceState.allocateCalibrated(
      Vector<1>{1.3}, SquareMatrix<1>{SquareMatrix<1>::Constant(2.7)});
  referenceState.setProjectorSubspaceIndices(
      std::array<std::uint8_t, 1>{eBoundLoc0});

  const double expected = Acts::calculatePredictedChi2(referenceState);

  // `chi2` is stored as float, so this is the exact round trip
  BOOST_CHECK_EQUAL(trackState.chi2(), static_cast<float>(expected));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
