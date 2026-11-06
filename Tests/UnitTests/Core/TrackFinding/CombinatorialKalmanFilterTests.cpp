// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CkfTrackStateCreatorBase.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsTests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <random>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::detail::Test;
using namespace Acts::UnitLiterals;

namespace ActsTests {

static const auto pion = ParticleHypothesis::pion();

using TrackContainer =
    TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                   detail::ValueHolder>;
using TrackStateContainerBackend =
    typename TrackContainer::TrackStateContainerBackend;

struct Detector {
  // expected number of measurements for the given detector
  std::size_t numMeasurements = 6u;

  // geometry
  CubicTrackingGeometry store;
  std::shared_ptr<const TrackingGeometry> geometry;

  // resolutions
  MeasurementResolution resPixel = {MeasurementType::eLoc01, {25_um, 50_um}};
  MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
  MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
  MeasurementResolutionMap resolutions = {
      {GeometryIdentifier().withVolume(2), resPixel},
      {GeometryIdentifier().withVolume(3).withLayer(2), resStrip0},
      {GeometryIdentifier().withVolume(3).withLayer(4), resStrip1},
      {GeometryIdentifier().withVolume(3).withLayer(6), resStrip0},
      {GeometryIdentifier().withVolume(3).withLayer(8), resStrip1},
  };

  explicit Detector(const GeometryContext& geoCtx)
      : store(geoCtx), geometry(store()) {}
};

struct Fixture {
  using StraightPropagator = Propagator<StraightLineStepper, Navigator>;
  using ConstantFieldStepper = EigenStepper<>;
  using ConstantFieldPropagator = Propagator<ConstantFieldStepper, Navigator>;

  using KalmanUpdater = GainMatrixUpdater;
  using KalmanSmoother = GainMatrixSmoother;
  using TestCombinatorialKalmanFilter =
      CombinatorialKalmanFilter<ConstantFieldPropagator, TrackContainer>;
  using TestSourceLinkContainer =
      std::unordered_multimap<GeometryIdentifier, TestSourceLink>;
  using TestCombinatorialKalmanFilterOptions =
      CombinatorialKalmanFilterOptions<TrackContainer>;

  KalmanUpdater kfUpdater;
  KalmanSmoother kfSmoother;

  GeometryContext geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext magCtx;
  CalibrationContext calCtx;

  Detector detector;

  // track parameters before and after the detector
  std::vector<BoundTrackParameters> startParameters;
  std::vector<BoundTrackParameters> endParameters;

  // generated measurements
  TestSourceLinkContainer sourceLinks;

  // CKF implementation to be tested
  TestCombinatorialKalmanFilter ckf;
  CombinatorialKalmanFilterExtensions<TrackContainer> getExtensions() const {
    CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
    extensions.updater.template connect<
        &KalmanUpdater::operator()<TrackStateContainerBackend>>(&kfUpdater);
    return extensions;
  }

  std::unique_ptr<const Logger> logger;

  explicit Fixture(double bz)
      : detector(geoCtx),
        ckf(makeConstantFieldPropagator(detector.geometry, bz)),
        logger(getDefaultLogger("CkfTest", Logging::INFO)) {
    // construct initial parameters
    // create common covariance matrix from reasonable standard deviations
    BoundVector stddev;
    stddev[eBoundLoc0] = 100_um;
    stddev[eBoundLoc1] = 100_um;
    stddev[eBoundTime] = 25_ns;
    stddev[eBoundPhi] = 2_degree;
    stddev[eBoundTheta] = 2_degree;
    stddev[eBoundQOverP] = 1 / 100_GeV;
    BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
    // all tracks close to the transverse plane along the x axis w/ small
    // variations in position, direction.
    Vector4 mStartPos0(-3_m, 0.0, 0.0, 1_ns);
    Vector4 mStartPos1(-3_m, -15_mm, -15_mm, 2_ns);
    Vector4 mStartPos2(-3_m, 15_mm, 15_mm, -1_ns);
    startParameters = {
        BoundTrackParameters::createCurvilinear(mStartPos0, 0_degree, 90_degree,
                                                1_e / 1_GeV, cov, pion),
        BoundTrackParameters::createCurvilinear(
            mStartPos1, -1_degree, 91_degree, 1_e / 1_GeV, cov, pion),
        BoundTrackParameters::createCurvilinear(mStartPos2, 1_degree, 89_degree,
                                                -1_e / 1_GeV, cov, pion),
    };
    Vector4 mEndPos0(3_m, 0.0, 0.0, 1_ns);
    Vector4 mEndPos1(3_m, -100_mm, -100_mm, 2_ns);
    Vector4 mEndPos2(3_m, 100_mm, 100_mm, -1_ns);
    endParameters = {
        BoundTrackParameters::createCurvilinear(mEndPos0, 0_degree, 90_degree,
                                                1_e / 1_GeV, cov * 100, pion),
        BoundTrackParameters::createCurvilinear(mEndPos1, -1_degree, 91_degree,
                                                1_e / 1_GeV, cov * 100, pion),
        BoundTrackParameters::createCurvilinear(mEndPos2, 1_degree, 89_degree,
                                                -1_e / 1_GeV, cov * 100, pion),
    };

    // create some measurements
    auto measPropagator = makeStraightPropagator(detector.geometry);
    std::default_random_engine rng(421235);
    for (std::size_t trackId = 0u; trackId < startParameters.size();
         ++trackId) {
      auto measurements = createMeasurements(
          measPropagator, geoCtx, magCtx, startParameters[trackId],
          detector.resolutions, rng, trackId);
      for (auto& sl : measurements.sourceLinks) {
        sourceLinks.emplace(sl.m_geometryId, std::move(sl));
      }
    }
  }

  // Construct a straight-line propagator.
  static StraightPropagator makeStraightPropagator(
      std::shared_ptr<const TrackingGeometry> geo) {
    Navigator::Config cfg{std::move(geo)};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Navigator navigator{cfg};
    StraightLineStepper stepper;
    return StraightPropagator(stepper, std::move(navigator));
  }

  // Construct a propagator using a constant magnetic field along z.
  static ConstantFieldPropagator makeConstantFieldPropagator(
      std::shared_ptr<const TrackingGeometry> geo, double bz) {
    Navigator::Config cfg{std::move(geo)};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Navigator navigator{cfg};
    auto field = std::make_shared<ConstantBField>(Vector3(0.0, 0.0, bz));
    ConstantFieldStepper stepper(std::move(field));
    return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
  }

  TestCombinatorialKalmanFilterOptions makeCkfOptions() const {
    // leave the accessor empty, this will have to be set before running the CKF
    return CombinatorialKalmanFilterOptions(
        geoCtx, magCtx, calCtx, getExtensions(),
        PropagatorPlainOptions(geoCtx, magCtx));
  }
};

/// A calibrated measurement reading straight out of a @ref TestSourceLink,
/// which stores its values inline in a 2d buffer. Its dimension is only known
/// at runtime, i.e. this exercises the dynamically sized side of
/// @c Acts::detail::CalibratedMeasurementConcept.
class CalibratedTestSourceLink {
 public:
  explicit CalibratedTestSourceLink(const TestSourceLink& sourceLink)
      : m_sourceLink(&sourceLink) {
    for (const BoundIndices index : sourceLink.indices) {
      if (index == eBoundSize) {
        break;
      }
      m_indices[m_size] = static_cast<std::uint8_t>(index);
      ++m_size;
    }
    if (m_size == 0) {
      throw std::runtime_error(
          "Tried to extract measurement from invalid TestSourceLink");
    }
  }

  std::size_t size() const { return m_size; }

  std::span<const std::uint8_t> subspaceIndices() const {
    return {m_indices.data(), m_size};
  }

  Eigen::Map<const Eigen::VectorXd> parameters() const {
    return {m_sourceLink->parameters.data(), static_cast<Eigen::Index>(m_size)};
  }

  Eigen::Map<const Eigen::MatrixXd> covariance() const {
    return {m_sourceLink->covariance.data(), static_cast<Eigen::Index>(m_size),
            static_cast<Eigen::Index>(m_size)};
  }

 private:
  const TestSourceLink* m_sourceLink{};
  std::size_t m_size{0};
  std::array<std::uint8_t, 2> m_indices{};
};

static_assert(
    Acts::detail::CalibratedMeasurementConcept<CalibratedTestSourceLink> &&
        !Acts::detail::StaticCalibratedMeasurementConcept<
            CalibratedTestSourceLink>,
    "The calibrated test source link should exercise the runtime "
    "dispatch on the measurement dimension");

/// The same thing built on @ref CkfTrackStateCreatorBase, i.e. without pushing
/// every measurement through the track EDM first.
struct TestTrackStateCreator final
    : public CkfTrackStateCreatorBase<TestTrackStateCreator, TrackContainer> {
  using Base = CkfTrackStateCreatorBase<TestTrackStateCreator, TrackContainer>;
  using State = typename Base::State;

  const Fixture::TestSourceLinkContainer* container = nullptr;

  /// Cuts the selection below applies
  std::size_t numMeasurements = 1;
  double chi2Measurement = std::numeric_limits<double>::max();
  double chi2Outlier = std::numeric_limits<double>::max();

  auto measurementRange(const State& state) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(state.surface->geometryId());
    return std::ranges::subrange(begin, end) | std::views::values;
  }

  SourceLink sourceLink(const State& /*state*/,
                        const TestSourceLink& measurement) const {
    return SourceLink{measurement};
  }

  CalibratedTestSourceLink calibrate(const State& /*state*/,
                                     const TestSourceLink& measurement) const {
    return CalibratedTestSourceLink{measurement};
  }

  /// Same policy as the legacy creator, but selecting on the measurements
  /// directly instead of on temporary track states.
  template <typename range_t, typename candidates_t>
  Result<void> selectMeasurements(const State& state, range_t&& range,
                                  candidates_t& candidates) const {
    for (const TestSourceLink& measurement : range) {
      candidates.push_back(
          {measurement,
           calculatePredictedChi2(state, calibrate(state, measurement)),
           false});
    }

    const auto best = std::ranges::min_element(
        candidates, {}, [](const auto& candidate) { return candidate.chi2; });
    if (best == candidates.end()) {
      return Result<void>::success();
    }

    const bool isOutlier = best->chi2 >= chi2Measurement;
    if (isOutlier && best->chi2 >= chi2Outlier) {
      candidates.clear();
      return Result<void>::success();
    }
    if (isOutlier || numMeasurements == 1) {
      const auto selected = *best;
      candidates.clear();
      candidates.push_back({selected.measurement, selected.chi2, isOutlier});
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

/// Returns a track state without any flags on every surface with
/// measurements, which the CKF cannot process.
struct FlaglessTrackStateCreator {
  const Fixture::TestSourceLinkContainer* container = nullptr;

  Result<CkfTypes::BranchVector<TrackIndexType>> createTrackStates(
      const CkfTypes::TrackStateCreatorInput<TrackContainer>& input) const {
    CkfTypes::BranchVector<TrackIndexType> trackStates;
    if (!container->contains(input.surface->geometryId())) {
      return trackStates;
    }

    auto trackState = input.trajectory->makeTrackState(
        TrackStatePropMask::Predicted | TrackStatePropMask::Jacobian,
        input.prevTip);
    const auto& params = std::get<BoundTrackParameters>(*input.boundState);
    trackState.predicted() = params.parameters();
    trackState.predictedCovariance() = *params.covariance();
    trackState.jacobian() = std::get<BoundMatrix>(*input.boundState);
    trackState.setReferenceSurface(input.surface->getSharedPtr());

    trackStates.push_back(trackState.index());
    return trackStates;
  }
};

/// Hides the measurements on some surfaces, so that the CKF passes surfaces
/// without a measurement. It uses the default selection of the base.
struct HoleTestCreator final
    : public CkfTrackStateCreatorBase<HoleTestCreator, TrackContainer> {
  using Base = CkfTrackStateCreatorBase<HoleTestCreator, TrackContainer>;
  using State = typename Base::State;

  /// Who creates the track state on a surface without a measurement
  enum class HoleMode {
    /// the creator, through the default of the base
    Creator,
    /// the CKF, because the creator returns no track state
    Ckf,
    /// the CKF, because the creator returns `NoMeasurementExpected`
    CkfNoMeasurementExpected,
  };

  const Fixture::TestSourceLinkContainer* container = nullptr;
  /// The measurements on these surfaces are hidden
  std::vector<GeometryIdentifier> hiddenSurfaces;
  HoleMode holeMode = HoleMode::Creator;
  bool measurementExpected = true;

  auto measurementRange(const State& state) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(state.surface->geometryId());
    if (std::ranges::find(hiddenSurfaces, state.surface->geometryId()) !=
        hiddenSurfaces.end()) {
      begin = end;
    }
    return std::ranges::subrange(begin, end) | std::views::values;
  }

  SourceLink sourceLink(const State& /*state*/,
                        const TestSourceLink& measurement) const {
    return SourceLink{measurement};
  }

  CalibratedTestSourceLink calibrate(const State& /*state*/,
                                     const TestSourceLink& measurement) const {
    return CalibratedTestSourceLink{measurement};
  }

  Result<void> createHoleTrackStates(
      State& state, CkfTypes::BranchVector<TrackIndexType>& trackStates) const {
    switch (holeMode) {
      case HoleMode::Creator:
        return Base::createHoleTrackStates(state, trackStates);
      case HoleMode::Ckf:
        return Result<void>::success();
      case HoleMode::CkfNoMeasurementExpected:
        return CombinatorialKalmanFilterError::NoMeasurementExpected;
    }
    return Result<void>::success();
  }

  bool isMeasurementExpected(const State& /*state*/) const {
    return measurementExpected;
  }
};

}  // namespace ActsTests

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFindingSuite)

BOOST_AUTO_TEST_CASE(ZeroFieldForward) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  // this is the default option. set anyway for consistency
  options.propagatorPlainOptions.direction = Direction::Forward();
  // Construct a plane surface as the target surface
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(Vector3{-3_m, 0., 0.}, Vector3{1., 0., 0})
          .planeSurface();

  TestTrackStateCreator trackStateCreator;
  trackStateCreator.container = &f.sourceLinks;

  options.extensions.trackStateCreator
      .template connect<&TestTrackStateCreator::createTrackStates>(
          &trackStateCreator);

  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

  // run the CKF for all initial track states
  for (std::size_t trackId = 0u; trackId < f.startParameters.size();
       ++trackId) {
    auto res = f.ckf.findTracks(f.startParameters.at(trackId), options, tc);
    if (!res.ok()) {
      BOOST_TEST_INFO(res.error() << " " << res.error().message());
    }
    BOOST_REQUIRE(res.ok());
  }

  // There should be three track finding results with three initial track states
  BOOST_CHECK_EQUAL(tc.size(), 3u);

  // check the found tracks
  for (std::size_t trackId = 0u; trackId < f.startParameters.size();
       ++trackId) {
    const auto track = tc.getTrack(trackId);
    const auto& params = f.startParameters[trackId];
    BOOST_TEST_INFO("initial parameters before detector:\n" << params);

    BOOST_CHECK_EQUAL(track.nTrackStates(), f.detector.numMeasurements);

    // check purity of first found track
    // find the number of hits not originating from the right track
    std::size_t numHits = 0u;
    std::size_t nummismatchedHits = 0u;
    for (const auto trackState : track.trackStatesReversed()) {
      numHits += 1u;
      auto sl =
          trackState.getUncalibratedSourceLink().template get<TestSourceLink>();
      if (trackId != sl.sourceId) {
        nummismatchedHits++;
      }
    }

    BOOST_CHECK_EQUAL(numHits, f.detector.numMeasurements);
    BOOST_CHECK_EQUAL(nummismatchedHits, 0u);
  }
}

/// A track state which the CKF cannot process must stop the track finding.
/// Before, the CKF dropped the branch without an error.
BOOST_AUTO_TEST_CASE(UnsupportedTrackStateFails) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();

  FlaglessTrackStateCreator trackStateCreator;
  trackStateCreator.container = &f.sourceLinks;

  options.extensions.trackStateCreator
      .template connect<&FlaglessTrackStateCreator::createTrackStates>(
          &trackStateCreator);

  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

  auto res = f.ckf.findTracks(f.startParameters.front(), options, tc);

  BOOST_REQUIRE(!res.ok());
  BOOST_CHECK_EQUAL(res.error(),
                    CombinatorialKalmanFilterError::UnsupportedTrackState);
}

namespace {

/// Run the CKF for all start parameters. If `hide` is set, the measurements on
/// the second and the fourth surface with measurements are hidden.
TrackContainer findWithHiddenSurfaces(Fixture& f,
                                      HoleTestCreator::HoleMode mode,
                                      bool measurementExpected,
                                      bool hide = true) {
  std::vector<GeometryIdentifier> surfaces;
  for (const auto& [id, sl] : f.sourceLinks) {
    if (std::ranges::find(surfaces, id) == surfaces.end()) {
      surfaces.push_back(id);
    }
  }
  std::sort(surfaces.begin(), surfaces.end());
  BOOST_REQUIRE_GE(surfaces.size(), 4u);

  HoleTestCreator trackStateCreator;
  trackStateCreator.container = &f.sourceLinks;
  if (hide) {
    trackStateCreator.hiddenSurfaces = {surfaces[1], surfaces[3]};
  }
  trackStateCreator.holeMode = mode;
  trackStateCreator.measurementExpected = measurementExpected;

  auto options = f.makeCkfOptions();
  options.extensions.trackStateCreator
      .template connect<&HoleTestCreator::createTrackStates>(
          &trackStateCreator);

  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  for (const auto& startParameters : f.startParameters) {
    auto res = f.ckf.findTracks(startParameters, options, tc);
    BOOST_REQUIRE(res.ok());
  }
  return tc;
}

/// The track states of both containers must be the same
void checkSameTracks(const TrackContainer& a, const TrackContainer& b) {
  BOOST_REQUIRE_EQUAL(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto ta = a.getTrack(i);
    const auto tb = b.getTrack(i);
    BOOST_CHECK_EQUAL(ta.nMeasurements(), tb.nMeasurements());
    BOOST_CHECK_EQUAL(ta.nOutliers(), tb.nOutliers());
    BOOST_CHECK_EQUAL(ta.nHoles(), tb.nHoles());
    BOOST_CHECK_EQUAL(ta.chi2(), tb.chi2());
    BOOST_REQUIRE_EQUAL(ta.nTrackStates(), tb.nTrackStates());

    auto sa = ta.trackStatesReversed().begin();
    auto sb = tb.trackStatesReversed().begin();
    for (; sa != ta.trackStatesReversed().end(); ++sa, ++sb) {
      const auto stateA = *sa;
      const auto stateB = *sb;
      BOOST_CHECK_EQUAL(stateA.typeFlags().raw(), stateB.typeFlags().raw());
      BOOST_CHECK_EQUAL(&stateA.referenceSurface(), &stateB.referenceSurface());
      BOOST_CHECK_EQUAL(stateA.pathLength(), stateB.pathLength());
      BOOST_CHECK(stateA.predicted() == stateB.predicted());
      BOOST_CHECK(stateA.predictedCovariance() == stateB.predictedCovariance());
      BOOST_CHECK(stateA.filtered() == stateB.filtered());
      BOOST_CHECK(stateA.filteredCovariance() == stateB.filteredCovariance());
      BOOST_CHECK(stateA.jacobian() == stateB.jacobian());
    }
  }
}

}  // namespace

/// A hole from the creator gives the same tracks as a hole from the CKF.
BOOST_AUTO_TEST_CASE(CreatorHolesMatchCkfHoles) {
  Fixture f(0_T);

  auto fromCreator =
      findWithHiddenSurfaces(f, HoleTestCreator::HoleMode::Creator, true);
  auto fromCkf =
      findWithHiddenSurfaces(f, HoleTestCreator::HoleMode::Ckf, true);
  auto reference = findWithHiddenSurfaces(f, HoleTestCreator::HoleMode::Creator,
                                          true, false);

  // the two hidden surfaces add two holes to each track
  BOOST_REQUIRE_EQUAL(fromCreator.size(), reference.size());
  for (std::size_t i = 0; i < fromCreator.size(); ++i) {
    BOOST_CHECK_EQUAL(fromCreator.getTrack(i).nHoles(),
                      reference.getTrack(i).nHoles() + 2u);
  }
  checkSameTracks(fromCreator, fromCkf);
}

/// If no measurement is expected, the creator flags the track state with no
/// expected hit, which does not count as a hole. This is the same as the
/// `NoMeasurementExpected` error.
BOOST_AUTO_TEST_CASE(CreatorNoExpectedHitMatchesCkf) {
  Fixture f(0_T);

  auto fromCreator =
      findWithHiddenSurfaces(f, HoleTestCreator::HoleMode::Creator, false);
  auto fromCkf = findWithHiddenSurfaces(
      f, HoleTestCreator::HoleMode::CkfNoMeasurementExpected, false);

  auto reference = findWithHiddenSurfaces(f, HoleTestCreator::HoleMode::Creator,
                                          true, false);

  // the hidden surfaces and the holes without them all have no expected hit
  BOOST_REQUIRE_EQUAL(fromCreator.size(), reference.size());
  for (std::size_t i = 0; i < fromCreator.size(); ++i) {
    const auto track = fromCreator.getTrack(i);
    BOOST_CHECK_EQUAL(track.nHoles(), 0u);
    std::size_t nNoExpectedHit = 0;
    for (const auto trackState : track.trackStatesReversed()) {
      if (trackState.typeFlags().hasNoExpectedHit()) {
        ++nNoExpectedHit;
      }
    }
    BOOST_CHECK_EQUAL(nNoExpectedHit, reference.getTrack(i).nHoles() + 2u);
  }
  checkSameTracks(fromCreator, fromCkf);
}

/// After a branch switch, the stepper counts its path length from zero again.
/// The path length of the track states must still grow along each track.
BOOST_AUTO_TEST_CASE(PathLengthGrowsAfterBranchSwitch) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  // two branches on every surface make about 64 branches
  options.propagatorPlainOptions.maxSteps = 100000;

  TestTrackStateCreator trackStateCreator;
  trackStateCreator.container = &f.sourceLinks;
  trackStateCreator.numMeasurements = 2;
  options.extensions.trackStateCreator
      .template connect<&TestTrackStateCreator::createTrackStates>(
          &trackStateCreator);

  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto res = f.ckf.findTracks(f.startParameters.front(), options, tc);
  BOOST_REQUIRE(res.ok());
  BOOST_REQUIRE_GT(res->size(), 1u);

  for (const auto& track : *res) {
    BOOST_CHECK_EQUAL(track.nTrackStates(), f.detector.numMeasurements);
    double next = std::numeric_limits<double>::infinity();
    for (const auto trackState : track.trackStatesReversed()) {
      // the closest planes are 2 mm apart along the track
      BOOST_CHECK_LT(trackState.pathLength(), next - 1_mm);
      next = trackState.pathLength();
    }
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldBackward) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  options.propagatorPlainOptions.direction = Direction::Backward();
  // Construct a plane surface as the target surface
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(Vector3{3_m, 0., 0.}, Vector3{1., 0., 0})
          .planeSurface();

  TestTrackStateCreator trackStateCreator;
  trackStateCreator.container = &f.sourceLinks;

  options.extensions.trackStateCreator
      .template connect<&TestTrackStateCreator::createTrackStates>(
          &trackStateCreator);

  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};

  // run the CKF for all initial track states
  for (std::size_t trackId = 0u; trackId < f.startParameters.size();
       ++trackId) {
    auto res = f.ckf.findTracks(f.endParameters.at(trackId), options, tc);
    if (!res.ok()) {
      BOOST_TEST_INFO(res.error() << " " << res.error().message());
    }
    BOOST_REQUIRE(res.ok());
  }
  // There should be three found tracks with three initial track states
  BOOST_CHECK_EQUAL(tc.size(), 3u);

  // check the found tracks
  for (std::size_t trackId = 0u; trackId < f.endParameters.size(); ++trackId) {
    const auto track = tc.getTrack(trackId);
    const auto& params = f.endParameters[trackId];
    BOOST_TEST_INFO("initial parameters after detector:\n" << params);

    BOOST_CHECK_EQUAL(track.nTrackStates(), f.detector.numMeasurements);

    // check purity of first found track
    // find the number of hits not originating from the right track
    std::size_t numHits = 0u;
    std::size_t nummismatchedHits = 0u;
    for (const auto trackState : track.trackStatesReversed()) {
      numHits += 1u;
      auto sl =
          trackState.getUncalibratedSourceLink().template get<TestSourceLink>();
      if (trackId != sl.sourceId) {
        nummismatchedHits++;
      }
    }

    BOOST_CHECK_EQUAL(numHits, f.detector.numMeasurements);
    BOOST_CHECK_EQUAL(nummismatchedHits, 0u);
  }
}

/// Puts material on the layer approach and sensitive surfaces
class ApproachMaterialDecorator final : public IMaterialDecorator {
 public:
  void decorate(Surface& surface) const override {
    const GeometryIdentifier geoId = surface.geometryId();
    if (geoId.approach() == 0 && geoId.sensitive() == 0) {
      return;
    }
    surface.assignSurfaceMaterial(m_material);
  }
  void decorate(TrackingVolume& /*volume*/) const override {}

 private:
  std::shared_ptr<const ISurfaceMaterial> m_material =
      std::make_shared<const HomogeneousSurfaceMaterial>(
          MaterialSlab(makeSilicon(), 5_mm));
};

/// Telescope along x whose layer approach surfaces are material-only surfaces,
/// which `CubicTrackingGeometry` does not have
struct MaterialTelescope {
  static constexpr std::size_t nSensitive = 4u;

  std::shared_ptr<const TrackingGeometry> geometry;

  explicit MaterialTelescope(const GeometryContext& geoCtx) {
    const double halfSize = 1_m;
    auto rBounds = std::make_shared<const RectangleBounds>(halfSize, halfSize);

    // rotate the plane normals onto the x axis
    RotationMatrix3 rotation = RotationMatrix3::Identity();
    rotation.col(0) = Vector3(0., 0., 1.);
    rotation.col(1) = Vector3(0., 1., 0.);
    rotation.col(2) = Vector3(-1., 0., 0.);

    std::vector<CuboidVolumeBuilder::LayerConfig> layerConfigs;
    for (std::size_t i = 1; i <= nSensitive; ++i) {
      CuboidVolumeBuilder::SurfaceConfig sCfg;
      sCfg.position = {i * 1_m, 0., 0.};
      sCfg.rotation = rotation;
      sCfg.rBounds = rBounds;
      sCfg.thickness = 1_um;
      sCfg.detElementConstructor =
          [](const Transform3& trafo,
             const std::shared_ptr<const RectangleBounds>& bounds,
             double thickness) {
            return new DetectorElementStub(trafo, bounds, thickness);
          };

      CuboidVolumeBuilder::LayerConfig lCfg;
      lCfg.surfaceCfg = {sCfg};
      lCfg.active = true;
      lCfg.envelopeX = {-1_mm, 1_mm};
      lCfg.envelopeY = {-0.1_mm, 0.1_mm};
      lCfg.envelopeZ = {-0.1_mm, 0.1_mm};
      layerConfigs.push_back(lCfg);
    }

    const double length = (nSensitive + 1) * 1_m;
    CuboidVolumeBuilder::VolumeConfig vCfg;
    vCfg.length = {length, 2 * halfSize, 2 * halfSize};
    vCfg.position = {length / 2, 0., 0.};
    vCfg.layerCfg = layerConfigs;
    vCfg.name = "MaterialTelescope";

    CuboidVolumeBuilder cvb;
    CuboidVolumeBuilder::Config cfg;
    cfg.length = vCfg.length;
    cfg.position = vCfg.position;
    cfg.volumeCfg = {vCfg};
    cvb.setConfig(cfg);

    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.materialDecorator = std::make_shared<ApproachMaterialDecorator>();
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return cvb.trackingVolume(context, inner, nullptr);
        });
    geometry = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);
  }
};

/// Fixture running the CKF over `MaterialTelescope` with a single track
struct MaterialFixture {
  using StraightPropagator = Fixture::StraightPropagator;
  using ConstantFieldStepper = Fixture::ConstantFieldStepper;
  using ConstantFieldPropagator = Fixture::ConstantFieldPropagator;
  using KalmanUpdater = Fixture::KalmanUpdater;
  using TestCombinatorialKalmanFilter = Fixture::TestCombinatorialKalmanFilter;
  using TestSourceLinkContainer = Fixture::TestSourceLinkContainer;
  using BremPropagator =
      Propagator<MultiStepperLoop<ConstantFieldStepper>, Navigator>;
  using BremCombinatorialKalmanFilter =
      CombinatorialKalmanFilter<BremPropagator, TrackContainer>;

  GeometryContext geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext magCtx;
  CalibrationContext calCtx;

  MaterialTelescope detector{geoCtx};

  BoundTrackParameters startParameters =
      BoundTrackParameters::createCurvilinear(
          Vector4(0.1_m, 0., 0., 0.), 0_degree, 90_degree, 1_e / 1_GeV,
          []() {
            BoundVector stddev;
            stddev[eBoundLoc0] = 100_um;
            stddev[eBoundLoc1] = 100_um;
            stddev[eBoundTime] = 25_ns;
            stddev[eBoundPhi] = 2_degree;
            stddev[eBoundTheta] = 2_degree;
            stddev[eBoundQOverP] = 1 / 100_GeV;
            return BoundMatrix(stddev.cwiseProduct(stddev).asDiagonal());
          }(),
          pion);

  TestSourceLinkContainer sourceLinks;

  KalmanUpdater kfUpdater;

  TestCombinatorialKalmanFilter ckf;
  BremCombinatorialKalmanFilter bremCkf;

  static Navigator makeNavigator(std::shared_ptr<const TrackingGeometry> geo) {
    Navigator::Config cfg{std::move(geo)};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    return Navigator{cfg};
  }

  static std::shared_ptr<ConstantBField> makeField() {
    return std::make_shared<ConstantBField>(Vector3(0., 0., 0.));
  }

  MaterialFixture()
      : ckf(ConstantFieldPropagator(ConstantFieldStepper(makeField()),
                                    makeNavigator(detector.geometry))),
        bremCkf(BremPropagator(BremPropagator::Stepper(makeField()),
                               makeNavigator(detector.geometry))) {
    MeasurementResolution res = {MeasurementType::eLoc01, {25_um, 50_um}};
    MeasurementResolutionMap resolutions = {{GeometryIdentifier(), res}};

    Navigator::Config navCfg{detector.geometry};
    navCfg.resolvePassive = false;
    navCfg.resolveMaterial = true;
    navCfg.resolveSensitive = true;
    StraightPropagator measPropagator(StraightLineStepper(), Navigator{navCfg});

    std::default_random_engine rng(421235);
    auto measurements = createMeasurements(
        measPropagator, geoCtx, magCtx, startParameters, resolutions, rng, 0);
    for (auto& sl : measurements.sourceLinks) {
      sourceLinks.emplace(sl.m_geometryId, std::move(sl));
    }
  }

  CombinatorialKalmanFilterExtensions<TrackContainer> getExtensions() {
    CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
    extensions.updater.template connect<
        &KalmanUpdater::operator()<TrackStateContainerBackend>>(&kfUpdater);
    return extensions;
  }

  /// Run the CKF once and return the container holding the single found track
  TrackContainer find(bool recordMaterialStates) {
    auto options = CombinatorialKalmanFilterOptions(
        geoCtx, magCtx, calCtx, getExtensions(),
        PropagatorPlainOptions(geoCtx, magCtx));
    options.recordMaterialStates = recordMaterialStates;

    // the creator defaults keep one measurement per surface without a chi2
    // cut, which is all this material test needs
    TestTrackStateCreator trackStateCreator;
    trackStateCreator.container = &sourceLinks;
    options.extensions.trackStateCreator
        .template connect<&TestTrackStateCreator::createTrackStates>(
            &trackStateCreator);

    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto res = ckf.findTracks(startParameters, options, tc);
    BOOST_REQUIRE(res.ok());
    BOOST_REQUIRE_EQUAL(tc.size(), 1u);
    return tc;
  }

  /// Run the CKF with bremsstrahlung recovery once and return the container
  /// holding the single found track
  TrackContainer findBrem(bool recordMaterialStates) {
    auto extensions = getExtensions();
    extensions.mixtureReducer.connect<&reduceMixtureWithKLDistance>();

    BremCombinatorialKalmanFilterOptions<TrackContainer> options(
        geoCtx, magCtx, calCtx, extensions,
        PropagatorPlainOptions(geoCtx, magCtx));
    options.recordMaterialStates = recordMaterialStates;
    options.betheHeitlerApprox = std::make_shared<PolynomialBetheHeitlerApprox>(
        makeDefaultBetheHeitlerApprox());

    TestTrackStateCreator trackStateCreator;
    trackStateCreator.container = &sourceLinks;
    options.extensions.trackStateCreator
        .template connect<&TestTrackStateCreator::createTrackStates>(
            &trackStateCreator);

    const BoundTrackParameters electronParameters(
        startParameters.referenceSurface().getSharedPtr(),
        startParameters.parameters(), startParameters.covariance(),
        ParticleHypothesis::electron());

    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto res = bremCkf.findTracks(electronParameters, options, tc);
    BOOST_REQUIRE(res.ok());
    BOOST_REQUIRE_EQUAL(tc.size(), 1u);
    // the found track must keep the (non-default) start hypothesis
    BOOST_CHECK_EQUAL(tc.getTrack(0).particleHypothesis(),
                      electronParameters.particleHypothesis());
    return tc;
  }
};

/// Product of the stored jacobians up to the last measurement, which is the
/// last surface both runs have in common
static BoundMatrix jacobianToLastMeasurement(const auto& track) {
  std::vector<BoundMatrix> jacobians;
  bool seenMeasurement = false;
  for (const auto state : track.trackStatesReversed()) {
    seenMeasurement = seenMeasurement || state.typeFlags().hasMeasurement();
    if (seenMeasurement) {
      jacobians.emplace_back(state.jacobian());
    }
  }
  std::ranges::reverse(jacobians);

  BoundMatrix product = BoundMatrix::Identity();
  for (const auto& jacobian : jacobians) {
    product = jacobian * product;
  }
  return product;
}

BOOST_AUTO_TEST_CASE(MaterialTelescopeGeometry) {
  MaterialFixture f;
  std::size_t nSensitive = 0;
  std::size_t nMaterialOnly = 0;
  f.detector.geometry->visitSurfaces(
      [&](const Surface* surface) {
        BOOST_TEST_MESSAGE("surface " << surface->geometryId()
                                      << " sensitive=" << surface->isSensitive()
                                      << " material=" << surface->hasMaterial()
                                      << " center="
                                      << surface->center(f.geoCtx).transpose());
        if (surface->isSensitive()) {
          ++nSensitive;
        } else if (surface->hasMaterial()) {
          ++nMaterialOnly;
        }
      },
      /*restrictToSensitives=*/false);
  BOOST_CHECK_EQUAL(nSensitive, MaterialTelescope::nSensitive);
  BOOST_CHECK_GT(nMaterialOnly, 0u);
}

BOOST_AUTO_TEST_CASE(MaterialStatesRecordedByDefault) {
  MaterialFixture f;
  auto tc = f.find(true);
  const auto track = tc.getTrack(0);

  BOOST_CHECK_EQUAL(track.nMeasurements(), MaterialTelescope::nSensitive);

  std::size_t nMaterialOnly = 0;
  for (const auto state : track.trackStatesReversed()) {
    const auto flags = state.typeFlags();
    if (!flags.hasMeasurement()) {
      BOOST_CHECK(flags.isMaterial());
      BOOST_CHECK(!flags.isHole());
      BOOST_CHECK(!flags.isOutlier());
      ++nMaterialOnly;
    }
  }
  BOOST_CHECK_GT(nMaterialOnly, 0u);
  BOOST_CHECK_EQUAL(track.nTrackStates(),
                    MaterialTelescope::nSensitive + nMaterialOnly);
}

BOOST_AUTO_TEST_CASE(MaterialStatesSkipped) {
  MaterialFixture f;
  auto tc = f.find(false);
  const auto track = tc.getTrack(0);

  BOOST_CHECK_EQUAL(track.nMeasurements(), MaterialTelescope::nSensitive);
  BOOST_CHECK_EQUAL(track.nTrackStates(), MaterialTelescope::nSensitive);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  for (const auto state : track.trackStatesReversed()) {
    BOOST_CHECK(state.typeFlags().hasMeasurement());
  }
}

// The transport must be the same whether or not the surfaces in between were
// recorded
BOOST_AUTO_TEST_CASE(MaterialStatesSkippedPreserveJacobianChain) {
  MaterialFixture fRecorded;
  MaterialFixture fSkipped;
  auto tcRecorded = fRecorded.find(true);
  auto tcSkipped = fSkipped.find(false);

  const BoundMatrix recorded =
      jacobianToLastMeasurement(tcRecorded.getTrack(0));
  const BoundMatrix skipped = jacobianToLastMeasurement(tcSkipped.getTrack(0));

  CHECK_CLOSE_ABS(recorded, skipped, 1e-6);
}

// Smoothing at the measurement surfaces must be unaffected
BOOST_AUTO_TEST_CASE(MaterialStatesSkippedSmoothingAgrees) {
  auto check = [](auto smoother) {
    MaterialFixture fRecorded;
    MaterialFixture fSkipped;
    auto tcRecorded = fRecorded.find(true);
    auto tcSkipped = fSkipped.find(false);

    auto trackRecorded = tcRecorded.getTrack(0);
    auto trackSkipped = tcSkipped.getTrack(0);
    BOOST_REQUIRE(smoothTrack(fRecorded.geoCtx, trackRecorded,
                              *getDefaultLogger("smoother", Logging::INFO),
                              smoother)
                      .ok());
    BOOST_REQUIRE(smoothTrack(fSkipped.geoCtx, trackSkipped,
                              *getDefaultLogger("smoother", Logging::INFO),
                              smoother)
                      .ok());

    // only the measurement states exist in both tracks
    struct Smoothed {
      GeometryIdentifier geoId{};
      BoundVector parameters = BoundVector::Zero();
      BoundMatrix covariance = BoundMatrix::Zero();
    };
    auto measurementStates = [](const auto& track) {
      std::vector<Smoothed> states;
      for (const auto state : track.trackStatesReversed()) {
        if (state.typeFlags().hasMeasurement()) {
          states.push_back({state.referenceSurface().geometryId(),
                            state.smoothed(), state.smoothedCovariance()});
        }
      }
      return states;
    };

    const auto statesRecorded = measurementStates(trackRecorded);
    const auto statesSkipped = measurementStates(trackSkipped);
    BOOST_REQUIRE_EQUAL(statesRecorded.size(), statesSkipped.size());
    BOOST_REQUIRE_EQUAL(statesRecorded.size(), MaterialTelescope::nSensitive);

    for (std::size_t i = 0; i < statesRecorded.size(); ++i) {
      BOOST_CHECK_EQUAL(statesRecorded[i].geoId, statesSkipped[i].geoId);
      CHECK_CLOSE_ABS(statesRecorded[i].parameters, statesSkipped[i].parameters,
                      1e-8);
      CHECK_CLOSE_ABS(statesRecorded[i].covariance, statesSkipped[i].covariance,
                      1e-8);
    }
  };

  check(GainMatrixSmoother());
  check(MbfSmoother());
}

// With bremsstrahlung recovery, every state carries the transport jacobian of
// the dominant component
BOOST_AUTO_TEST_CASE(BremStatesCarryJacobian) {
  MaterialFixture f;
  auto tc = f.findBrem(true);
  const auto track = tc.getTrack(0);

  BOOST_CHECK_EQUAL(track.nMeasurements(), MaterialTelescope::nSensitive);
  for (const auto state : track.trackStatesReversed()) {
    BOOST_CHECK(!state.jacobian().isZero());
  }
}

BOOST_AUTO_TEST_CASE(BremMaterialStatesSkippedPreserveJacobianChain) {
  MaterialFixture fRecorded;
  MaterialFixture fSkipped;
  auto tcRecorded = fRecorded.findBrem(true);
  auto tcSkipped = fSkipped.findBrem(false);

  const auto trackSkipped = tcSkipped.getTrack(0);
  BOOST_CHECK_EQUAL(trackSkipped.nMeasurements(),
                    MaterialTelescope::nSensitive);
  BOOST_CHECK_EQUAL(trackSkipped.nTrackStates(), MaterialTelescope::nSensitive);

  const BoundMatrix recorded =
      jacobianToLastMeasurement(tcRecorded.getTrack(0));
  const BoundMatrix skipped = jacobianToLastMeasurement(trackSkipped);

  CHECK_CLOSE_ABS(recorded, skipped, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
