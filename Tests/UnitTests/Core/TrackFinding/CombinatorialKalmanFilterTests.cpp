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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsTests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"

#include <cassert>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <random>
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

/// The map(-like) container accessor
template <typename container_t>
struct TestContainerAccessor {
  using Container = container_t;
  using Key = typename container_t::key_type;
  using Value = typename container_t::mapped_type;

  /// This iterator adapter is needed to have the deref operator return a single
  /// source link instead of the map pair <GeometryIdentifier,SourceLink>
  struct Iterator {
    using BaseIterator = typename container_t::const_iterator;

    using iterator_category = typename BaseIterator::iterator_category;
    using value_type = typename BaseIterator::value_type;
    using difference_type = typename BaseIterator::difference_type;
    using pointer = typename BaseIterator::pointer;
    using reference = typename BaseIterator::reference;

    Iterator& operator++() {
      ++m_iterator;
      return *this;
    }

    bool operator==(const Iterator& other) const {
      return m_iterator == other.m_iterator;
    }

    SourceLink operator*() const {
      const auto& sl = m_iterator->second;
      return SourceLink{sl};
    }

    BaseIterator m_iterator;
  };

  // pointer to the container
  const Container* container = nullptr;

  // get the range of elements with requested key
  std::pair<Iterator, Iterator> range(const Surface& surface) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    return {Iterator{begin}, Iterator{end}};
  }
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
  using TestSourceLinkAccessor = TestContainerAccessor<TestSourceLinkContainer>;
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
  // configuration for the measurement selector
  MeasurementSelector::Config measurementSelectorCfg = {
      // global default: no chi2 cut, only one measurement per surface
      {GeometryIdentifier(), {{}, {std::numeric_limits<double>::max()}, {1u}}},
  };

  MeasurementSelector measSel{measurementSelectorCfg};

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
    BoundSquareMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
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

// set up composable track state creator from these components:
//  - source link accessor,
//  - measurement selector
//  - track  state candidate creator
template <typename source_link_accessor_t>
inline auto makeTrackStateCreator(const source_link_accessor_t& slAccessor,
                                  const MeasurementSelector& measSel) {
  using TrackStateCreatorType =
      TrackStateCreator<typename source_link_accessor_t::Iterator,
                        TrackContainer>;
  TrackStateCreatorType trackStateCreator;
  trackStateCreator.sourceLinkAccessor
      .template connect<&source_link_accessor_t::range>(&slAccessor);
  trackStateCreator.calibrator.template connect<
      &testSourceLinkCalibrator<TrackStateContainerBackend>>();
  trackStateCreator.measurementSelector.template connect<
      &MeasurementSelector::select<TrackStateContainerBackend>>(&measSel);
  return trackStateCreator;
}

}  // namespace ActsTests

namespace Acts {

// somehow this is not automatically instantiated
template Result<std::pair<
    std::vector<
        ActsTests::TrackStateContainerBackend::TrackStateProxy>::iterator,
    std::vector<
        ActsTests::TrackStateContainerBackend::TrackStateProxy>::iterator>>
MeasurementSelector::select<ActsTests::TrackStateContainerBackend>(
    std::vector<ActsTests::TrackStateContainerBackend::TrackStateProxy>&, bool&,
    const Logger&) const;
}  // namespace Acts

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

  Fixture::TestSourceLinkAccessor slAccessor;
  slAccessor.container = &f.sourceLinks;

  auto trackStateCreator = makeTrackStateCreator(slAccessor, f.measSel);

  options.extensions.createTrackStates
      .template connect<&decltype(trackStateCreator)::createTrackStates>(
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

BOOST_AUTO_TEST_CASE(ZeroFieldBackward) {
  Fixture f(0_T);

  auto options = f.makeCkfOptions();
  options.propagatorPlainOptions.direction = Direction::Backward();
  // Construct a plane surface as the target surface
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(Vector3{3_m, 0., 0.}, Vector3{1., 0., 0})
          .planeSurface();

  Fixture::TestSourceLinkAccessor slAccessor;
  slAccessor.container = &f.sourceLinks;

  auto trackStateCreator = makeTrackStateCreator(slAccessor, f.measSel);
  options.extensions.createTrackStates
      .template connect<&decltype(trackStateCreator)::createTrackStates>(
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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
