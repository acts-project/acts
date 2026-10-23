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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/BetheHeitlerApprox.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/TrackFitting/GsfMixtureReduction.hpp"
#include "Acts/TrackFitting/GsfOptions.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"

#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "FitterTestsCommon.hpp"

using namespace Acts;
using namespace detail::Test;
using namespace UnitLiterals;

namespace ActsTests {

static const auto electron = ParticleHypothesis::electron();

GainMatrixUpdater kfUpdater;

FitterTester tester;

GsfExtensions<VectorMultiTrajectory> getExtensions() {
  GsfExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibratorStrict<VectorMultiTrajectory>>();
  extensions.updater
      .connect<&GainMatrixUpdater::operator()<VectorMultiTrajectory>>(
          &kfUpdater);
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(
          &tester.surfaceAccessor);
  extensions.mixtureReducer.connect<&reduceMixtureWithKLDistance>();
  return extensions;
}

using Stepper = MultiEigenStepperLoop<>;
using Propagator = Propagator<Stepper, Navigator>;
using GSF = GaussianSumFitter<Propagator, VectorMultiTrajectory>;

const GSF gsfZero(makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T),
                  std::make_shared<PolynomialBetheHeitlerApprox>(
                      makeDefaultBetheHeitlerApprox()));

std::default_random_engine rng(42);

auto makeDefaultGsfOptions() {
  GsfOptions<VectorMultiTrajectory> opts{tester.geoCtx, tester.magCtx,
                                         tester.calCtx};
  opts.extensions = getExtensions();
  opts.propagatorPlainOptions =
      PropagatorPlainOptions(tester.geoCtx, tester.magCtx);
  opts.transverseMomentumCut =
      0.0;  // Explicitly disable cut for existing tests
  return opts;
}

// A Helper type to allow us to put the MultiComponentBoundTrackParameters into
// the function so that it can also be used as GenericBoundTrackParameters for
// the MeasurementsCreator
struct MultiCmpsParsInterface : public BoundTrackParameters {
  MultiComponentBoundTrackParameters multi_pars;

  explicit MultiCmpsParsInterface(const MultiComponentBoundTrackParameters& p)
      : BoundTrackParameters(p.merge(ComponentMergeMethod::eMean)),
        multi_pars(p) {}

  explicit operator MultiComponentBoundTrackParameters() const {
    return multi_pars;
  }
};

auto makeParameters() {
  // create covariance matrix from reasonable standard deviations
  BoundVector stddev;
  stddev[eBoundLoc0] = 100_um;
  stddev[eBoundLoc1] = 100_um;
  stddev[eBoundTime] = 25_ns;
  stddev[eBoundPhi] = 2_degree;
  stddev[eBoundTheta] = 2_degree;
  stddev[eBoundQOverP] = 1 / 100_GeV;
  BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();

  // define a track in the transverse plane along x
  Vector4 mPos4(-3_m, 0., 0., 42_ns);
  BoundTrackParameters cp = BoundTrackParameters::createCurvilinear(
      mPos4, 0_degree, 90_degree, 1_e / 1_GeV, cov, electron);

  // Construct bound multi component parameters from curvilinear ones
  BoundVector deltaLOC0 = BoundVector::Zero();
  deltaLOC0[eBoundLoc0] = 0.5_mm;

  BoundVector deltaLOC1 = BoundVector::Zero();
  deltaLOC1[eBoundLoc1] = 0.5_mm;

  BoundVector deltaQOP = BoundVector::Zero();
  deltaQOP[eBoundQOverP] = 1 / 100_GeV;

  std::vector<std::tuple<double, BoundVector, BoundMatrix>> cmps = {
      {0.2, cp.parameters(), cov},
      {0.2, cp.parameters() + deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() + deltaLOC0 - deltaLOC1 - deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 - deltaLOC1 - deltaQOP, cov}};

  return MultiCmpsParsInterface(MultiComponentBoundTrackParameters(
      cp.referenceSurface().getSharedPtr(), cmps,
      [](const auto& cmp) {
        return std::tie(std::get<0>(cmp), std::get<1>(cmp), std::get<2>(cmp));
      },
      electron));
}

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

struct MomentumCutUpdater {
  mutable std::size_t calls = 0;
  std::size_t rejectOnCall = 0;

  Result<void> update(const GeometryContext& geoContext,
                      VectorMultiTrajectory::TrackStateProxy state,
                      const Logger& logger) const {
    auto result =
        kfUpdater.operator()<VectorMultiTrajectory>(geoContext, state, logger);
    if (result.ok() && ++calls == rejectOnCall) {
      state.filtered()[eBoundQOverP] = 1 / 100_MeV;
    }
    return result;
  }
};

struct SecondMeasurementOutlier {
  mutable std::size_t calls = 0;
  bool operator()(VectorMultiTrajectory::ConstTrackStateProxy /*state*/) const {
    return ++calls == 2;
  }
};

class ExhaustingBetheHeitlerApprox final : public BetheHeitlerApprox {
 public:
  explicit ExhaustingBetheHeitlerApprox(std::size_t exhaustOnCall)
      : m_exhaustOnCall(exhaustOnCall) {}

  std::size_t maxComponents() const override { return 1; }
  bool validXOverX0(double /*xOverX0*/) const override { return true; }
  std::span<Component> mixture(double /*xOverX0*/,
                               std::span<Component> cache) const override {
    cache[0] = {++m_calls == m_exhaustOnCall ? 0.0 : 1.0, 1.0, 0.0};
    return cache.first(1);
  }

 private:
  std::size_t m_exhaustOnCall;
  mutable std::size_t m_calls = 0;
};

BOOST_AUTO_TEST_CASE(
    ForwardMomentumExhaustionStartsSmoothingAtLastMeasurement) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  const BoundTrackParameters start = makeParameters();
  std::default_random_engine localRng(42);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, localRng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();
  options.disableAllMaterialHandling = true;
  options.transverseMomentumCut = 500_MeV;
  MomentumCutUpdater updater{0, 3};
  options.extensions.updater.connect<&MomentumCutUpdater::update>(&updater);

  auto result = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                            options, tracks);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->nMeasurements(), 2u);
  for (const auto state : result->trackStatesReversed()) {
    if (state.typeFlags().isMeasurement()) {
      BOOST_CHECK(state.hasSmoothed());
      BOOST_CHECK(state.smoothed().allFinite());
    }
  }
}

BOOST_AUTO_TEST_CASE(MomentumExhaustionBeforeFirstMeasurementFails) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  const BoundTrackParameters start = makeParameters();
  std::default_random_engine localRng(42);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, localRng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();
  options.disableAllMaterialHandling = true;
  options.transverseMomentumCut = 500_MeV;
  MomentumCutUpdater updater{0, 1};
  options.extensions.updater.connect<&MomentumCutUpdater::update>(&updater);

  auto result = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                            options, tracks);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == GsfError::NoComponentsLeft);
}

BOOST_AUTO_TEST_CASE(ForwardExhaustionAfterOutlierKeepsTrailingOutlier) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  const BoundTrackParameters start = makeParameters();
  std::default_random_engine localRng(42);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, localRng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();
  options.disableAllMaterialHandling = true;
  options.transverseMomentumCut = 500_MeV;
  MomentumCutUpdater updater{0, 2};
  SecondMeasurementOutlier outlier;
  options.extensions.updater.connect<&MomentumCutUpdater::update>(&updater);
  options.extensions.outlierFinder
      .connect<&SecondMeasurementOutlier::operator()>(&outlier);
  auto result = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), start,
                            options, tracks);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(result->nMeasurements(), 1u);
  BOOST_CHECK_EQUAL(result->nOutliers(), 1u);
}

BOOST_AUTO_TEST_CASE(MomentumCutRemovesOnlyRejectedComponentsWithoutMaterial) {
  using Actor = detail::Gsf::GsfActor<VectorMultiTrajectory>;
  VectorMultiTrajectory trajectory;
  const auto start = makeParameters();
  std::default_random_engine localRng(42);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, localRng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();
  options.disableAllMaterialHandling = true;
  options.weightCutoff = 0.0;
  options.transverseMomentumCut = 500_MeV;
  MomentumCutUpdater updater{0, 1};
  options.extensions.updater.connect<&MomentumCutUpdater::update>(&updater);

  auto propagator = makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T);
  Propagator::Options<ActorList<Actor>> propagationOptions(tester.geoCtx,
                                                           tester.magCtx);
  auto& actor = propagationOptions.actorList.get<Actor>();
  actor.setOptions(options);
  auto logger = getDefaultLogger("SelectiveMomentumCut", Logging::INFO);
  actor.m_cfg.logger = logger.get();
  std::map<GeometryIdentifier, SourceLink> measurementMap;
  for (const auto& link : sourceLinks) {
    measurementMap.emplace(
        options.extensions.surfaceAccessor(link)->geometryId(), link);
  }
  actor.m_cfg.inputMeasurements = &measurementMap;
  actor.m_cfg.numberMeasurements = measurementMap.size();
  auto state = propagator.makeState(propagationOptions);
  BOOST_REQUIRE(propagator.initialize(state, start.multi_pars, nullptr).ok());
  state.get<Actor::result_type>().fittedStates = &trajectory;
  auto result = propagator.propagate(state);
  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(state.stepping.components.size(), 4u);
}

BOOST_AUTO_TEST_CASE(MaterialExhaustionStopsForwardButFailsBackward) {
  const BoundTrackParameters start = makeParameters();
  std::default_random_engine localRng(42);
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         start, tester.resolutions, localRng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), tester.nMeasurements);

  for (const bool backward : {false, true}) {
    BOOST_TEST_CONTEXT("backward=" << backward) {
      TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
      auto options = makeDefaultGsfOptions();
      options.maxComponents = 1;
      GSF fitter(makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T),
                 std::make_shared<ExhaustingBetheHeitlerApprox>(
                     backward ? tester.nMeasurements + 1 : 3));
      auto result = fitter.fit(sourceLinks.begin(), sourceLinks.end(), start,
                               options, tracks);
      if (backward) {
        BOOST_REQUIRE(!result.ok());
        BOOST_CHECK(result.error() == GsfError::NoComponentsLeft);
      } else {
        BOOST_REQUIRE(result.ok());
        BOOST_CHECK_EQUAL(result->nMeasurements(), 3u);
        for (const auto state : result->trackStatesReversed()) {
          if (state.typeFlags().isMeasurement()) {
            BOOST_CHECK(state.hasSmoothed());
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(MaterialCutUsesTransverseMomentum) {
  const BoundTrackParameters start = makeParameters();
  ExhaustingBetheHeitlerApprox approximation(0);
  auto logger = getDefaultLogger("TransverseMomentumCut", Logging::INFO);
  for (const double theta : {1_degree, 90_degree}) {
    for (const double cut : {0_GeV, 500_MeV}) {
      auto parameters = start.parameters();
      parameters[eBoundQOverP] = 1 / 10_GeV;
      parameters[eBoundTheta] = theta;
      BoundTrackParameters bound(start.referenceSurface().getSharedPtr(),
                                 parameters, start.covariance(), electron);
      std::vector<BetheHeitlerApprox::Component> cache;
      std::vector<GsfComponent> components;
      std::size_t invalid = 0;
      double maxMaterial = 0;
      detail::Gsf::applyBetheHeitler(tester.geoCtx, start.referenceSurface(),
                                     Direction::Forward(), bound, 1.0,
                                     approximation, cache, 0.0, components,
                                     invalid, maxMaterial, *logger, cut);
      BOOST_CHECK_EQUAL(components.size(),
                        cut > 0 && theta == 1_degree ? 0u : 1u);
    }
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldNoSurfaceForward(gsfZero, options, multi_pars, rng, true,
                                        false, false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceForward(gsfZero, options, multi_pars, rng,
                                          true, false, false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceBackward(gsfZero, options, multi_pars, rng,
                                           true, false, false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceBackward(gsfZero, options, multi_pars, rng,
                                           true, false, false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldShuffled(gsfZero, options, multi_pars, rng, true, false,
                                false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto options = makeDefaultGsfOptions();
  auto multi_pars = makeParameters();

  tester.test_ZeroFieldWithHole(gsfZero, options, multi_pars, rng, true, false,
                                false);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  TestOutlierFinder tof{5_mm};
  auto options = makeDefaultGsfOptions();
  options.extensions.outlierFinder
      .connect<&TestOutlierFinder::operator()<VectorMultiTrajectory>>(&tof);

  auto multi_pars = makeParameters();

  tester.test_ZeroFieldWithOutliers(gsfZero, options, multi_pars, rng, true,
                                    false, false);
}

// Smoothed parameters must be present on exactly the measurement states
BOOST_AUTO_TEST_CASE(SmoothedParametersOnMeasurementStates) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};

  auto multi_pars = makeParameters();
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         multi_pars, tester.resolutions, rng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();

  auto res = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), multi_pars,
                         options, tracks);
  BOOST_REQUIRE(res.ok());

  std::size_t nMeasurementStates = 0;
  for (const auto ts : res->trackStatesReversed()) {
    if (!ts.typeFlags().isMeasurement()) {
      BOOST_CHECK(!ts.hasSmoothed());
      continue;
    }

    ++nMeasurementStates;

    BOOST_CHECK(ts.hasSmoothed());
    BOOST_CHECK(ts.smoothed().allFinite());
    BOOST_CHECK(ts.smoothedCovariance().allFinite());
    // the best available parameters must be the smoothed ones
    BOOST_CHECK_EQUAL(ts.parameters().data(), ts.smoothed().data());
  }

  BOOST_CHECK_EQUAL(nMeasurementStates, sourceLinks.size());
}

// Every state carries the jacobian of the segment that ends at it
BOOST_AUTO_TEST_CASE(JacobianOnAllStates) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};

  auto multi_pars = makeParameters();
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         multi_pars, tester.resolutions, rng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();

  auto res = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), multi_pars,
                         options, tracks);
  BOOST_REQUIRE(res.ok());

  for (const auto ts : res->trackStatesReversed()) {
    BOOST_REQUIRE(ts.hasJacobian());
    BOOST_CHECK(ts.jacobian().allFinite());
    BOOST_CHECK_NE(ts.jacobian().determinant(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(WithFinalMultiComponentState) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  using namespace GsfConstants;
  std::string key(kFinalMultiComponentStateColumn);
  tracks.template addColumn<FinalMultiComponentState>(key);

  auto multi_pars = makeParameters();
  auto measurements =
      createMeasurements(tester.simPropagator, tester.geoCtx, tester.magCtx,
                         multi_pars, tester.resolutions, rng);
  auto sourceLinks = tester.prepareSourceLinks(measurements.sourceLinks);
  auto options = makeDefaultGsfOptions();

  // create a boundless target surface near the tracker exit
  Vector3 center(-3._m, 0., 0.);
  Vector3 normal(1., 0., 0.);
  std::shared_ptr<PlaneSurface> targetSurface =
      CurvilinearSurface(center, normal).planeSurface();

  options.referenceSurface = targetSurface.get();

  auto res = gsfZero.fit(sourceLinks.begin(), sourceLinks.end(), multi_pars,
                         options, tracks);

  BOOST_REQUIRE(res.ok());
  BOOST_CHECK(res->template component<FinalMultiComponentState>(
                     kFinalMultiComponentStateColumn)
                  .has_value());
}

BOOST_AUTO_TEST_CASE(OutputTrackParticleHypothesis) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ParticleHypothesis(gsfZero, options, multi_pars, rng);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
