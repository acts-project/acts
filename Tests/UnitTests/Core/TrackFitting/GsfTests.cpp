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
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/MultiComponentTrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"
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
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
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

const GSF gsfZero(
    makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T),
    std::make_shared<AtlasBetheHeitlerApprox>(makeDefaultBetheHeitlerApprox()));

std::default_random_engine rng(42);

auto makeDefaultGsfOptions() {
  GsfOptions<VectorMultiTrajectory> opts{tester.geoCtx, tester.magCtx,
                                         tester.calCtx};
  opts.extensions = getExtensions();
  opts.propagatorPlainOptions =
      PropagatorPlainOptions(tester.geoCtx, tester.magCtx);
  return opts;
}

// A Helper type to allow us to put the MultiComponentBoundTrackParameters into
// the function so that it can also be used as GenericBoundTrackParameters for
// the MeasurementsCreator
struct MultiCmpsParsInterface : public BoundTrackParameters {
  MultiComponentBoundTrackParameters multi_pars;

  explicit MultiCmpsParsInterface(const MultiComponentBoundTrackParameters &p)
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
  deltaQOP[eBoundQOverP] = 0.01_GeV;

  std::vector<std::tuple<double, BoundVector, BoundMatrix>> cmps = {
      {0.2, cp.parameters(), cov},
      {0.2, cp.parameters() + deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() + deltaLOC0 - deltaLOC1 - deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 - deltaLOC1 - deltaQOP, cov}};

  return MultiCmpsParsInterface(MultiComponentBoundTrackParameters(
      cp.referenceSurface().getSharedPtr(), cmps, electron));
}

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

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

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
