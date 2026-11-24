// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <array>
#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::Test;

namespace ActsTests {
namespace {

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mctx;
const CalibrationContext cctx;
const MaterialSlab silicon(makeSilicon(), 300_um);
const auto pion = ParticleHypothesis::pion();
constexpr double qOverP = 1. / 2_GeV;

RotationMatrix3 telescopeToCanonical() {
  RotationMatrix3 rotation;
  rotation << 0., 1., 0., 0., 0., 1., 1., 0., 0.;
  return rotation;
}

auto plane(const RotationMatrix3& rotation, double z,
           const BoundParametrization& chart) {
  Transform3 transform(rotation);
  transform.translation() = rotation * Vector3(0., 0., z);
  return Surface::makeShared<PlaneSurface>(transform, nullptr, chart);
}

struct TelescopePlanes {
  std::vector<const Surface*> sequence;

  const Surface* find(const SourceLink& sl) const {
    return sequence.at(sl.get<TestSourceLink>().sourceId);
  }
};

auto fitTelescope(const RotationMatrix3& rotation,
                  const BoundParametrization& chart, bool material,
                  bool alternateCharts = false) {
  const std::array<Vector2, 4> hits{
      Vector2(15_um, -25_um), Vector2(-10_um, -5_um), Vector2(22_um, 30_um),
      Vector2(7_um, 12_um)};
  std::vector<std::shared_ptr<PlaneSurface>> planes;
  std::vector<SourceLink> sourceLinks;
  TelescopePlanes telescope;
  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto sensorChart = alternateCharts && i % 2 == 1
                                 ? BoundParametrization::telescope()
                                 : chart;
    auto surface = plane(rotation, (i + 1) * 20_mm, sensorChart);
    if (material && i == 0) {
      // Deposit the entire slab before the first hit, while the track is axial.
      surface->assignSurfaceMaterial(
          std::make_shared<HomogeneousSurfaceMaterial>(silicon, 0.));
    }
    sourceLinks.emplace_back(TestSourceLink(
        eBoundLoc0, eBoundLoc1, hits[i],
        SquareMatrix2::Identity() * 20_um * 20_um, GeometryIdentifier{}, i));
    telescope.sequence.push_back(surface.get());
    planes.push_back(std::move(surface));
  }

  BoundVector sigma;
  sigma << 100_um, 150_um, 2e-3, 3e-3, 0.02, 1_ns;
  BoundMatrix cov = sigma.cwiseAbs2().asDiagonal();
  cov(eBoundLoc0, eBoundDir1) = cov(eBoundDir1, eBoundLoc0) =
      0.25 * sigma[eBoundLoc0] * sigma[eBoundDir1];
  auto seed = BoundTrackParameters::create(
      gctx, plane(rotation, 0., chart), Vector4::Zero(),
      rotation * Vector3::UnitZ(), qOverP, cov, pion);
  BOOST_REQUIRE(seed.ok());

  GainMatrixUpdater updater;
  GainMatrixSmoother smoother;
  KalmanFitterExtensions<VectorMultiTrajectory> extensions;
  extensions.surfaceAccessor.connect<&TelescopePlanes::find>(&telescope);
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  extensions.updater
      .connect<&GainMatrixUpdater::operator()<VectorMultiTrajectory>>(&updater);
  extensions.smoother
      .connect<&GainMatrixSmoother::operator()<VectorMultiTrajectory>>(
          &smoother);
  KalmanFitterOptions options(gctx, mctx, cctx, extensions,
                              PropagatorPlainOptions(gctx, mctx));
  options.energyLoss = false;

  using Prop = Propagator<StraightLineStepper, DirectNavigator>;
  KalmanFitter<Prop, VectorMultiTrajectory> fitter{
      Prop{StraightLineStepper{}, DirectNavigator{}}};
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto result = fitter.fit(sourceLinks.begin(), sourceLinks.end(), *seed,
                           options, telescope.sequence, tracks);
  if (!result.ok()) {
    BOOST_FAIL(result.error().message());
  }
  BOOST_REQUIRE_EQUAL(result->nMeasurements(), hits.size());
  BOOST_REQUIRE_EQUAL(tracks.trackStateContainer().size(), hits.size());
  BOOST_CHECK_GT(result->chi2(), 0.);
  return tracks;
}

}  // namespace

BOOST_AUTO_TEST_CASE(AxialTelescopeKalmanFit) {
  const auto rotation = telescopeToCanonical();
  const auto chart = BoundParametrization::telescope();
  for (double sign : {-1., 1.}) {
    // Turn the entire telescope around its local x axis for the -Z beam.
    const RotationMatrix3 placement = Vector3(1., sign, sign).asDiagonal();
    BoundMatrix vacuumPrediction = BoundMatrix::Zero();
    for (bool material : {false, true}) {
      const auto telescope = fitTelescope(placement, chart, material);
      // A rigid rotation of the complete experiment takes the beam to global
      // ±Y. The reference chart angles and the two local sensor axes then
      // correspond exactly to those of the telescope, including covariance
      // cross terms.
      const auto reference =
          fitTelescope(rotation * placement,
                       BoundParametrization::globalSpherical(), material);
      CHECK_CLOSE_ABS(telescope.getTrack(0).chi2(),
                      reference.getTrack(0).chi2(), 1e-7);
      const auto& states = telescope.trackStateContainer();
      for (TrackIndexType i = 0; i < states.size(); ++i) {
        const auto state = states.getTrackState(i);
        const auto ref = reference.trackStateContainer().getTrackState(i);
        BOOST_REQUIRE(state.hasFiltered() && state.hasSmoothed());
        BOOST_REQUIRE(ref.hasFiltered() && ref.hasSmoothed());
        const auto compare =
            [&](const BoundVector& parameters, const BoundMatrix& cov,
                const BoundVector& refParameters, const BoundMatrix& refCov) {
              const FreeVector free = transformBoundToFreeParameters(
                  state.referenceSurface(), gctx, parameters);
              const FreeVector refFree = transformBoundToFreeParameters(
                  ref.referenceSurface(), gctx, refParameters);
              CHECK_CLOSE_ABS(rotation * free.segment<3>(eFreePos0),
                              refFree.segment<3>(eFreePos0), 1e-9);
              CHECK_CLOSE_ABS(rotation * free.segment<3>(eFreeDir0),
                              refFree.segment<3>(eFreeDir0), 1e-11);
              CHECK_CLOSE_ABS(parameters, refParameters, 1e-10);
              BOOST_REQUIRE(cov.allFinite() && refCov.allFinite());
              CHECK_CLOSE_COVARIANCE(cov, refCov, 1e-9);
            };
        compare(state.filtered(), state.filteredCovariance(), ref.filtered(),
                ref.filteredCovariance());
        compare(state.smoothed(), state.smoothedCovariance(), ref.smoothed(),
                ref.smoothedCovariance());
      }
      const auto first = states.getTrackState(0);
      if (!material) {
        vacuumPrediction = first.predictedCovariance();
      } else {
        // Direct check of the full increment, independent of the fit
        // comparison.
        const double theta0 = computeMultipleScatteringTheta0(
            silicon, pion.absolutePdg(), pion.mass(), qOverP,
            pion.absoluteCharge());
        BoundMatrix increment = BoundMatrix::Zero();
        increment(eBoundDir0, eBoundDir0) = increment(eBoundDir1, eBoundDir1) =
            theta0 * theta0;
        CHECK_CLOSE_ABS(first.predictedCovariance() - vacuumPrediction,
                        increment, 1e-12);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(MixedChartKalmanFit) {
  // Both spherical charts are regular around +Y. Alternate the measurement
  // chart while keeping the seed, geometry, material and hit values identical.
  const auto canonical = BoundParametrization::globalSpherical();
  const auto rotation = telescopeToCanonical();
  const auto mixed = fitTelescope(rotation, canonical, true, true);
  const auto reference = fitTelescope(rotation, canonical, true);
  CHECK_CLOSE_ABS(mixed.getTrack(0).chi2(), reference.getTrack(0).chi2(), 1e-5);
  const auto& states = mixed.trackStateContainer();
  for (TrackIndexType i = 0; i < states.size(); ++i) {
    const auto state = states.getTrackState(i);
    const auto ref = reference.trackStateContainer().getTrackState(i);
    const auto chart = state.referenceSurface().boundParametrization(gctx);
    BOOST_REQUIRE_EQUAL(chart.isGlobalSpherical(), i % 2 == 0);
    const auto compare = [&](const BoundVector& parameters,
                             const BoundMatrix& cov,
                             const BoundVector& referenceParameters,
                             const BoundMatrix& referenceCov) {
      const Vector3 d = chart.direction(parameters);
      BoundVector inCanonical = parameters;
      inCanonical.segment<2>(eBoundDir0) = canonical.directionCoordinates(d);
      // Small, nonzero hit residuals give second-order differences between
      // Kalman corrections applied in the two different nonlinear charts.
      CHECK_CLOSE_ABS(inCanonical, referenceParameters, 1e-7);
      const auto inverse =
          canonical.directionCoordinatesJacobianFromUnitDirection(d);
      BOOST_REQUIRE(inverse.ok());
      BoundMatrix basis = BoundMatrix::Identity();
      basis.block<2, 2>(eBoundDir0, eBoundDir0) =
          *inverse * chart.directionJacobian(parameters);
      const BoundMatrix inCanonicalCov = basis * cov * basis.transpose();
      BOOST_REQUIRE(inCanonicalCov.allFinite());
      CHECK_CLOSE_COVARIANCE(inCanonicalCov, referenceCov, 1e-5);
    };
    compare(state.filtered(), state.filteredCovariance(), ref.filtered(),
            ref.filteredCovariance());
    compare(state.smoothed(), state.smoothedCovariance(), ref.smoothed(),
            ref.smoothedCovariance());
  }
}
}  // namespace ActsTests
