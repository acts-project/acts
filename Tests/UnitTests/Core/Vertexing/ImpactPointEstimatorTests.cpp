// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace {

namespace bd = boost::unit_test::data;

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

using MagneticField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper>;
using Estimator =
    Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

const Acts::GeometryContext geoContext;
const Acts::MagneticFieldContext magFieldContext;

Acts::MagneticFieldProvider::Cache magFieldCache() {
  return NullBField{}.makeCache(magFieldContext);
}

// perigee track parameters dataset
// only non-zero distances are tested
auto d0s = bd::make({-25_um, 25_um});
auto l0s = bd::make({-1_mm, 1_mm});
auto t0s = bd::make({-2_ns, 2_ns});
auto phis = bd::make({0_degree, -45_degree, 45_degree});
auto thetas = bd::make({90_degree, 20_degree, 160_degree});
auto ps = bd::make({0.4_GeV, 1_GeV, 10_GeV});
auto qs = bd::make({-1_e, 1_e});
// Cartesian products over all parameters
auto tracks = d0s * l0s * t0s * phis * thetas * ps * qs;

// vertex parameters dataset
auto vx0s = bd::make({0_um, -10_um, 10_um});
auto vy0s = bd::make({0_um, -10_um, 10_um});
auto vz0s = bd::make({0_mm, -25_mm, 25_mm});
auto vt0s = bd::make({0_ns, -2_ns, 2_ns});
// Cartesian products over all parameters
auto vertices = vx0s * vy0s * vz0s * vt0s;

// Construct an impact point estimator for a constant bfield along z.
Estimator makeEstimator(double bZ) {
  auto field = std::make_shared<MagneticField>(Vector3(0, 0, bZ));
  Stepper stepper(field);
  Estimator::Config cfg(field,
                        std::make_shared<Propagator>(std::move(stepper)));
  return Estimator(cfg);
}

// Construct a diagonal track covariance w/ reasonable values.
Acts::BoundSymMatrix makeBoundParametersCovariance() {
  BoundVector stddev;
  stddev[eBoundLoc0] = 15_um;
  stddev[eBoundLoc1] = 100_um;
  stddev[eBoundTime] = 5_ns;
  stddev[eBoundPhi] = 1_degree;
  stddev[eBoundTheta] = 1_degree;
  stddev[eBoundQOverP] = 1_e / 100_GeV;
  return stddev.cwiseProduct(stddev).asDiagonal();
}

// Construct a diagonal vertex covariance w/ reasonable values.
Acts::SymMatrix4 makeVertexCovariance() {
  Vector4 stddev;
  stddev[ePos0] = 10_um;
  stddev[ePos1] = 10_um;
  stddev[ePos2] = 75_um;
  stddev[eTime] = 1_ns;
  return stddev.cwiseProduct(stddev).asDiagonal();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(VertexingImpactPointEstimator)

// Check `calculated3dDistance`, `estimate3DImpactParameters`, and
// `get3dVertexCompatibility`.
BOOST_DATA_TEST_CASE(SingleTrackDistanceParametersCompatibility3d, tracks, d0,
                     l0, t0, phi, theta, p, q) {
  BoundVector par;
  par[eBoundLoc0] = d0;
  par[eBoundLoc1] = l0;
  par[eBoundTime] = t0;
  par[eBoundPhi] = phi;
  par[eBoundTheta] = theta;
  par[eBoundQOverP] = q / p;

  Estimator ipEstimator = makeEstimator(2_T);
  Estimator::State state(magFieldCache());
  // reference position and corresponding perigee surface
  Vector3 refPosition(0., 0., 0.);
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(refPosition);
  // create the track
  BoundTrackParameters myTrack(perigeeSurface, par,
                               makeBoundParametersCovariance());

  // initial distance to the reference position in the perigee frame
  double distT = std::hypot(d0, l0);
  double dist3 =
      ipEstimator.calculate3dDistance(geoContext, myTrack, refPosition, state)
          .value();
  // estimated 3D distance should be less than the 2d distance in the perigee
  // frame. it should be equal if the track is a transverse track w/ theta =
  // 90deg. in that case there might be numerical deviations and we need to
  // check that it is less or equal within the numerical tolerance.
  BOOST_CHECK((dist3 < distT) or (std::abs(dist3 - distT) < 1_um));

  // estimate parameters at the closest point in 3d
  auto res = ipEstimator.estimate3DImpactParameters(
      geoContext, magFieldContext, myTrack, refPosition, state);
  BoundTrackParameters trackAtIP3d = *res;
  const auto& atPerigee = myTrack.parameters();
  const auto& atIp3d = trackAtIP3d.parameters();

  // all parameters except the helix invariants theta, q/p shoud be changed
  BOOST_CHECK_NE(atPerigee[eBoundLoc0], atIp3d[eBoundLoc0]);
  BOOST_CHECK_NE(atPerigee[eBoundLoc1], atIp3d[eBoundLoc1]);
  // BOOST_CHECK_NE(atPerigee[eBoundTime], atIp3d[eBoundTime]);
  // BOOST_CHECK_NE(atPerigee[eBoundPhi], atIp3d[eBoundPhi]);
  CHECK_CLOSE_ABS(atPerigee[eBoundTheta], atIp3d[eBoundTheta], 0.01_mrad);
  CHECK_CLOSE_REL(atPerigee[eBoundQOverP], atIp3d[eBoundQOverP],
                  std::numeric_limits<ActsScalar>::epsilon());

  // check that we get sensible compatibility scores
  // this is a chi2-like value and should always be positive
  auto compatibility =
      ipEstimator
          .get3dVertexCompatibility(geoContext, &trackAtIP3d, refPosition)
          .value();
  BOOST_CHECK_GT(compatibility, 0);
}

// Compare calculations w/ known good values from Athena.
//
// Checks the results for a single track with the same test values as in Athena
// unit test algorithm
//
//   Tracking/TrkVertexFitter/TrkVertexFitterUtils/test/ImpactPointEstimator_test
//
BOOST_AUTO_TEST_CASE(SingleTrackDistanceParametersAthenaRegression) {
  Estimator ipEstimator = makeEstimator(1.9971546939_T);
  Estimator::State state(magFieldCache());

  // Use same values as in Athena unit test
  Vector4 pos1(2_mm, 1_mm, -10_mm, 0_ns);
  Vector3 mom1(400_MeV, 600_MeV, 200_MeV);
  Vector3 vtxPos(1.2_mm, 0.8_mm, -7_mm);

  // Start creating some track parameters
  auto perigeeSurface =
      Surface::makeShared<PerigeeSurface>(pos1.segment<3>(ePos0));
  // Some fixed track parameter values
  auto params1 = BoundTrackParameters::create(
                     perigeeSurface, geoContext, pos1, mom1, mom1.norm(), 1_e,
                     BoundTrackParameters::CovarianceMatrix::Identity())
                     .value();

  // Compare w/ desired result from Athena unit test
  auto distance =
      ipEstimator.calculate3dDistance(geoContext, params1, vtxPos, state)
          .value();
  CHECK_CLOSE_ABS(distance, 3.10391_mm, 0.00001_mm);

  auto res2 = ipEstimator.estimate3DImpactParameters(
      geoContext, magFieldContext, params1, vtxPos, state);
  BOOST_CHECK(res2.ok());
  BoundTrackParameters endParams = *res2;
  Vector3 surfaceCenter = endParams.referenceSurface().center(geoContext);

  BOOST_CHECK_EQUAL(surfaceCenter, vtxPos);
}

// Test the Impact3d Point estimator 2d and 3d lifetimes sign
// on a single track.

BOOST_AUTO_TEST_CASE(Lifetimes2d3d) {
  Estimator ipEstimator = makeEstimator(2_T);

  // Create a track from a decay
  BoundVector trk_par;
  trk_par[eBoundLoc0] = 200_um;
  trk_par[eBoundLoc1] = 300_um;
  trk_par[eBoundTime] = 1_ns;
  trk_par[eBoundPhi] = 45_degree;
  trk_par[eBoundTheta] = 45_degree;
  trk_par[eBoundQOverP] = 1_e / 10_GeV;

  Vector4 ip_pos{0., 0., 0., 0.};
  Vertex<BoundTrackParameters> ip_vtx(ip_pos, makeVertexCovariance(), {});

  // Form the bound track parameters at the ip
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(ip_pos.head<3>());
  BoundTrackParameters track(perigeeSurface, trk_par,
                             makeBoundParametersCovariance());

  Vector3 direction{0., 1., 0.};
  auto lifetimes_signs = ipEstimator.getLifetimesSignOfTrack(
      track, ip_vtx, direction, geoContext, magFieldContext);

  // Check if the result is OK
  BOOST_CHECK(lifetimes_signs.ok());

  // Check that d0 sign is positive
  BOOST_CHECK((*lifetimes_signs).first > 0.);

  // Check that z0 sign is negative
  BOOST_CHECK((*lifetimes_signs).second < 0.);

  // Check the 3d sign

  auto sign3d = ipEstimator.get3DLifetimeSignOfTrack(
      track, ip_vtx, direction, geoContext, magFieldContext);

  // Check result is OK
  BOOST_CHECK(sign3d.ok());

  // Check 3D sign (should be positive)
  BOOST_CHECK((*sign3d) > 0.);
}

// Check `.estimateImpactParameters`.
BOOST_DATA_TEST_CASE(SingeTrackImpactParameters, tracks* vertices, d0, l0, t0,
                     phi, theta, p, q, vx0, vy0, vz0, vt0) {
  BoundVector par;
  par[eBoundLoc0] = d0;
  par[eBoundLoc1] = l0;
  par[eBoundTime] = t0;
  par[eBoundPhi] = phi;
  par[eBoundTheta] = theta;
  par[eBoundQOverP] = q / p;
  Vector4 vtxPos;
  vtxPos[ePos0] = vx0;
  vtxPos[ePos1] = vy0;
  vtxPos[ePos2] = vz0;
  vtxPos[eTime] = vt0;

  Estimator ipEstimator = makeEstimator(1_T);
  Estimator::State state(magFieldCache());

  // reference position and corresponding perigee surface
  Vector3 refPosition(0., 0., 0.);
  auto perigeeSurface = Surface::makeShared<PerigeeSurface>(refPosition);
  // create track and vertex
  BoundTrackParameters track(perigeeSurface, par,
                             makeBoundParametersCovariance());
  Vertex<BoundTrackParameters> myConstraint(vtxPos, makeVertexCovariance(), {});

  // check that computed impact parameters are meaningful
  ImpactParametersAndSigma output =
      ipEstimator
          .estimateImpactParameters(track, myConstraint, geoContext,
                                    magFieldContext)
          .value();
  BOOST_CHECK_NE(output.IPd0, 0.);
  BOOST_CHECK_NE(output.IPz0, 0.);
  // TODO what about the other struct members? can the parameter space be
  // restricted further?
}

BOOST_AUTO_TEST_SUITE_END()
