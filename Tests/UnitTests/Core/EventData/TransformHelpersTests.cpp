// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

#include "TrackParametersDatasets.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();
}

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(CanonicalDirectionUnderAlignment) {
  class AlignedElement : public DetectorElementStub {
   public:
    using DetectorElementStub::DetectorElementStub;
    const Transform3& localToGlobalTransform(
        const GeometryContext& gctx) const override {
      return gctx.get<Transform3>();
    }
  };
  AlignedElement element(Transform3::Identity(),
                         std::make_shared<RectangleBounds>(10., 10.), 0.1);
  const auto& surface = element.surface();
  Transform3 aligned = Transform3::Identity();
  aligned.translate(Vector3(4., 5., 6.));
  aligned.rotate(AngleAxis3(0.7, Vector3::UnitY()));
  const BoundVector bound(1., 2., 0.3, 1.2, 0.01, 4.);
  const Vector3 direction = makeDirectionFromPhiTheta(0.3, 1.2);
  for (const Transform3& placement : {Transform3::Identity(), aligned}) {
    const GeometryContext gctx{placement};
    BOOST_CHECK(surface.boundParametrization(gctx).isGlobalSpherical());
    const auto free = transformBoundToFreeParameters(surface, gctx, bound);
    CHECK_CLOSE_ABS(free.segment<3>(eFreePos0), placement * Vector3(1., 2., 0.),
                    1e-14);
    CHECK_CLOSE_ABS(free.segment<3>(eFreeDir0), direction, 1e-14);
    const auto restored = transformFreeToBoundParameters(free, surface, gctx);
    BOOST_REQUIRE(restored.ok());
    CHECK_CLOSE_ABS(*restored, bound, 1e-14);
    const auto restoredComponents = transformFreeToBoundParameters(
        free.segment<3>(eFreePos0), free[eFreeTime], direction,
        free[eFreeQOverP], surface, gctx);
    BOOST_REQUIRE(restoredComponents.ok());
    CHECK_CLOSE_ABS(*restoredComponents, bound, 1e-14);
  }
}

BOOST_DATA_TEST_CASE(Parameters,
                     surfaces * posSymmetric * posSymmetric * ts * phis *
                         thetas * ps * qsNonZero,
                     surface, l0, l1, time, phi, theta, p, q) {
  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();

  Vector2 loc(l0, l1);
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);
  // transform reference position
  Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);

  const auto qOverP = q / p;

  // construct bound parameters
  BoundVector bv = BoundVector::Zero();
  bv[eBoundLoc0] = l0;
  bv[eBoundLoc1] = l1;
  bv[eBoundTime] = time;
  bv[eBoundPhi] = phi;
  bv[eBoundTheta] = theta;
  bv[eBoundQOverP] = qOverP;

  // convert to free parameters
  FreeVector fv = transformBoundToFreeParameters(*surface, geoCtx, bv);

  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreePos0), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeTime], bv[eBoundTime], eps, eps);
  CHECK_CLOSE_REL(fv.segment<3>(eFreeDir0).norm(), 1, eps);
  CHECK_CLOSE_OR_SMALL(fv.segment<3>(eFreeDir0), dir, eps, eps);
  CHECK_CLOSE_OR_SMALL(fv[eFreeQOverP], bv[eBoundQOverP], eps, eps);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TransformFreeToBound)

BOOST_DATA_TEST_CASE(GlobalToBoundTrackParameters,
                     surfaces * posSymmetric * posSymmetric * ts * phis *
                         thetas * ps * qsNonZero,
                     surface, l0, l1, time, phiInput, theta, p, q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) && (theta < std::numbers::pi)) ? phiInput : 0.;
  const auto qOverP = q / p;

  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Vector2 loc(l0, l1);
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);
  // transform reference position
  Vector3 pos = surface->localToGlobal(geoCtx, loc, dir);

  // convert free parameters to bound parameters
  {
    BOOST_TEST_INFO("Transform free parameters vector onto surface "
                    << surface->name());

    FreeVector fv = FreeVector::Zero();
    fv[eFreePos0] = pos[ePos0];
    fv[eFreePos1] = pos[ePos1];
    fv[eFreePos2] = pos[ePos2];
    fv[eFreeTime] = time;
    fv[eFreeDir0] = dir[eMom0];
    fv[eFreeDir1] = dir[eMom1];
    fv[eFreeDir2] = dir[eMom2];
    fv[eFreeQOverP] = qOverP;
    BoundVector bv =
        transformFreeToBoundParameters(fv, *surface, geoCtx).value();
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }

  // Assert failure when trying to convert a position that is not on-surface.
  {
    Vector3 posOff = pos + surface->normal(geoCtx, loc) * 0.5;
    BOOST_TEST_INFO("Transform free parameters vector onto surface "
                    << surface->name());

    FreeVector fv = FreeVector::Zero();
    fv[eFreePos0] = posOff[ePos0];
    fv[eFreePos1] = posOff[ePos1];
    fv[eFreePos2] = posOff[ePos2];
    fv[eFreeTime] = time;
    fv[eFreeDir0] = dir[eMom0];
    fv[eFreeDir1] = dir[eMom1];
    fv[eFreeDir2] = dir[eMom2];
    fv[eFreeQOverP] = qOverP;
    auto res = transformFreeToBoundParameters(fv, *surface, geoCtx);
    BOOST_CHECK(!res.ok());
  }

  // convert separate components to bound parameters
  {
    BOOST_TEST_INFO("Transform free parameters components onto surface "
                    << surface->name());

    BoundVector bv =
        transformFreeToBoundParameters(pos, time, dir, qOverP, *surface, geoCtx)
            .value();
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc0], l0, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundLoc1], l1, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }

  // Assert failure when trying to convert a position that is not on-surface.
  {
    BOOST_TEST_INFO("Transform free parameters components onto surface "
                    << surface->name());

    Vector3 posOff = pos + surface->normal(geoCtx, loc) * 0.5;
    auto res = transformFreeToBoundParameters(posOff, time, dir, qOverP,
                                              *surface, geoCtx);
    BOOST_CHECK(!res.ok());
  }
}

BOOST_DATA_TEST_CASE(GlobalToCurvilinearParameters,
                     ts * phis * thetas * ps * qsNonZero, time, phiInput, theta,
                     p, q) {
  // phi is ill-defined in forward/backward tracks
  const auto phi = ((0 < theta) && (theta < std::numbers::pi)) ? phiInput : 0.;
  const auto qOverP = q / p;

  auto geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  Vector3 dir = makeDirectionFromPhiTheta(phi, theta);

  // convert w/ direction
  {
    BoundVector bv = transformFreeToCurvilinearParameters(time, dir, qOverP);
    CHECK_SMALL(bv[eBoundLoc0], eps);
    CHECK_SMALL(bv[eBoundLoc1], eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
  // convert w/ angles
  {
    BoundVector bv =
        transformFreeToCurvilinearParameters(time, phi, theta, qOverP);
    CHECK_SMALL(bv[eBoundLoc0], eps);
    CHECK_SMALL(bv[eBoundLoc1], eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTime], time, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundPhi], phi, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundTheta], theta, eps, eps);
    CHECK_CLOSE_OR_SMALL(bv[eBoundQOverP], qOverP, eps, eps);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
