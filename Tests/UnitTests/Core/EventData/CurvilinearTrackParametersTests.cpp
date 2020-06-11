// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ParametersTestHelper.hpp"

namespace Acts {
using VectorHelpers::phi;
using VectorHelpers::theta;
namespace Test {

/// @brief Unit test for Curvilinear parameters
///
BOOST_AUTO_TEST_CASE(curvilinear_initialization) {
  using namespace Acts::UnitLiterals;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // some position and momentum
  Vector3D pos(1.34_mm, 2.34_mm, 3.45_mm);
  Vector3D mom(1000_GeV, 1000_GeV, -0.100_GeV);
  Vector3D dir(mom.normalized());
  Vector3D z_axis_global(0., 0., 1.);
  /// create curvilinear parameters without covariance +1/-1 charge
  CurvilinearParameters curvilinear_pos(std::nullopt, pos, mom, 1_e, 1_s);
  CurvilinearParameters curvilinear_neg(std::nullopt, pos, mom, -1_e, 2.5_s);
  NeutralCurvilinearTrackParameters curvilinear_neut(std::nullopt, pos, mom,
                                                     33.33_s);

  /// check local coordinates
  const auto fphi = phi(mom);
  const auto ftheta = theta(mom);
  const double oOp = 1. / mom.norm();

  // check parameters
  consistencyCheck(curvilinear_pos, pos, mom, +1_e, 1_s,
                   {0., 0., fphi, ftheta, oOp, 1_s});
  consistencyCheck(curvilinear_neg, pos, mom, -1_e, 2.5_s,
                   {0., 0., fphi, ftheta, -oOp, 2.5_s});
  consistencyCheck(curvilinear_neut, pos, mom, 0., 33.33_s,
                   {0., 0., fphi, ftheta, oOp, 33.33_s});

  // check that the created surface is at the position
  CHECK_CLOSE_REL(curvilinear_pos.referenceSurface().center(tgContext), pos,
                  1e-6);
  CHECK_CLOSE_REL(curvilinear_neg.referenceSurface().center(tgContext), pos,
                  1e-6);
  CHECK_CLOSE_REL(curvilinear_neut.referenceSurface().center(tgContext), pos,
                  1e-6);

  // check that the z-axis of the created surface is along momentum direction
  CHECK_CLOSE_REL(curvilinear_pos.referenceSurface().normal(tgContext, pos),
                  dir, 1e-6);
  CHECK_CLOSE_REL(curvilinear_neg.referenceSurface().normal(tgContext, pos),
                  dir, 1e-6);
  CHECK_CLOSE_REL(curvilinear_neut.referenceSurface().normal(tgContext, pos),
                  dir, 1e-6);

  // check the reference frame of curvilinear parameters
  // it is the x-y frame of the created surface
  RotationMatrix3D mFrame = RotationMatrix3D::Zero();
  Vector3D tAxis = dir;
  Vector3D uAxis = (z_axis_global.cross(tAxis)).normalized();
  Vector3D vAxis = tAxis.cross(uAxis);
  mFrame.col(0) = uAxis;
  mFrame.col(1) = vAxis;
  mFrame.col(2) = tAxis;
  CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_pos.referenceFrame(tgContext), 1e-6,
                       1e-9);
  CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_neg.referenceFrame(tgContext), 1e-6,
                       1e-9);
  CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_neut.referenceFrame(tgContext), 1e-6,
                       1e-9);

  /// copy construction test
  CurvilinearParameters curvilinear_pos_copy(curvilinear_pos);
  CurvilinearParameters curvilinear_neg_copy(curvilinear_neg);
  NeutralCurvilinearTrackParameters curvilinear_neut_copy(curvilinear_neut);

  BOOST_CHECK_EQUAL(curvilinear_pos_copy, curvilinear_pos);
  BOOST_CHECK_EQUAL(curvilinear_neg_copy, curvilinear_neg);
  BOOST_CHECK_EQUAL(curvilinear_neut_copy, curvilinear_neut);

  /// modification test with set methods
  double ux = 0.1_mm;
  double uy = 0.5_mm;
  curvilinear_pos_copy.set<eLOC_0>(tgContext, ux);
  curvilinear_pos_copy.set<eLOC_1>(tgContext, uy);
  // the local parameter should still be (0,0) for Curvilinear
  BOOST_CHECK_EQUAL(curvilinear_pos_copy.parameters()[eLOC_0], 0u);
  BOOST_CHECK_EQUAL(curvilinear_pos_copy.parameters()[eLOC_1], 0u);
  // the position should be updated though
  Vector3D uposition =
      curvilinear_neg_copy.referenceSurface().transform(tgContext) *
      Vector3D(ux, uy, 0.);
  // the position should be updated
  CHECK_CLOSE_REL(curvilinear_pos_copy.position(), uposition, 1e-6);
  // it should be the position of the surface
  CHECK_CLOSE_REL(curvilinear_pos_copy.referenceSurface().center(tgContext),
                  uposition, 1e-6);

  double uphi = 1.2;
  double utheta = 0.2;
  double uqop = 0.025_e / 1_GeV;
  double ut = 1337_s;

  curvilinear_pos_copy.set<ePHI>(tgContext, uphi);
  curvilinear_pos_copy.set<eTHETA>(tgContext, utheta);
  curvilinear_pos_copy.set<eQOP>(tgContext, uqop);
  curvilinear_pos_copy.set<eT>(tgContext, ut);
  // we should have a new updated momentum
  Vector3D umomentum =
      40. *
      Vector3D(cos(uphi) * sin(utheta), sin(uphi) * sin(utheta), cos(utheta)) *
      1_GeV;
  CHECK_CLOSE_REL(umomentum, curvilinear_pos_copy.momentum(), 1e-6);
  // the updated momentum should be the col(2) of the transform
  CHECK_CLOSE_REL(umomentum.normalized(),
                  curvilinear_pos_copy.referenceSurface()
                      .transform(tgContext)
                      .rotation()
                      .col(2),
                  1e-6);
}
}  // namespace Test
}  // namespace Acts
