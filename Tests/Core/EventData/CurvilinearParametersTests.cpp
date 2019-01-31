// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE CurvilinearParameters Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "ParametersTestHelper.hpp"

namespace Acts {
using VectorHelpers::phi;
using VectorHelpers::theta;
namespace Test {

  /// @brief Unit test for Curvilinear parameters
  ///
  BOOST_AUTO_TEST_CASE(curvilinear_initialization)
  {

    // some position and momentum
    Vector3D pos(1.34, 2.34, 3.45);
    Vector3D mom(1000., 1000., -0.100);
    Vector3D dir(mom.normalized());
    Vector3D z_axis_global(0., 0., 1.);
    /// create curvilinear parameters without covariance +1/-1 charge
    CurvilinearParameters        curvilinear_pos(nullptr, pos, mom, 1.);
    CurvilinearParameters        curvilinear_neg(nullptr, pos, mom, -1.);
    NeutralCurvilinearParameters curvilinear_neut(nullptr, pos, mom);

    /// check local coordinates
    const auto   fphi   = phi(mom);
    const auto   ftheta = theta(mom);
    const double oOp    = 1. / mom.norm();

    // check parameters
    consistencyCheck(
        curvilinear_pos, pos, mom, +1., {{0., 0., fphi, ftheta, oOp}});
    consistencyCheck(
        curvilinear_neg, pos, mom, -1., {{0., 0., fphi, ftheta, -oOp}});
    consistencyCheck(
        curvilinear_neut, pos, mom, 0., {{0., 0., fphi, ftheta, oOp}});

    // check that the created surface is at the position
    CHECK_CLOSE_REL(curvilinear_pos.referenceSurface().center(), pos, 1e-6);
    CHECK_CLOSE_REL(curvilinear_neg.referenceSurface().center(), pos, 1e-6);
    CHECK_CLOSE_REL(curvilinear_neut.referenceSurface().center(), pos, 1e-6);

    // check that the z-axis of the created surface is along momentum direction
    CHECK_CLOSE_REL(curvilinear_pos.referenceSurface().normal(pos), dir, 1e-6);
    CHECK_CLOSE_REL(curvilinear_neg.referenceSurface().normal(pos), dir, 1e-6);
    CHECK_CLOSE_REL(curvilinear_neut.referenceSurface().normal(pos), dir, 1e-6);

    // check the reference frame of curvilinear parameters
    // it is the x-y frame of the created surface
    RotationMatrix3D mFrame = RotationMatrix3D::Zero();
    Vector3D         tAxis  = dir;
    Vector3D         uAxis  = (z_axis_global.cross(tAxis)).normalized();
    Vector3D         vAxis  = tAxis.cross(uAxis);
    mFrame.col(0)           = uAxis;
    mFrame.col(1)           = vAxis;
    mFrame.col(2)           = tAxis;
    CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_pos.referenceFrame(), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_neg.referenceFrame(), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(mFrame, curvilinear_neut.referenceFrame(), 1e-6, 1e-9);

    /// copy construction test
    CurvilinearParameters        curvilinear_pos_copy(curvilinear_pos);
    CurvilinearParameters        curvilinear_neg_copy(curvilinear_neg);
    NeutralCurvilinearParameters curvilinear_neut_copy(curvilinear_neut);

    BOOST_CHECK_EQUAL(curvilinear_pos_copy, curvilinear_pos);
    BOOST_CHECK_EQUAL(curvilinear_neg_copy, curvilinear_neg);
    BOOST_CHECK_EQUAL(curvilinear_neut_copy, curvilinear_neut);

    /// modification test with set methods
    double ux = 0.1;
    double uy = 0.5;
    curvilinear_pos_copy.set<eLOC_0>(ux);
    curvilinear_pos_copy.set<eLOC_1>(uy);
    // the local parameter should still be (0,0) for Curvilinear
    BOOST_CHECK_EQUAL(curvilinear_pos_copy.parameters()[eLOC_0], 0);
    BOOST_CHECK_EQUAL(curvilinear_pos_copy.parameters()[eLOC_1], 0);
    // the position should be updated though
    Vector3D uposition = curvilinear_neg_copy.referenceSurface().transform()
        * Vector3D(ux, uy, 0.);
    // the position should be updated
    CHECK_CLOSE_REL(curvilinear_pos_copy.position(), uposition, 1e-6);
    // it should be the position of the surface
    CHECK_CLOSE_REL(
        curvilinear_pos_copy.referenceSurface().center(), uposition, 1e-6);

    double uphi   = 1.2;
    double utheta = 0.2;
    double uqop   = 0.025;

    curvilinear_pos_copy.set<ePHI>(uphi);
    curvilinear_pos_copy.set<eTHETA>(utheta);
    curvilinear_pos_copy.set<eQOP>(uqop);
    // we should have a new updated momentum
    Vector3D umomentum = 40. * Vector3D(cos(uphi) * sin(utheta),
                                        sin(uphi) * sin(utheta),
                                        cos(utheta));
    CHECK_CLOSE_REL(umomentum, curvilinear_pos_copy.momentum(), 1e-6);
    // the updated momentum should be the col(2) of the transform
    CHECK_CLOSE_REL(
        umomentum.normalized(),
        curvilinear_pos_copy.referenceSurface().transform().rotation().col(2),
        1e-6);
  }
}  // namespace Test
}  // namespace Acts
