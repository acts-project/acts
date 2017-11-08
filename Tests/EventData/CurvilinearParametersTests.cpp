// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE CurvilinearParameters Tests
#include <boost/test/included/unit_test.hpp>
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ParametersTestHelper.hpp"

namespace Acts {
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
    const auto   fphi   = mom.phi();
    const auto   ftheta = mom.theta();
    const double oOp    = 1. / mom.mag();

    // check parameters
    consistencyCheck(
        curvilinear_pos, pos, mom, +1., {{0., 0., fphi, ftheta, oOp}});
    consistencyCheck(
        curvilinear_neg, pos, mom, -1., {{0., 0., fphi, ftheta, -oOp}});
    consistencyCheck(
        curvilinear_neut, pos, mom, 0., {{0., 0., fphi, ftheta, oOp}});

    // check that the created surface is at the position
    BOOST_CHECK_EQUAL(curvilinear_pos.referenceSurface().center(), pos);
    BOOST_CHECK_EQUAL(curvilinear_neg.referenceSurface().center(), pos);
    BOOST_CHECK_EQUAL(curvilinear_neut.referenceSurface().center(), pos);

    // check that the z-axis of the created surface is along momentum direction
    BOOST_CHECK_EQUAL(curvilinear_pos.referenceSurface().normal(pos), dir);
    BOOST_CHECK_EQUAL(curvilinear_neg.referenceSurface().normal(pos), dir);
    BOOST_CHECK_EQUAL(curvilinear_neut.referenceSurface().normal(pos), dir);

    // check the reference frame of curvilinear parameters
    // it is the x-y frame of the created surface
    RotationMatrix3D mFrame = RotationMatrix3D::Zero();
    Vector3D         tAxis  = dir;
    Vector3D         uAxis  = (z_axis_global.cross(tAxis)).normalized();
    Vector3D         vAxis  = tAxis.cross(uAxis);
    mFrame.col(0)           = uAxis;
    mFrame.col(1)           = vAxis;
    mFrame.col(2)           = tAxis;
    BOOST_CHECK_EQUAL(mFrame, curvilinear_pos.referenceFrame());
    BOOST_CHECK_EQUAL(mFrame, curvilinear_neg.referenceFrame());
    BOOST_CHECK_EQUAL(mFrame, curvilinear_neut.referenceFrame());

    /// copy construction test
    CurvilinearParameters        curvilinear_pos_copy(curvilinear_pos);
    CurvilinearParameters        curvilinear_neg_copy(curvilinear_neg);
    NeutralCurvilinearParameters curvilinear_neut_copy(curvilinear_neut);

    BOOST_CHECK_EQUAL(curvilinear_pos_copy, curvilinear_pos);
    BOOST_CHECK_EQUAL(curvilinear_neg_copy, curvilinear_neg);
    BOOST_CHECK_EQUAL(curvilinear_neut_copy, curvilinear_neut);

    /// modification test with set methods
    double uphi   = 1.2;
    double utheta = 0.2;
    double uqop   = 0.025;

    curvilinear_pos_copy.set<Acts::ePHI>(uphi);
    curvilinear_pos_copy.set<Acts::eTHETA>(utheta);
    curvilinear_pos_copy.set<Acts::eQOP>(uqop);
    // we should have a new updated momentum
    Vector3D umomentum = 40. * Vector3D(cos(uphi) * sin(utheta),
                                        sin(uphi) * sin(utheta),
                                        cos(utheta));

    BOOST_CHECK(umomentum.isApprox(curvilinear_pos_copy.momentum()));
  }
}  // end of namespace Test
}  // end of namespace Acts
