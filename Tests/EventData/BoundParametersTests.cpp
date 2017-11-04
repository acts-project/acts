// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE BoundParameters Tests
#include <boost/test/included/unit_test.hpp>
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ParametersTestHelper.hpp"

namespace Acts {
namespace Test {

  /// @brief Unit test for parameters at a plane
  ///
  BOOST_AUTO_TEST_CASE(bound_to_plane_test)
  {
    // create the plane first
    auto quaternion = Eigen::AngleAxisd{0.2343, Acts::Vector3D{0, 1, 0}};
    auto rotation   = quaternion.toRotationMatrix();
    auto transform  = std::make_shared<Transform3D>(quaternion);
    auto bounds     = std::make_shared<RectangleBounds>(100., 100.);
    PlaneSurface pSurface(transform, bounds);

    // now create parameters on this surface
    // l_x, l_y, phi, theta, q/p (1/p)
    std::array<double, 5> pars_array = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    const double phi   = pars_array[2];
    const double theta = pars_array[3];
    double       p     = fabs(1. / pars_array[4]);
    Vector3D     direction(
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    Vector3D mom = p * direction;
    Vector3D pos
        = pars_array[0] * rotation.col(0) + pars_array[1] * rotation.col(1);
    // constructor from parameter vector
    BoundParameters ataPlane_from_pars(nullptr, pars, pSurface);
    consistencyCheck(ataPlane_from_pars, pos, mom, 1., pars_array);
    // constructor from global parameters
    BoundParameters ataPlane_from_global(nullptr, pos, mom, 1., pSurface);
    consistencyCheck(ataPlane_from_global, pos, mom, 1., pars_array);
    // constructor for neutral parameters
    NeutralBoundParameters n_ataPlane_from_pars(nullptr, pars, pSurface);
    consistencyCheck(n_ataPlane_from_pars, pos, mom, 0., pars_array);
    // constructor for neutral global parameters
    NeutralBoundParameters n_ataPlane_from_global(nullptr, pars, pSurface);
    consistencyCheck(n_ataPlane_from_global, pos, mom, 0., pars_array);

    // check that indeed the surfaces are copied
    BOOST_CHECK(&(ataPlane_from_pars.referenceSurface())
                != &(ataPlane_from_global.referenceSurface()));

    // check that the reference frame is the rotation matrix
    BOOST_CHECK(ataPlane_from_pars.referenceFrame().isApprox(rotation));
  }

  /// @brief Unit test to disc
  ///
  BOOST_AUTO_TEST_CASE(bound_to_disc_test)
  {
    // create the plane first
    double z_position = 150.;
    auto   transform
        = std::make_shared<Transform3D>(Translation3D(0., 0., z_position));
    auto        bounds = std::make_shared<RadialBounds>(100., 1200.);
    DiscSurface dSurface(transform, bounds);

    // now create parameters on this surface
    // r, phi, phi, theta, q/p (1/p)
    std::array<double, 5> pars_array = {{125., 0.345, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    const double phi   = pars_array[2];
    const double theta = pars_array[3];
    double       p     = fabs(1. / pars_array[4]);
    Vector3D     direction(
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    Vector3D mom = p * direction;
    Vector3D pos(pars_array[0] * cos(pars_array[1]),
                 pars_array[0] * sin(pars_array[1]),
                 z_position);
    // constructor from parameter vector
    BoundParameters ataDisc_from_pars(nullptr, pars, dSurface);
    consistencyCheck(ataDisc_from_pars, pos, mom, 1., pars_array);
    // constructor from global parameters
    BoundParameters ataDisc_from_global(nullptr, pos, mom, 1., dSurface);
    consistencyCheck(ataDisc_from_global, pos, mom, 1., pars_array);
    // constructor for neutral parameters
    NeutralBoundParameters n_ataDisc_from_pars(nullptr, pars, dSurface);
    consistencyCheck(n_ataDisc_from_pars, pos, mom, 0., pars_array);
    // constructor for neutral global parameters
    NeutralBoundParameters n_ataDisc_from_global(nullptr, pars, dSurface);
    consistencyCheck(n_ataDisc_from_global, pos, mom, 0., pars_array);

    // check that indeed the surfaces are copied
    BOOST_CHECK(&(ataDisc_from_pars.referenceSurface())
                != &(ataDisc_from_global.referenceSurface()));
  }

  /// @brief Unit test to cylinder
  ///
  BOOST_AUTO_TEST_CASE(bound_to_cylinder_test) {}

  /// @brief Unit test to perigee
  ///
  BOOST_AUTO_TEST_CASE(bound_to_perigee_test) {}

  /// @brief Unit test to line
  ///
  BOOST_AUTO_TEST_CASE(bound_to_line_test) {}
}
}