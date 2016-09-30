// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE Propagation Tests
#include <boost/test/included/unit_test.hpp>
#include <cmath>

#include <boost/test/data/test_case.hpp>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AtlasStepper.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using units::Nat2SI;

namespace Test {

  BOOST_DATA_TEST_CASE(constant_bfield_forward_propagation,
                       bdata::random(0.4 * units::_GeV, 10. * units::_GeV)
                           ^ bdata::random(0., 2 * M_PI)
                           ^ bdata::random(0., M_PI)
                           ^ bdata::random(-1, 1)
                           ^ bdata::xrange(1000),
                       pT,
                       phi,
                       theta,
                       charge,
                       index)
  {
    typedef ConstantFieldSvc          BField_type;
    typedef AtlasStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    // setup propagator with constant B-field
    const double        Bz = 2 * units::_T;
    BField_type::Config c;
    c.field = {0, 0, Bz};
    BField_type     bField(std::move(c));
    Stepper_type    stepper(std::move(bField));
    Propagator_type propagator(std::move(stepper));

    // setup propagation options
    Propagator_type::Options<> options;
    options.max_path_length = 5 * units::_m;
    options.max_step_size   = 1 * units::_cm;

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = (charge != 0) ? charge : +1;
    Vector3D              pos(x, y, z);
    Vector3D              mom(px, py, pz);
    CurvilinearParameters pars(nullptr, pos, mom, q);

    // do propagation
    const auto& tp = propagator.propagate(pars, options).endParameters;

    // test propagation invariants
    // clang-format off
    BOOST_TEST((pT - tp.momentum().perp()) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((pz - tp.momentum()(2)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((theta - tp.momentum().theta()) == 0., tt::tolerance(1e-4));
    // clang-format on

    // calculate bending radius
    double r = fabs(Nat2SI<units::MOMENTUM>(pT) / (q * Bz));
    // calculate number of turns of helix
    double turns = options.max_path_length / (2 * M_PI * r) * sin(theta);
    // respect direction of curl
    turns = (q * Bz < 0) ? turns : -turns;

    // calculate expected final momentum direction in phi in [-pi,pi]
    double exp_phi = std::fmod(phi + turns * 2 * M_PI, 2 * M_PI);
    if (exp_phi < -M_PI) exp_phi += 2 * M_PI;
    if (exp_phi > M_PI) exp_phi -= 2 * M_PI;

    // calculate expected position
    double exp_z = z + pz / pT * 2 * M_PI * r * fabs(turns);

    // calculate center of bending circle in transverse plane
    double xc, yc;
    // offset with respect to starting point
    double dx = r * cos(M_PI / 2 - phi);
    double dy = r * sin(M_PI / 2 - phi);
    if (q * Bz < 0) {
      xc = x - dx;
      yc = y + dy;
    } else {
      xc = x + dx;
      yc = y - dy;
    }
    // phi position of starting point in bending circle
    double phi0 = std::atan2(y - yc, x - xc);

    // calculated expected position in transverse plane
    double exp_x = xc + r * cos(phi0 + turns * 2 * M_PI);
    double exp_y = yc + r * sin(phi0 + turns * 2 * M_PI);

    // clang-format off
    BOOST_TEST((exp_phi - tp.momentum().phi()) == 0., tt::tolerance(1e-4));
    BOOST_TEST((exp_x - tp.position()(0)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((exp_y - tp.position()(1)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((exp_z - tp.position()(2)) == 0., tt::tolerance(0.1 * units::_um));
    // clang-format on
  }

  BOOST_DATA_TEST_CASE(forward_backward_propagation,
                       bdata::random(0.4 * units::_GeV, 10. * units::_GeV)
                           ^ bdata::random(0., 2 * M_PI)
                           ^ bdata::random(0., M_PI)
                           ^ bdata::random(-1, 1)
                           ^ bdata::xrange(1000),
                       pT,
                       phi,
                       theta,
                       charge,
                       index)
  {
    typedef ConstantFieldSvc          BField_type;
    typedef AtlasStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    // setup propagator with constant B-field
    const double        Bz = 2 * units::_T;
    BField_type::Config c;
    c.field = {0, 0, Bz};
    BField_type     bField(std::move(c));
    Stepper_type    stepper(std::move(bField));
    Propagator_type propagator(std::move(stepper));

    // setup propagation options
    Propagator_type::Options<forward> fwd_options;
    fwd_options.max_path_length = 5 * units::_m;
    fwd_options.max_step_size   = 1 * units::_cm;

    Propagator_type::Options<backward> back_options;
    back_options.max_path_length = 5 * units::_m;
    back_options.max_step_size   = 1 * units::_cm;

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = (charge != 0) ? charge : +1;
    Vector3D              pos(x, y, z);
    Vector3D              mom(px, py, pz);
    CurvilinearParameters start(nullptr, pos, mom, q);

    // do forward-backward propagation
    const auto& tp1 = propagator.propagate(start, fwd_options).endParameters;
    const auto& tp2 = propagator.propagate(tp1, back_options).endParameters;

    // test propagation invariants
    // clang-format off
    BOOST_TEST((x - tp2.position()(0)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((y - tp2.position()(1)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((z - tp2.position()(2)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((px - tp2.momentum()(0)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((py - tp2.momentum()(1)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((pz - tp2.momentum()(2)) == 0., tt::tolerance(1 * units::_keV));
    // clang-format on
  }
}  // namespace Test

}  // namespace Acts
