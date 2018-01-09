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
#include "ACTS/ACTSVersion.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/AtlasStepper.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "covariance_validation_fixture.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using units::Nat2SI;
using namespace propagation;

namespace IntegrationTest {

  /// test forward propagation in constant magnetic field
  BOOST_DATA_TEST_CASE(
      constant_bfield_forward_propagation,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 2 * M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_int_distribution<>(-1, 1)))
          ^ bdata::xrange(1000),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    typedef ConstantBField            BField_type;
    typedef EigenStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    // setup propagator with constant B-field
    const double    Bz = 2 * units::_T;
    BField_type     bField(0, 0, Bz);
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
    BOOST_TEST((pT - tp->momentum().perp()) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((pz - tp->momentum()(2)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((theta - tp->momentum().theta()) == 0., tt::tolerance(1e-4));
    // clang-format on

    // calculate bending radius
    double r = std::abs(Nat2SI<units::MOMENTUM>(pT) / (q * Bz));
    // calculate number of turns of helix
    double turns = options.max_path_length / (2 * M_PI * r) * sin(theta);
    // respect direction of curl
    turns = (q * Bz < 0) ? turns : -turns;

    // calculate expected final momentum direction in phi in [-pi,pi]
    double exp_phi = std::fmod(phi + turns * 2 * M_PI, 2 * M_PI);
    if (exp_phi < -M_PI) exp_phi += 2 * M_PI;
    if (exp_phi > M_PI) exp_phi -= 2 * M_PI;

    // calculate expected position
    double exp_z = z + pz / pT * 2 * M_PI * r * std::abs(turns);

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
    BOOST_TEST((exp_phi - tp->momentum().phi()) == 0., tt::tolerance(1e-4));
    BOOST_TEST((exp_x - tp->position()(0)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((exp_y - tp->position()(1)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((exp_z - tp->position()(2)) == 0., tt::tolerance(0.1 * units::_um));
    // clang-format on
  }

  /// test consistency of forward-backward propagation in constant field
  BOOST_DATA_TEST_CASE(
      forward_backward_propagation,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 2 * M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_int_distribution<>(-1, 1)))
          ^ bdata::xrange(1000),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    typedef ConstantBField            BField_type;
    typedef EigenStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    // setup propagator with constant B-field
    const double    Bz = 2 * units::_T;
    BField_type     bField(0, 0, Bz);
    Stepper_type    stepper(std::move(bField));
    Propagator_type propagator(std::move(stepper));

    // setup propagation options
    Propagator_type::Options<> fwd_options;
    fwd_options.max_path_length = 5 * units::_m;
    fwd_options.max_step_size   = 1 * units::_cm;

    Propagator_type::Options<> back_options;
    back_options.direction       = backward;
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
    const auto& tp2 = propagator.propagate(*tp1, back_options).endParameters;

    // test propagation invariants
    // clang-format off
    BOOST_TEST((x - tp2->position()(0)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((y - tp2->position()(1)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((z - tp2->position()(2)) == 0., tt::tolerance(0.1 * units::_um));
    BOOST_TEST((px - tp2->momentum()(0)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((py - tp2->momentum()(1)) == 0., tt::tolerance(1 * units::_keV));
    BOOST_TEST((pz - tp2->momentum()(2)) == 0., tt::tolerance(1 * units::_keV));
    // clang-format on
  }

  /// test correct covariance transport for curvilinear parameters
  BOOST_DATA_TEST_CASE(
      covariance_transport_curvilinear,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(2. * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 2 * M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 0,
                           bdata::distribution
                           = std::uniform_int_distribution<>(-1, 1)))
          ^ bdata::xrange(1000),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    typedef ConstantBField            BField_type;
    typedef EigenStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    // setup propagator with constant B-field
    const double    Bz = 2 * units::_T;
    BField_type     bField(0, 0, Bz);
    Stepper_type    stepper(std::move(bField));
    Propagator_type propagator(std::move(stepper));

    covariance_validation_fixture<Propagator_type> fixture(propagator);

    // setup propagation options
    Propagator_type::Options<> options;
    options.max_path_length = 5 * units::_m;
    options.max_step_size   = 1 * units::_cm;

    // define start parameters
    double            x  = 0;
    double            y  = 0;
    double            z  = 0;
    double            px = pT * cos(phi);
    double            py = pT * sin(phi);
    double            pz = pT / tan(theta);
    double            q  = (charge != 0) ? charge : +1;
    Vector3D          pos(x, y, z);
    Vector3D          mom(px, py, pz);
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);

    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);

    // do propagation
    const auto& tp = propagator.propagate(start, options).endParameters;

    // get numerically propagated covariance matrix
    ActsSymMatrixD<5> calculated_cov
        = fixture.calculateCovariance(start, options);

    if ((calculated_cov - *tp->covariance()).norm()
            / std::min(calculated_cov.norm(), tp->covariance()->norm())
        > 2e-7) {
      std::cout << "initial parameters = " << tp->parameters() << std::endl;
      std::cout << "calculated = " << calculated_cov << std::endl << std::endl;
      std::cout << "obtained = " << *tp->covariance() << std::endl;
    }

    BOOST_TEST(
        (calculated_cov - *tp->covariance()).norm()
                / std::min(calculated_cov.norm(), tp->covariance()->norm())
            == 0.,
        tt::tolerance(2e-7));
  }
}  // namespace Test

}  // namespace Acts
