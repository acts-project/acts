// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef INTEGRATIONTEST_PROPAGATION_HELPER_H
#define INTEGRATIONTEST_PROPAGATION_HELPER_H

#include "ACTS/Propagator/Propagator.hpp"
#include "covariance_validation_fixture.hpp"

namespace tt = boost::test_tools;

namespace Acts {

using namespace propagation;
using units::Nat2SI;

namespace IntegrationTest {

  template <typename Propagator_type>
  void
  constant_field_propagation(const Propagator_type& propagator,
                             double                 pT,
                             double                 phi,
                             double                 theta,
                             double                 charge,
                             int                    index,
                             double                 Bz,
                             double                 disttol = 0.1 * units::_um)
  {

    // setup propagation options
    typename Propagator_type::template Options<> options;
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
    BOOST_TEST((exp_x - tp->position()(0)) == 0., tt::tolerance(disttol));
    BOOST_TEST((exp_y - tp->position()(1)) == 0., tt::tolerance(disttol));
    BOOST_TEST((exp_z - tp->position()(2)) == 0., tt::tolerance(disttol));
    // clang-format on
  }

  template <typename Propagator_type>
  void
  foward_backward(const Propagator_type& propagator,
                  double                 pT,
                  double                 phi,
                  double                 theta,
                  double                 charge,
                  int                    index,
                  double                 disttol = 0.1 * units::_um,
                  double                 momtol  = 1 * units::_keV)
  {

    // setup propagation options
    typename Propagator_type::template Options<> fwd_options;
    fwd_options.max_path_length = 5 * units::_m;
    fwd_options.max_step_size   = 1 * units::_cm;

    typename Propagator_type::template Options<> back_options;
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
    BOOST_TEST((x - tp2->position()(0)) == 0.,  tt::tolerance(disttol));
    BOOST_TEST((y - tp2->position()(1)) == 0.,  tt::tolerance(disttol));
    BOOST_TEST((z - tp2->position()(2)) == 0.,  tt::tolerance(disttol));
    BOOST_TEST((px - tp2->momentum()(0)) == 0., tt::tolerance(momtol));
    BOOST_TEST((py - tp2->momentum()(1)) == 0., tt::tolerance(momtol));
    BOOST_TEST((pz - tp2->momentum()(2)) == 0., tt::tolerance(momtol));
    // clang-format on
  }

  template <typename Propagator_type>
  void
  covaraiance_check_curvilinear(const Propagator_type& propagator,
                                double                 pT,
                                double                 phi,
                                double                 theta,
                                double                 charge,
                                int                    index,
                                double                 reltol = 2e-7)
  {

    covariance_validation_fixture<Propagator_type> fixture(propagator);
    // setup propagation options
    typename Propagator_type::template Options<> options;
    // setup propagation options
    options.max_step_size   = 1 * units::_cm;
    options.max_path_length = 5 * units::_m;

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

    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);

    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);

    // do propagation
    const auto  result = propagator.propagate(start, options);
    const auto& tp     = result.endParameters;

    // get numerically propagated covariance matrix

    ActsSymMatrixD<5> calculated_cov
        = fixture.calculateCovariance(start, *tp, options);

    BOOST_TEST(
        (calculated_cov - *tp->covariance()).norm()
                / std::min(calculated_cov.norm(), tp->covariance()->norm())
            == 0.,
        tt::tolerance(reltol));
  }
}
}

#endif
