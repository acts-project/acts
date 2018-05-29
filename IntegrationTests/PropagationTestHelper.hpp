// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "covariance_validation_fixture.hpp"

namespace tt = boost::test_tools;

namespace Acts {

using namespace propagation;
using units::Nat2SI;

namespace IntegrationTest {

  /// Helper method to create a transform for a plane
  /// to mimic detector situations, the plane is roughly
  /// perpendicular to the track
  ///
  /// @param nnomal The nominal normal direction
  /// @param angleT Rotation around the norminal normal
  /// @param angleU Roation around the original U axis
  std::shared_ptr<Transform3D>
  createPlanarTransform(const Vector3D& nposition,
                        const Vector3D& nnormal,
                        double          angleT,
                        double          angleU)
  {
    // the rotation of the destination surface
    Vector3D T = nnormal.normalized();
    Vector3D U = std::abs(T.dot(Vector3D::UnitZ())) < 0.99
        ? Vector3D::UnitZ().cross(T).normalized()
        : Vector3D::UnitX().cross(T).normalized();
    Vector3D V = T.cross(U);
    // that's the plane curvilinear Rotation
    RotationMatrix3D curvilinearRotation;
    curvilinearRotation.col(0) = U;
    curvilinearRotation.col(1) = V;
    curvilinearRotation.col(2) = T;
    // curvilinear surfaces are boundless
    Transform3D ctransform{curvilinearRotation};
    ctransform.pretranslate(nposition);
    ctransform.prerotate(AngleAxis3D(angleT, T));
    ctransform.prerotate(AngleAxis3D(angleU, U));
    //
    return std::make_shared<Transform3D>(ctransform);
  }

  /// Helper method to create a transform for a plane
  /// to mimic detector situations, the plane is roughly
  /// perpendicular to the track
  ///
  /// @param nnomal The nominal normal direction
  /// @param angleT Rotation around the norminal normal
  /// @param angleU Roation around the original U axis
  std::shared_ptr<Transform3D>
  createCylindricTransform(const Vector3D& nposition,
                           double          angleX,
                           double          angleY)
  {
    Transform3D ctransform;
    ctransform.setIdentity();
    ctransform.pretranslate(nposition);
    ctransform.prerotate(AngleAxis3D(angleX, Vector3D::UnitX()));
    ctransform.prerotate(AngleAxis3D(angleY, Vector3D::UnitY()));
    return std::make_shared<Transform3D>(ctransform);
  }

  template <typename Propagator_type>
  Vector3D
  constant_field_propagation(const Propagator_type& propagator,
                             double                 pT,
                             double                 phi,
                             double                 theta,
                             double                 charge,
                             int /*index*/,
                             double Bz,
                             double disttol = 0.1 * units::_um)
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
    double                q  = charge;
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
    return tp->position();
  }

  template <typename Propagator_type>
  void
  foward_backward(const Propagator_type& propagator,
                  double                 pT,
                  double                 phi,
                  double                 theta,
                  double                 charge,
                  double                 plimit,
                  int /*index*/,
                  double disttol = 1. * units::_um,
                  double momtol  = 10. * units::_keV)
  {

    // setup propagation options
    typename Propagator_type::template Options<> fwd_options;
    fwd_options.max_path_length = plimit;
    fwd_options.max_step_size   = 1 * units::_cm;

    typename Propagator_type::template Options<> back_options;
    back_options.direction       = backward;
    back_options.max_path_length = plimit;
    back_options.max_step_size   = 1 * units::_cm;

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = charge;
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

  // test propagation to cylinder
  template <typename Propagator_type>
  std::pair<Vector3D, double>
  to_cylinder(const Propagator_type& propagator,
              double                 pT,
              double                 phi,
              double                 theta,
              double                 charge,
              double                 plimit,
              double                 rand1,
              double                 rand2,
              double /*rand3*/,
              bool covtransport = false)
  {
    // setup propagation options
    typename Propagator_type::template Options<> options;
    // setup propagation options
    options.max_step_size   = plimit;
    options.max_path_length = plimit;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = charge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);

    std::unique_ptr<const ActsSymMatrixD<5>> cov_ptr = nullptr;
    if (covtransport) {
      ActsSymMatrixD<5> cov;
      // take some major correlations (off-diagonals)
      cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
          0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
          1. / (10 * units::_GeV);
      cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    }
    // do propagation of the start parameters
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);

    // The transform at the destination
    auto seTransform = createCylindricTransform(
        Vector3D(0., 0., 0.), 0.05 * rand1, 0.05 * rand2);
    CylinderSurface endSurface(
        seTransform, plimit * units::_m, std::numeric_limits<double>::max());

    // Increase the path limit - to be safe hitting the surface
    options.max_path_length *= 2;
    const auto  result = propagator.propagate(start, endSurface, options);
    const auto& tp     = result.endParameters;
    // check for null pointer
    BOOST_CHECK(tp != nullptr);
    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), result.pathLength);
  }

  // test propagation to most surfaces
  template <typename Propagator_type, typename Surface_type>
  std::pair<Vector3D, double>
  to_surface(const Propagator_type& propagator,
             double                 pT,
             double                 phi,
             double                 theta,
             double                 charge,
             double                 plimit,
             double                 rand1,
             double                 rand2,
             double                 rand3,
             bool                   planar       = true,
             bool                   covtransport = false)
  {
    // setup propagation options
    typename Propagator_type::template Options<> options;
    // setup propagation options
    options.max_step_size   = plimit;
    options.max_path_length = plimit;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = charge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);

    std::unique_ptr<const ActsSymMatrixD<5>> cov_ptr = nullptr;
    if (covtransport) {
      ActsSymMatrixD<5> cov;
      // take some major correlations (off-diagonals)
      cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
          0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
          1. / (10 * units::_GeV);
      cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    }
    // Create curvilinear start parameters
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);
    const auto            result_s = propagator.propagate(start, options);
    const auto&           tp_s     = result_s.endParameters;

    // The transform at the destination
    auto seTransform = planar ? createPlanarTransform(tp_s->position(),
                                                      tp_s->momentum().unit(),
                                                      0.1 * rand3,
                                                      0.1 * rand1)
                              : createCylindricTransform(tp_s->position(),
                                                         0.05 * rand1,
                                                         0.05 * rand2);

    Surface_type endSurface(seTransform, nullptr);
    // Increase the path limit - to be safe hitting the surface
    options.max_path_length *= 2;
    const auto  result = propagator.propagate(start, endSurface, options);
    const auto& tp     = result.endParameters;
    // check the result for nullptr
    BOOST_CHECK(tp != nullptr);
    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), result.pathLength);
  }

  template <typename Propagator_type>
  void
  covariance_curvilinear(const Propagator_type& propagator,
                         double                 pT,
                         double                 phi,
                         double                 theta,
                         double                 charge,
                         double                 plimit,
                         int /*index*/,
                         double reltol = 1e-4)
  {
    covariance_validation_fixture<Propagator_type> fixture(propagator);
    // setup propagation options
    typename Propagator_type::template Options<> options;
    // setup propagation options
    options.max_step_size   = plimit;
    options.max_path_length = plimit;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = charge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);

    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // do propagation of the start parameters
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);
    CurvilinearParameters start_wo_c(nullptr, pos, mom, q);

    const auto  result = propagator.propagate(start, options);
    const auto& tp     = result.endParameters;

    // get numerically propagated covariance matrix
    ActsSymMatrixD<5> calculated_cov = fixture.calculateCovariance(
        start_wo_c, *(start.covariance()), *tp, options);
    ActsSymMatrixD<5> obtained_cov = (*(tp->covariance()));
    bool cov_similar = calculated_cov.isApprox(obtained_cov, reltol);
    BOOST_CHECK(cov_similar);
    if (!cov_similar) {
      BOOST_CHECK_EQUAL(calculated_cov, obtained_cov);
    }
  }

  template <typename Propagator_type,
            typename StartSurface_type,
            typename DestSurface_type>
  void
  covariance_bound(const Propagator_type& propagator,
                   double                 pT,
                   double                 phi,
                   double                 theta,
                   double                 charge,
                   double                 plimit,
                   double                 rand1,
                   double                 rand2,
                   double                 rand3,
                   int /*index*/,
                   bool   startPlanar = true,
                   bool   destPlanar  = true,
                   double reltol      = 1e-3)
  {
    covariance_validation_fixture<Propagator_type> fixture(propagator);
    // setup propagation options
    typename Propagator_type::template Options<> options;
    // setup propagation options
    options.max_step_size   = plimit;
    options.max_path_length = plimit;

    // define start parameters
    double            x  = 0;
    double            y  = 0;
    double            z  = 0;
    double            px = pT * cos(phi);
    double            py = pT * sin(phi);
    double            pz = pT / tan(theta);
    double            q  = charge;
    Vector3D          pos(x, y, z);
    Vector3D          mom(px, py, pz);
    ActsSymMatrixD<5> cov;

    // take some major correlations (off-diagonals)
    // cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162,
    // 0,
    //     0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
    //     1. / (10 * units::_GeV);

    cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 0.1,
        0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);

    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // create curvilinear start parameters
    CurvilinearParameters start_c(nullptr, pos, mom, q);
    const auto            result_c = propagator.propagate(start_c, options);
    const auto&           tp_c     = result_c.endParameters;

    auto ssTransform = startPlanar
        ? createPlanarTransform(pos, mom.unit(), 0.1 * rand1, 0.1 * rand2)
        : createCylindricTransform(pos, 0.05 * rand1, 0.05 * rand2);
    auto seTransform = destPlanar
        ? createPlanarTransform(tp_c->position(),
                                tp_c->momentum().unit(),
                                0.1 * rand3,
                                0.1 * rand1)
        : createCylindricTransform(
              tp_c->position(), 0.05 * rand1, 0.05 * rand2);

    StartSurface_type startSurface(ssTransform, nullptr);
    BoundParameters   start(std::move(cov_ptr), pos, mom, q, startSurface);
    BoundParameters   start_wo_c(nullptr, pos, mom, q, startSurface);

    // increase the path limit - to be safe hitting the surface
    options.max_path_length *= 2;

    DestSurface_type endSurface(seTransform, nullptr);
    const auto       result = propagator.propagate(start, endSurface, options);
    const auto&      tp     = result.endParameters;

    // get obtained covariance matrix
    ActsSymMatrixD<5> obtained_cov = (*(tp->covariance()));

    // get numerically propagated covariance matrix
    ActsSymMatrixD<5> calculated_cov = fixture.calculateCovariance(
        start_wo_c, *(start.covariance()), *tp, options);

    bool cov_similar = calculated_cov.isApprox(obtained_cov, reltol);
    BOOST_CHECK(cov_similar);
    if (!cov_similar) {
      // the second check is for screen output
      BOOST_CHECK_EQUAL(calculated_cov, obtained_cov);
    }
  }
}
}