// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include "covariance_validation_fixture.hpp"

namespace tt = boost::test_tools;

namespace Acts {

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
                             double disttol = 0.1 * units::_um,
                             bool   debug   = false)
  {

    namespace VH = VectorHelpers;

    // setup propagation options
    PropagatorOptions<> options;
    options.pathLimit   = 5 * units::_m;
    options.maxStepSize = 1 * units::_cm;
    options.debug       = debug;

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
    CHECK_CLOSE_ABS(pT, VH::perp(tp->momentum()), 1 * units::_keV);
    CHECK_CLOSE_ABS(pz, tp->momentum()(2), 1 * units::_keV);
    CHECK_CLOSE_ABS(theta, VH::theta(tp->momentum()), 1e-4);
    // clang-format on

    double r = std::abs(Nat2SI<units::MOMENTUM>(pT) / (q * Bz));

    // calculate number of turns of helix
    double turns = options.pathLimit / (2 * M_PI * r) * sin(theta);
    // respect direction of curl
    turns = (q * Bz < 0) ? turns : -turns;

    // calculate expected final momentum direction in phi in [-pi,pi]
    double exp_phi = std::fmod(phi + turns * 2 * M_PI, 2 * M_PI);
    if (exp_phi < -M_PI) {
      exp_phi += 2 * M_PI;
    }
    if (exp_phi > M_PI) {
      exp_phi -= 2 * M_PI;
    }

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
    CHECK_CLOSE_ABS(exp_phi, VH::phi(tp->momentum()), 1e-4);
    CHECK_CLOSE_ABS(exp_x, tp->position()(0), disttol);
    CHECK_CLOSE_ABS(exp_y, tp->position()(1), disttol);
    CHECK_CLOSE_ABS(exp_z, tp->position()(2), disttol);
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
                  double momtol  = 10. * units::_keV,
                  bool   debug   = false)
  {

    // setup propagation options
    // Action list and abort list
    using DebugOutput = Acts::detail::DebugOutputActor;
    using ActionList  = Acts::ActionList<DebugOutput>;

    PropagatorOptions<ActionList> fwdOptions;
    fwdOptions.pathLimit   = plimit;
    fwdOptions.maxStepSize = 1 * units::_cm;
    fwdOptions.debug       = debug;

    PropagatorOptions<ActionList> bwdOptions;
    bwdOptions.direction   = backward;
    bwdOptions.pathLimit   = -plimit;
    bwdOptions.maxStepSize = 1 * units::_cm;
    bwdOptions.debug       = debug;

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
    const auto& fwdResult = propagator.propagate(start, fwdOptions);
    const auto& bwdResult
        = propagator.propagate(*fwdResult.endParameters, bwdOptions);

    const Vector3D& bwdPosition = bwdResult.endParameters->position();
    const Vector3D& bwdMomentum = bwdResult.endParameters->momentum();

    // test propagation invariants
    // clang-format off
    CHECK_CLOSE_ABS(x, bwdPosition(0), disttol);
    CHECK_CLOSE_ABS(y, bwdPosition(1), disttol);
    CHECK_CLOSE_ABS(z, bwdPosition(2), disttol);
    CHECK_CLOSE_ABS(px, bwdMomentum(0), momtol);
    CHECK_CLOSE_ABS(py, bwdMomentum(1), momtol);
    CHECK_CLOSE_ABS(pz, bwdMomentum(2), momtol);
    // clang-format on

    if (debug) {
      auto fwdOutput = fwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>>>> Output for forward propagation " << std::endl;
      std::cout << fwdOutput.debugString << std::endl;
      std::cout << " - resulted at position : "
                << fwdResult.endParameters->position() << std::endl;

      auto bwdOutput = fwdResult.template get<DebugOutput::result_type>();
      std::cout << ">>>>> Output for backward propagation " << std::endl;
      std::cout << bwdOutput.debugString << std::endl;
      std::cout << " - resulted at position : "
                << bwdResult.endParameters->position() << std::endl;
    }
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
              bool covtransport = false,
              bool debug        = false)
  {
    // setup propagation options
    PropagatorOptions<> options;
    // setup propagation options
    options.maxStepSize = plimit;
    options.pathLimit   = plimit;
    options.debug       = debug;

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

    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    if (covtransport) {
      ActsSymMatrixD<5> cov;
      // take some major correlations (off-diagonals)
      cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
          0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
          1. / (10 * units::_GeV);
      covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    }
    // do propagation of the start parameters
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    // The transform at the destination
    auto seTransform = createCylindricTransform(
        Vector3D(0., 0., 0.), 0.05 * rand1, 0.05 * rand2);
    auto endSurface = Surface::makeShared<CylinderSurface>(
        seTransform, plimit * units::_m, std::numeric_limits<double>::max());

    // Increase the path limit - to be safe hitting the surface
    options.pathLimit *= 2;
    const auto  result = propagator.propagate(start, *endSurface, options);
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
             bool                   covtransport = false,
             bool                   debug        = false)
  {

    using DebugOutput = detail::DebugOutputActor;

    // setup propagation options
    PropagatorOptions<ActionList<DebugOutput>> options;
    // setup propagation options
    options.maxStepSize = plimit;
    options.pathLimit   = plimit;
    options.debug       = debug;

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

    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    if (covtransport) {
      ActsSymMatrixD<5> cov;
      // take some major correlations (off-diagonals)
      cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
          0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
          1. / (10 * units::_GeV);
      covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    }
    // Create curvilinear start parameters
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);
    const auto            result_s = propagator.propagate(start, options);
    const auto&           tp_s     = result_s.endParameters;

    // The transform at the destination
    auto seTransform = planar
        ? createPlanarTransform(tp_s->position(),
                                tp_s->momentum().normalized(),
                                0.1 * rand3,
                                0.1 * rand1)
        : createCylindricTransform(
              tp_s->position(), 0.05 * rand1, 0.05 * rand2);

    auto endSurface = Surface::makeShared<Surface_type>(seTransform, nullptr);
    // Increase the path limit - to be safe hitting the surface
    options.pathLimit *= 2;

    if (debug) {
      std::cout << ">>> Path limit for this propgation is set to: "
                << options.pathLimit << std::endl;
    }

    const auto  result = propagator.propagate(start, *endSurface, options);
    const auto& tp     = result.endParameters;
    // check the result for nullptr
    BOOST_CHECK(tp != nullptr);

    // screen output in case you are running in debug mode
    if (debug) {
      const auto& debugOutput = result.template get<DebugOutput::result_type>();
      std::cout << ">>> Debug output of this propagation " << std::endl;
      std::cout << debugOutput.debugString << std::endl;
      std::cout << ">>> Propagation status is : " << int(result.status)
                << std::endl;
    }

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
                         double reltol = 1e-3,
                         bool   debug  = false)
  {
    covariance_validation_fixture<Propagator_type> fixture(propagator);
    // setup propagation options
    DenseStepperPropagatorOptions<> options;
    // setup propagation options
    options.maxStepSize = plimit;
    options.pathLimit   = plimit;
    options.debug       = debug;
    options.tolerance   = 1e-7;

    // define start parameters
    double   x  = 1.;
    double   y  = 0.;
    double   z  = 0.;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = charge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);

    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10. * units::_mm, 0, 0.123, 0, 0.5, 0, 10. * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10. * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // do propagation of the start parameters
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);
    CurvilinearParameters start_wo_c(nullptr, pos, mom, q);

    const auto  result = propagator.propagate(start, options);
    const auto& tp     = result.endParameters;

    // get numerically propagated covariance matrix
    ActsSymMatrixD<5> calculated_cov = fixture.calculateCovariance(
        start_wo_c, *(start.covariance()), *tp, options);
    ActsSymMatrixD<5> obtained_cov = (*(tp->covariance()));
    CHECK_CLOSE_COVARIANCE(calculated_cov, obtained_cov, reltol);
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
                   double reltol      = 1e-3,
                   bool   debug       = false)
  {
    covariance_validation_fixture<Propagator_type> fixture(propagator);
    // setup propagation options
    DenseStepperPropagatorOptions<> options;
    options.maxStepSize = plimit;
    options.pathLimit   = plimit;
    options.debug       = debug;

    // define start parameters
    double            x  = 1.;
    double            y  = 0.;
    double            z  = 0.;
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

    cov << 10. * units::_mm, 0, 0, 0, 0, 0, 10. * units::_mm, 0, 0, 0, 0, 0,
        0.1, 0, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 1. / (10. * units::_GeV);

    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // create curvilinear start parameters
    CurvilinearParameters start_c(nullptr, pos, mom, q);
    const auto            result_c = propagator.propagate(start_c, options);
    const auto&           tp_c     = result_c.endParameters;

    auto ssTransform = startPlanar
        ? createPlanarTransform(pos, mom.normalized(), 0.1 * rand1, 0.1 * rand2)
        : createCylindricTransform(pos, 0.05 * rand1, 0.05 * rand2);
    auto seTransform = destPlanar
        ? createPlanarTransform(tp_c->position(),
                                tp_c->momentum().normalized(),
                                0.1 * rand3,
                                0.1 * rand1)
        : createCylindricTransform(
              tp_c->position(), 0.05 * rand1, 0.05 * rand2);

    auto startSurface
        = Surface::makeShared<StartSurface_type>(ssTransform, nullptr);
    BoundParameters start(std::move(covPtr), pos, mom, q, startSurface);
    BoundParameters start_wo_c(nullptr, pos, mom, q, startSurface);

    // increase the path limit - to be safe hitting the surface
    options.pathLimit *= 2;

    auto endSurface
        = Surface::makeShared<DestSurface_type>(seTransform, nullptr);
    const auto  result = propagator.propagate(start, *endSurface, options);
    const auto& tp     = result.endParameters;

    // get obtained covariance matrix
    ActsSymMatrixD<5> obtained_cov = (*(tp->covariance()));

    // get numerically propagated covariance matrix
    ActsSymMatrixD<5> calculated_cov = fixture.calculateCovariance(
        start_wo_c, *(start.covariance()), *tp, options);

    CHECK_CLOSE_COVARIANCE(calculated_cov, obtained_cov, reltol);
  }
}
}
