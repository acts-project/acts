// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace tt = boost::test_tools;

namespace Acts {
namespace IntegrationTest {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

/// Helper method to create a transform for a plane
/// to mimic detector situations, the plane is roughly
/// perpendicular to the track
///
/// @param nnomal The nominal normal direction
/// @param angleT Rotation around the norminal normal
/// @param angleU Roation around the original U axis
std::shared_ptr<Transform3D> createPlanarTransform(const Vector3D& nposition,
                                                   const Vector3D& nnormal,
                                                   double angleT,
                                                   double angleU) {
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
std::shared_ptr<Transform3D> createCylindricTransform(const Vector3D& nposition,
                                                      double angleX,
                                                      double angleY) {
  Transform3D ctransform;
  ctransform.setIdentity();
  ctransform.pretranslate(nposition);
  ctransform.prerotate(AngleAxis3D(angleX, Vector3D::UnitX()));
  ctransform.prerotate(AngleAxis3D(angleY, Vector3D::UnitY()));
  return std::make_shared<Transform3D>(ctransform);
}

template <typename Propagator_type>
Vector3D constant_field_propagation(const Propagator_type& propagator,
                                    double pT, double phi, double theta,
                                    double charge, double time, double Bz,
                                    double disttol = 0.1 *
                                                     Acts::UnitConstants::um,
                                    bool debug = false) {
  using namespace Acts::UnitLiterals;
  namespace VH = VectorHelpers;

  // setup propagation options
  PropagatorOptions<> options(tgContext, mfContext);
  options.pathLimit = 5_m;
  options.maxStepSize = 1_cm;
  options.debug = debug;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  CurvilinearParameters pars(std::nullopt, pos, mom, q, time);

  // do propagation
  const auto& tp = propagator.propagate(pars, options).value().endParameters;

  // test propagation invariants
  // clang-format off
    CHECK_CLOSE_ABS(pT, VH::perp(tp->momentum()), 1_keV);
    CHECK_CLOSE_ABS(pz, tp->momentum()(2), 1_keV);
    CHECK_CLOSE_ABS(theta, VH::theta(tp->momentum()), 1e-4);
  // clang-format on

  double r = (q * Bz != 0.) ? std::abs(pT / (q * Bz))
                            : std::numeric_limits<double>::max();

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
void foward_backward(const Propagator_type& propagator, double pT, double phi,
                     double theta, double charge, double plimit,
                     double disttol = 1 * Acts::UnitConstants::um,
                     double momtol = 10 * Acts::UnitConstants::keV,
                     bool debug = false) {
  using namespace Acts::UnitLiterals;

  // setup propagation options
  // Action list and abort list
  using DebugOutput = Acts::DebugOutputActor;
  using ActionList = Acts::ActionList<DebugOutput>;

  PropagatorOptions<ActionList> fwdOptions(tgContext, mfContext);
  fwdOptions.pathLimit = plimit;
  fwdOptions.maxStepSize = 1_cm;
  fwdOptions.debug = debug;

  PropagatorOptions<ActionList> bwdOptions(tgContext, mfContext);
  bwdOptions.direction = backward;
  bwdOptions.pathLimit = -plimit;
  bwdOptions.maxStepSize = 1_cm;
  bwdOptions.debug = debug;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  CurvilinearParameters start(std::nullopt, pos, mom, q, time);

  // do forward-backward propagation
  const auto& fwdResult = propagator.propagate(start, fwdOptions).value();
  const auto& bwdResult =
      propagator.propagate(*fwdResult.endParameters, bwdOptions).value();

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

    auto bwdOutput = bwdResult.template get<DebugOutput::result_type>();
    std::cout << ">>>>> Output for backward propagation " << std::endl;
    std::cout << bwdOutput.debugString << std::endl;
    std::cout << " - resulted at position : "
              << bwdResult.endParameters->position() << std::endl;
  }
}

// test propagation to cylinder
template <typename Propagator_type>
std::pair<Vector3D, double> to_cylinder(
    const Propagator_type& propagator, double pT, double phi, double theta,
    double charge, double plimit, double rand1, double rand2, double /*rand3*/,
    bool covtransport = false, bool debug = false) {
  using namespace Acts::UnitLiterals;

  // setup propagation options
  PropagatorOptions<> options(tgContext, mfContext);
  // setup propagation options
  options.maxStepSize = plimit * 0.1;
  options.pathLimit = plimit;
  options.debug = debug;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtransport) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }

  // The transform at the destination
  auto seTransform = createCylindricTransform(Vector3D(0., 0., 0.),
                                              0.04 * rand1, 0.04 * rand2);
  auto endSurface = Surface::makeShared<CylinderSurface>(
      seTransform, plimit, std::numeric_limits<double>::max());

  // Increase the path limit - to be safe hitting the surface
  options.pathLimit *= 2;

  if (q == 0.) {
    auto start = new NeutralCurvilinearTrackParameters(covOpt, pos, mom, time);

    const auto result =
        propagator.propagate(*start, *endSurface, options).value();
    const auto& tp = result.endParameters;
    // check for null pointer
    BOOST_CHECK(tp != nullptr);
    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), result.pathLength);
  } else {
    auto start = new CurvilinearParameters(covOpt, pos, mom, q, time);

    const auto result =
        propagator.propagate(*start, *endSurface, options).value();
    const auto& tp = result.endParameters;
    // check for null pointer
    BOOST_CHECK(tp != nullptr);
    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), result.pathLength);
  }
}

// test propagation to most surfaces
template <typename Propagator_type, typename SurfaceType>
std::pair<Vector3D, double> to_surface(
    const Propagator_type& propagator, double pT, double phi, double theta,
    double charge, double plimit, double rand1, double rand2, double rand3,
    bool planar = true, bool covtransport = false, bool debug = false) {
  using namespace Acts::UnitLiterals;
  using DebugOutput = DebugOutputActor;

  // setup propagation options
  PropagatorOptions<ActionList<DebugOutput>> options(tgContext, mfContext);
  // setup propagation options
  options.maxStepSize = plimit;
  options.pathLimit = plimit;
  options.debug = debug;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtransport) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }
  // Create curvilinear start parameters
  if (q == 0.) {
    auto start = new NeutralCurvilinearTrackParameters(covOpt, pos, mom, time);
    const auto result_s = propagator.propagate(*start, options).value();
    const auto& tp_s = result_s.endParameters;
    // The transform at the destination
    auto seTransform =
        planar ? createPlanarTransform(tp_s->position(),
                                       tp_s->momentum().normalized(),
                                       0.1 * rand3, 0.1 * rand1)
               : createCylindricTransform(tp_s->position(), 0.04 * rand1,
                                          0.04 * rand2);

    auto endSurface = Surface::makeShared<SurfaceType>(seTransform, nullptr);
    // Increase the path limit - to be safe hitting the surface
    options.pathLimit *= 2;

    if (debug) {
      std::cout << ">>> Path limit for this propgation is set to: "
                << options.pathLimit << std::endl;
    }

    auto result = propagator.propagate(*start, *endSurface, options);
    const auto& propRes = *result;
    const auto& tp = propRes.endParameters;
    // check the result for nullptr
    BOOST_CHECK(tp != nullptr);

    // screen output in case you are running in debug mode
    if (debug) {
      const auto& debugOutput =
          propRes.template get<DebugOutput::result_type>();
      std::cout << ">>> Debug output of this propagation " << std::endl;
      std::cout << debugOutput.debugString << std::endl;
      std::cout << ">>> Propagation status is : ";
      if (result.ok()) {
        std::cout << "success";
      } else {
        std::cout << result.error();
      }
      std::cout << std::endl;
    }

    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), propRes.pathLength);
  } else {
    auto start = new CurvilinearParameters(covOpt, pos, mom, q, time);
    const auto result_s = propagator.propagate(*start, options).value();
    const auto& tp_s = result_s.endParameters;
    // The transform at the destination
    auto seTransform =
        planar ? createPlanarTransform(tp_s->position(),
                                       tp_s->momentum().normalized(),
                                       0.1 * rand3, 0.1 * rand1)
               : createCylindricTransform(tp_s->position(), 0.04 * rand1,
                                          0.04 * rand2);

    auto endSurface = Surface::makeShared<SurfaceType>(seTransform, nullptr);
    // Increase the path limit - to be safe hitting the surface
    options.pathLimit *= 2;

    if (debug) {
      std::cout << ">>> Path limit for this propgation is set to: "
                << options.pathLimit << std::endl;
    }

    auto result = propagator.propagate(*start, *endSurface, options);
    const auto& propRes = *result;
    const auto& tp = propRes.endParameters;
    // check the result for nullptr
    BOOST_CHECK(tp != nullptr);

    // screen output in case you are running in debug mode
    if (debug) {
      const auto& debugOutput =
          propRes.template get<DebugOutput::result_type>();
      std::cout << ">>> Debug output of this propagation " << std::endl;
      std::cout << debugOutput.debugString << std::endl;
      std::cout << ">>> Propagation status is : ";
      if (result.ok()) {
        std::cout << "success";
      } else {
        std::cout << result.error();
      }
      std::cout << std::endl;
    }
    // The position and path length
    return std::pair<Vector3D, double>(tp->position(), propRes.pathLength);
  }
}

template <typename Propagator_type>
Covariance covariance_curvilinear(const Propagator_type& propagator, double pT,
                                  double phi, double theta, double charge,
                                  double plimit, bool debug = false) {
  using namespace Acts::UnitLiterals;

  // setup propagation options
  DenseStepperPropagatorOptions<> options(tgContext, mfContext);
  options.maxStepSize = plimit;
  options.pathLimit = 0.1 * plimit;
  options.debug = debug;
  options.tolerance = 1e-9;

  // define start parameters
  double x = 1.;
  double y = 0.;
  double z = 0.;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  Covariance cov;
  // take some major correlations (off-diagonals)
  // clang-format off
  cov <<
   10_mm, 0, 0.123, 0, 0.5, 0,
   0, 10_mm, 0, 0.162, 0, 0,
   0.123, 0, 0.1, 0, 0, 0,
   0, 0.162, 0, 0.1, 0, 0,
   0.5, 0, 0, 0, 1_e / 10_GeV, 0,
   0, 0, 0, 0, 0, 1_us;
  // clang-format on

  // do propagation of the start parameters
  CurvilinearParameters start(cov, pos, mom, q, time);

  const auto result = propagator.propagate(start, options).value();
  const auto& tp = result.endParameters;

  return *(tp->covariance());
}

template <typename Propagator_type, typename StartSurfaceType,
          typename DestSurfaceType>
Covariance covariance_bound(const Propagator_type& propagator, double pT,
                            double phi, double theta, double charge,
                            double plimit, double rand1, double rand2,
                            double rand3, bool startPlanar = true,
                            bool destPlanar = true, bool debug = false) {
  using namespace Acts::UnitLiterals;

  // setup propagation options
  DenseStepperPropagatorOptions<> options(tgContext, mfContext);
  options.maxStepSize = 0.1 * plimit;
  options.pathLimit = plimit;
  options.debug = debug;

  // define start parameters
  double x = 1.;
  double y = 0.;
  double z = 0.;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  Covariance cov;

  // take some major correlations (off-diagonals)
  // clang-format off
  cov <<
   10_mm, 0, 0.123, 0, 0.5, 0,
   0, 10_mm, 0, 0.162, 0, 0,
   0.123, 0, 0.1, 0, 0, 0,
   0, 0.162, 0, 0.1, 0, 0,
   0.5, 0, 0, 0, 1_e / 10_GeV, 0,
   0, 0, 0, 0, 0, 1_us;
  // clang-format on

  // create curvilinear start parameters
  CurvilinearParameters start_c(std::nullopt, pos, mom, q, time);
  const auto result_c = propagator.propagate(start_c, options).value();
  const auto& tp_c = result_c.endParameters;

  auto ssTransform =
      startPlanar ? createPlanarTransform(pos, mom.normalized(), 0.05 * rand1,
                                          0.05 * rand2)
                  : createCylindricTransform(pos, 0.01 * rand1, 0.01 * rand2);
  auto seTransform = destPlanar
                         ? createPlanarTransform(tp_c->position(),
                                                 tp_c->momentum().normalized(),
                                                 0.05 * rand3, 0.05 * rand1)
                         : createCylindricTransform(tp_c->position(),
                                                    0.01 * rand1, 0.01 * rand2);

  auto startSurface =
      Surface::makeShared<StartSurfaceType>(ssTransform, nullptr);
  BoundParameters start(tgContext, cov, pos, mom, q, time, startSurface);

  // increase the path limit - to be safe hitting the surface
  options.pathLimit *= 2;

  auto endSurface = Surface::makeShared<DestSurfaceType>(seTransform, nullptr);
  const auto result = propagator.propagate(start, *endSurface, options).value();
  const auto& tp = result.endParameters;

  // get obtained covariance matrix
  return *(tp->covariance());
}
}  // namespace IntegrationTest
}  // namespace Acts
