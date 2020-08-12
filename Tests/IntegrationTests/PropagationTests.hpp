// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <utility>

/// Construct (initial) curvilinear parameters.
inline Acts::CurvilinearParameters makeParametersCurvilinear(double phi,
                                                             double theta,
                                                             double absMom,
                                                             double charge) {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  Vector3D pos = Vector3D::Zero();
  double time = 0.0;
  Vector3D mom = absMom * makeDirectionUnitFromPhiTheta(phi, theta);
  CurvilinearParameters params(std::nullopt, pos, mom, charge, time);

  // ensure initial parameters are valid
  CHECK_CLOSE_ABS(params.position(), pos, 0.125_um);
  CHECK_CLOSE_ABS(params.time(), time, 1_ps);
  CHECK_CLOSE_ABS(params.momentum(), mom, 0.125_eV);
  // charge should be identical not just similar
  BOOST_CHECK_EQUAL(params.charge(), charge);

  return params;
}

/// Check that two parameters object are consistent within the tolerances.
template <typename charge_t>
inline void checkParametersConsistency(
    const Acts::SingleBoundTrackParameters<charge_t>& cmp,
    const Acts::SingleBoundTrackParameters<charge_t>& ref,
    const Acts::GeometryContext& geoCtx, double epsPos, double epsDir,
    double epsMom) {
  using namespace Acts;

  // check stored parameters
  CHECK_CLOSE_ABS(cmp.template get<eBoundLoc0>(),
                  ref.template get<eBoundLoc0>(), epsPos);
  CHECK_CLOSE_ABS(cmp.template get<eBoundLoc1>(),
                  ref.template get<eBoundLoc1>(), epsPos);
  CHECK_CLOSE_ABS(cmp.template get<eBoundTime>(),
                  ref.template get<eBoundTime>(), epsPos);
  // check phi equivalence with circularity
  CHECK_SMALL(detail::radian_sym(cmp.template get<eBoundPhi>() -
                                 ref.template get<eBoundPhi>()),
              epsDir);
  CHECK_CLOSE_ABS(cmp.template get<eBoundTheta>(),
                  ref.template get<eBoundTheta>(), epsDir);
  CHECK_CLOSE_ABS(cmp.template get<eBoundQOverP>(),
                  ref.template get<eBoundQOverP>(), epsMom);
  // check derived parameters
  CHECK_CLOSE_ABS(cmp.position(geoCtx), ref.position(geoCtx), epsPos);
  CHECK_CLOSE_ABS(cmp.time(), ref.time(), epsPos);
  CHECK_CLOSE_ABS(cmp.momentum(), ref.momentum(), epsMom);
  // charge should be identical not just similar
  BOOST_CHECK_EQUAL(cmp.charge(), ref.charge());
}

/// Propagate the initial parameters the given path length along its
/// trajectory and then propagate the final parameters back. Verify that the
/// propagated parameters match the initial ones.
template <typename propagator_t, template <typename, typename>
                                 class options_t = Acts::PropagatorOptions>
inline void runForwardBackwardTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::CurvilinearParameters& initialParams, double pathLength,
    double epsPos, double epsDir, double epsMom, bool showDebug) {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  using DebugOutput = DebugOutputActor;
  using Actions = ActionList<DebugOutput>;
  using Aborts = AbortList<>;

  // forward propagation
  options_t<Actions, Aborts> fwdOptions(geoCtx, magCtx);
  fwdOptions.direction = Acts::forward;
  fwdOptions.pathLimit = pathLength;
  fwdOptions.maxStepSize = 1_cm;
  fwdOptions.debug = showDebug;
  // backward propagation
  options_t<Actions, Aborts> bwdOptions(geoCtx, magCtx);
  bwdOptions.direction = Acts::backward;
  bwdOptions.pathLimit = -pathLength;
  bwdOptions.maxStepSize = 1_cm;
  bwdOptions.debug = showDebug;

  // propagate parameters forward
  auto fwdResult = propagator.propagate(initialParams, fwdOptions);
  BOOST_CHECK(fwdResult.ok());
  // propagate propagated parameters back
  auto bwdResult =
      propagator.propagate(*(fwdResult.value().endParameters), bwdOptions);
  BOOST_CHECK(bwdResult.ok());

  // check that initial and back-propagated parameters match
  checkParametersConsistency(initialParams, *(bwdResult.value().endParameters),
                             geoCtx, epsPos, epsDir, epsMom);

  if (showDebug) {
    auto fwdOutput = fwdResult.value().template get<DebugOutput::result_type>();
    auto fwdParams = *(fwdResult.value().endParameters);
    std::cout << ">>>>> Output for forward propagation " << std::endl;
    std::cout << fwdOutput.debugString << std::endl;
    std::cout << fwdParams << std::endl;
    auto bwdOutput = bwdResult.value().template get<DebugOutput::result_type>();
    auto bwdParams = *(bwdResult.value().endParameters);
    std::cout << ">>>>> Output for backward propagation " << std::endl;
    std::cout << bwdOutput.debugString << std::endl;
    std::cout << bwdParams << std::endl;
  }
}

/// Build a cylinder along z with the given radius.
std::shared_ptr<Acts::CylinderSurface> makeTargetCylinder(double radius) {
  using namespace Acts;

  auto transform = std::make_shared<Transform3D>(Transform3D::Identity());
  return Surface::makeShared<CylinderSurface>(
      std::move(transform), radius, std::numeric_limits<double>::max());
}

/// Propagate the initial parameters freely in space.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline std::pair<Acts::BoundParameters, double> transportFreely(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleBoundTrackParameters<charge_t>& initialParams,
    double pathLength, bool showDebug) {
  using namespace Acts::UnitLiterals;

  using DebugOutput = Acts::DebugOutputActor;
  using Actions = Acts::ActionList<DebugOutput>;
  using Aborts = Acts::AbortList<>;

  // setup propagation options
  options_t<Actions, Aborts> options(geoCtx, magCtx);
  options.direction = Acts::forward;
  options.pathLimit = pathLength;
  options.maxStepSize = 1_cm;
  options.debug = showDebug;

  auto result = propagator.propagate(initialParams, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  if (showDebug) {
    auto output = result.value().template get<DebugOutput::result_type>();
    auto params = *(result.value().endParameters);
    std::cout << ">>>>> Output for to-surface propagation " << std::endl;
    std::cout << output.debugString << std::endl;
    std::cout << params << std::endl;
  }

  return {*result.value().endParameters, result.value().pathLength};
}

/// Propagate the initial parameters to the target surface.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline std::pair<Acts::BoundParameters, double> transportToSurface(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleBoundTrackParameters<charge_t>& initialParams,
    const Acts::Surface& targetSurface, double pathLimit, bool showDebug) {
  using namespace Acts::UnitLiterals;

  using DebugOutput = Acts::DebugOutputActor;
  using Actions = Acts::ActionList<DebugOutput>;
  using Aborts = Acts::AbortList<>;

  // setup propagation options
  options_t<Actions, Aborts> options(geoCtx, magCtx);
  options.direction = Acts::forward;
  options.pathLimit = pathLimit;
  options.maxStepSize = 1_cm;
  options.debug = showDebug;

  auto result = propagator.propagate(initialParams, targetSurface, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  if (showDebug) {
    auto output = result.value().template get<DebugOutput::result_type>();
    auto params = *(result.value().endParameters);
    std::cout << ">>>>> Output for to-surface propagation " << std::endl;
    std::cout << output.debugString << std::endl;
    std::cout << params << std::endl;
  }

  return {*result.value().endParameters, result.value().pathLength};
}

// Propagate the initial parameters once for the given path length and
// use the propagated parameters to define a target surface. Propagate the
// initial parameters again to the target surface. Verify that the surface has
// been found and the parameters are consistent.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline void runToSurfaceTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleBoundTrackParameters<charge_t>& initialParams,
    double pathLength, double epsPos, double epsDir, double epsMom,
    bool showDebug) {
  // free propagation for the given path length
  auto [freeParams, freePathLength] = transportFreely(
      propagator, geoCtx, magCtx, initialParams, pathLength, showDebug);
  CHECK_CLOSE_ABS(freePathLength, pathLength, epsPos);

  // TODO produce surface
  auto surface = makeTargetCylinder(freeParams.position(geoCtx).norm());
  BOOST_CHECK(surface);

  // bound propagation onto the surface
  // increase path length limit to ensure the surface is reached
  auto [surfParams, surfPathLength] =
      transportToSurface(propagator, geoCtx, magCtx, initialParams, *surface,
                         2 * pathLength, showDebug);
  CHECK_CLOSE_ABS(surfPathLength, pathLength, epsPos);

  CHECK_CLOSE_ABS(freeParams.position(geoCtx), surfParams.position(geoCtx),
                  epsPos);
  CHECK_CLOSE_ABS(freeParams.time(), surfParams.time(), epsPos);
  CHECK_CLOSE_ABS(freeParams.momentum().normalized(),
                  surfParams.momentum().normalized(), epsDir);
  CHECK_CLOSE_ABS(freeParams.momentum().norm(), surfParams.momentum().norm(),
                  epsMom);
}
