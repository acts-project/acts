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
template <template <typename, typename> class options_t, typename propagator_t>
inline void runForwardBackwardTest(const propagator_t& propagator,
                                   const Acts::GeometryContext& geoCtx,
                                   const Acts::MagneticFieldContext& magCtx,
                                   const Acts::CurvilinearParameters& initial,
                                   double pathLength, double epsPos,
                                   double epsDir, double epsMom,
                                   bool debug = false) {
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
  fwdOptions.debug = debug;
  // backward propagation
  options_t<Actions, Aborts> bwdOptions(geoCtx, magCtx);
  bwdOptions.direction = Acts::backward;
  bwdOptions.pathLimit = -pathLength;
  bwdOptions.maxStepSize = 1_cm;
  bwdOptions.debug = debug;

  // propagate parameters forward
  auto fwdResult = propagator.propagate(initial, fwdOptions);
  BOOST_CHECK(fwdResult.ok());
  // propagate propagated parameters back
  auto bwdResult =
      propagator.propagate(*(fwdResult.value().endParameters), bwdOptions);
  BOOST_CHECK(bwdResult.ok());

  // check that initial and propagated parameters match
  checkParametersConsistency(initial, *(bwdResult.value().endParameters),
                             geoCtx, epsPos, epsDir, epsMom);

  if (debug) {
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

/// Propagate the initial parameters to the target surface.
template <typename propagator_t>
inline std::pair<Acts::BoundParameters, double> transportToSurface(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::CurvilinearParameters& initial,
    const Acts::Surface& targetSurface, double pathLimit, bool debug = false) {
  using namespace Acts::UnitLiterals;

  using DebugOutput = Acts::DebugOutputActor;
  using Actions = Acts::ActionList<DebugOutput>;
  using Aborts = Acts::AbortList<>;

  // setup propagation options
  Acts::PropagatorOptions<Actions, Aborts> options(geoCtx, magCtx);
  options.direction = Acts::forward;
  options.pathLimit = pathLimit;
  options.maxStepSize = 1_cm;
  options.debug = debug;

  auto result = propagator.propagate(initial, targetSurface, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  if (debug) {
    auto output = result.value().template get<DebugOutput::result_type>();
    auto params = *(result.value().endParameters);
    std::cout << ">>>>> Output for to-surface propagation " << std::endl;
    std::cout << output.debugString << std::endl;
    std::cout << params << std::endl;
  }

  return {*result.value().endParameters, result.value().pathLength};
}
