// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <utility>

// parameter construction helpers

/// Construct (initial) curvilinear parameters.
inline Acts::CurvilinearTrackParameters makeParametersCurvilinear(
    double phi, double theta, double absMom, double charge) {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  // phi is ill-defined in forward/backward tracks. normalize the value to
  // ensure parameter comparisons give correct answers.
  if (not((0 < theta) and (theta < M_PI))) {
    phi = 0;
  }

  Vector4 pos4 = Vector4::Zero();
  return CurvilinearTrackParameters(pos4, phi, theta, absMom, charge);
}

/// Construct (initial) curvilinear parameters with covariance.
inline Acts::CurvilinearTrackParameters makeParametersCurvilinearWithCovariance(
    double phi, double theta, double absMom, double charge) {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  // phi is ill-defined in forward/backward tracks. normalize the value to
  // ensure parameter comparisons give correct answers.
  if (not((0 < theta) and (theta < M_PI))) {
    phi = 0;
  }

  BoundVector stddev = BoundVector::Zero();
  // TODO use momentum-dependent resolutions
  stddev[eBoundLoc0] = 15_um;
  stddev[eBoundLoc1] = 80_um;
  stddev[eBoundTime] = 25_ns;
  stddev[eBoundPhi] = 1_degree;
  stddev[eBoundTheta] = 1.5_degree;
  stddev[eBoundQOverP] = 1_e / 10_GeV;
  BoundSymMatrix corr = BoundSymMatrix::Identity();
  corr(eBoundLoc0, eBoundLoc1) = corr(eBoundLoc1, eBoundLoc0) = 0.125;
  corr(eBoundLoc0, eBoundPhi) = corr(eBoundPhi, eBoundLoc0) = 0.25;
  corr(eBoundLoc1, eBoundTheta) = corr(eBoundTheta, eBoundLoc1) = -0.25;
  corr(eBoundTime, eBoundQOverP) = corr(eBoundQOverP, eBoundTime) = 0.125;
  corr(eBoundPhi, eBoundTheta) = corr(eBoundTheta, eBoundPhi) = -0.25;
  corr(eBoundPhi, eBoundQOverP) = corr(eBoundPhi, eBoundQOverP) = -0.125;
  corr(eBoundTheta, eBoundQOverP) = corr(eBoundTheta, eBoundQOverP) = 0.5;
  BoundSymMatrix cov = stddev.asDiagonal() * corr * stddev.asDiagonal();

  Vector4 pos4 = Vector4::Zero();
  return CurvilinearTrackParameters(pos4, phi, theta, absMom, charge, cov);
}

/// Construct (initial) neutral curvilinear parameters.
inline Acts::NeutralCurvilinearTrackParameters makeParametersCurvilinearNeutral(
    double phi, double theta, double absMom) {
  using namespace Acts;
  using namespace Acts::UnitLiterals;

  // phi is ill-defined in forward/backward tracks. normalize the value to
  // ensure parameter comparisons give correct answers.
  if (not((0 < theta) and (theta < M_PI))) {
    phi = 0;
  }

  Vector4 pos4 = Vector4::Zero();
  return NeutralCurvilinearTrackParameters(pos4, phi, theta, 1 / absMom);
}

// helpers to compare track parameters

/// Check that two parameters object are consistent within the tolerances.
///
/// \warning Does not check that they are defined on the same surface.
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
  CHECK_CLOSE_ABS(cmp.unitDirection(), ref.unitDirection(), epsDir);
  CHECK_CLOSE_ABS(cmp.absoluteMomentum(), ref.absoluteMomentum(), epsMom);
  // charge should be identical not just similar
  BOOST_CHECK_EQUAL(cmp.charge(), ref.charge());
}

/// Check that two parameters covariances are consistent within the tolerances.
///
/// \warning Does not check that the parameters value itself are consistent.
template <typename charge_t>
inline void checkCovarianceConsistency(
    const Acts::SingleBoundTrackParameters<charge_t>& cmp,
    const Acts::SingleBoundTrackParameters<charge_t>& ref,
    double relativeTolerance) {
  // either both or none have covariance set
  if (cmp.covariance().has_value()) {
    // comparison parameters have covariance but the reference does not
    BOOST_CHECK(ref.covariance().has_value());
  }
  if (ref.covariance().has_value()) {
    // reference parameters have covariance but the comparison does not
    BOOST_CHECK(cmp.covariance().has_value());
  }
  if (cmp.covariance().has_value() and ref.covariance().has_value()) {
    CHECK_CLOSE_COVARIANCE(cmp.covariance().value(), ref.covariance().value(),
                           relativeTolerance);
  }
}

// helpers to construct target surfaces from track states

/// Construct the transformation from the curvilinear to the global coordinates.
template <typename charge_t>
inline Acts::Transform3 makeCurvilinearTransform(
    const Acts::SingleBoundTrackParameters<charge_t>& params,
    const Acts::GeometryContext& geoCtx) {
  Acts::Vector3 unitW = params.unitDirection();
  auto [unitU, unitV] = Acts::makeCurvilinearUnitVectors(unitW);

  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Zero();
  rotation.col(0) = unitU;
  rotation.col(1) = unitV;
  rotation.col(2) = unitW;
  Acts::Translation3 offset(params.position(geoCtx));
  Acts::Transform3 toGlobal = offset * rotation;

  return toGlobal;
}

/// Construct a z-cylinder centered at zero with the track on its surface.
struct ZCylinderSurfaceBuilder {
  template <typename charge_t>
  std::shared_ptr<Acts::CylinderSurface> operator()(
      const Acts::SingleBoundTrackParameters<charge_t>& params,
      const Acts::GeometryContext& geoCtx) {
    auto radius = params.position(geoCtx).template head<2>().norm();
    auto halfz = std::numeric_limits<double>::max();
    return Acts::Surface::makeShared<Acts::CylinderSurface>(
        Acts::Transform3::Identity(), radius, halfz);
  }
};

/// Construct a disc at track position with plane normal along track tangent.
struct DiscSurfaceBuilder {
  template <typename charge_t>
  std::shared_ptr<Acts::DiscSurface> operator()(
      const Acts::SingleBoundTrackParameters<charge_t>& params,
      const Acts::GeometryContext& geoCtx) {
    using namespace Acts;
    using namespace Acts::UnitLiterals;

    auto cl = makeCurvilinearTransform(params, geoCtx);
    // shift the origin of the plane so the local particle position does not
    // sit directly at the rho=0,phi=undefined singularity
    // TODO this is a hack do avoid issues with the numerical covariance
    //      transport that does not work well at rho=0,
    Acts::Vector3 localOffset = Acts::Vector3::Zero();
    localOffset[Acts::ePos0] = 1_cm;
    localOffset[Acts::ePos1] = -1_cm;
    Acts::Vector3 globalOriginDelta = cl.linear() * localOffset;
    cl.pretranslate(globalOriginDelta);

    return Acts::Surface::makeShared<Acts::DiscSurface>(cl);
  }
};

/// Construct a plane at track position with plane normal along track tangent.
struct PlaneSurfaceBuilder {
  template <typename charge_t>
  std::shared_ptr<Acts::PlaneSurface> operator()(
      const Acts::SingleBoundTrackParameters<charge_t>& params,
      const Acts::GeometryContext& geoCtx) {
    return Acts::Surface::makeShared<Acts::PlaneSurface>(
        makeCurvilinearTransform(params, geoCtx));
  }
};

/// Construct a z-straw at the track position.
struct ZStrawSurfaceBuilder {
  template <typename charge_t>
  std::shared_ptr<Acts::StrawSurface> operator()(
      const Acts::SingleBoundTrackParameters<charge_t>& params,
      const Acts::GeometryContext& geoCtx) {
    return Acts::Surface::makeShared<Acts::StrawSurface>(
        Acts::Transform3(Acts::Translation3(params.position(geoCtx))));
  }
};

// helper functions to run the propagation with additional checks

/// Propagate the initial parameters for the given pathlength in space.
///
/// Use a negative path length to indicate backward propagation.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline std::pair<Acts::CurvilinearTrackParameters, double> transportFreely(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    double pathLength) {
  using namespace Acts::UnitLiterals;

  using Actions = Acts::ActionList<>;
  using Aborts = Acts::AbortList<>;

  // setup propagation options
  options_t<Actions, Aborts> options(geoCtx, magCtx);
  options.direction = (0 <= pathLength) ? Acts::NavigationDirection::Forward
                                        : Acts::NavigationDirection::Backward;
  options.pathLimit = pathLength;
  options.maxStepSize = 1_cm;

  auto result = propagator.propagate(initialParams, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  return {*result.value().endParameters, result.value().pathLength};
}

/// Propagate the initial parameters to the target surface.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline std::pair<Acts::BoundTrackParameters, double> transportToSurface(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    const Acts::Surface& targetSurface, double pathLimit) {
  using namespace Acts::UnitLiterals;

  using Actions = Acts::ActionList<>;
  using Aborts = Acts::AbortList<>;

  // setup propagation options
  options_t<Actions, Aborts> options(geoCtx, magCtx);
  options.direction = Acts::NavigationDirection::Forward;
  options.pathLimit = pathLimit;
  options.maxStepSize = 1_cm;

  auto result = propagator.propagate(initialParams, targetSurface, options);
  BOOST_CHECK(result.ok());
  BOOST_CHECK(result.value().endParameters);

  return {*result.value().endParameters, result.value().pathLength};
}

// self-consistency tests for a single propagator

/// Propagate the initial parameters the given path length along its
/// trajectory and then propagate the final parameters back. Verify that the
/// propagated parameters match the initial ones.
template <typename propagator_t, typename charge_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline void runForwardBackwardTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    double pathLength, double epsPos, double epsDir, double epsMom) {
  // propagate parameters NavigationDirection::Forward
  auto [fwdParams, fwdPathLength] =
      transportFreely<propagator_t, charge_t, options_t>(
          propagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(fwdPathLength, pathLength, epsPos);
  // propagate propagated parameters back again
  auto [bwdParams, bwdPathLength] =
      transportFreely<propagator_t, charge_t, options_t>(
          propagator, geoCtx, magCtx, fwdParams, -pathLength);
  CHECK_CLOSE_ABS(bwdPathLength, -pathLength, epsPos);
  // check that initial and back-propagated parameters match
  checkParametersConsistency(initialParams, bwdParams, geoCtx, epsPos, epsDir,
                             epsMom);
}

/// Propagate the initial parameters once for the given path length and
/// use the propagated parameters to define a target surface. Propagate the
/// initial parameters again to the target surface. Verify that the surface has
/// been found and the parameters are consistent.
template <typename propagator_t, typename charge_t, typename surface_builder_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline void runToSurfaceTest(
    const propagator_t& propagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    double pathLength, surface_builder_t&& buildTargetSurface, double epsPos,
    double epsDir, double epsMom) {
  // free propagation for the given path length
  auto [freeParams, freePathLength] =
      transportFreely<propagator_t, charge_t, options_t>(
          propagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(freePathLength, pathLength, epsPos);
  // build a target surface at the propagated position
  auto surface = buildTargetSurface(freeParams, geoCtx);
  BOOST_CHECK(surface);

  // bound propagation onto the target surface
  // increase path length limit to ensure the surface can be reached
  auto [surfParams, surfPathLength] =
      transportToSurface<propagator_t, charge_t, options_t>(
          propagator, geoCtx, magCtx, initialParams, *surface,
          1.5 * pathLength);
  CHECK_CLOSE_ABS(surfPathLength, pathLength, epsPos);

  // check that the to-surface propagation matches the defining free parameters
  CHECK_CLOSE_ABS(surfParams.position(geoCtx), freeParams.position(geoCtx),
                  epsPos);
  CHECK_CLOSE_ABS(surfParams.time(), freeParams.time(), epsPos);
  CHECK_CLOSE_ABS(surfParams.unitDirection(), freeParams.unitDirection(),
                  epsDir);
  CHECK_CLOSE_ABS(surfParams.absoluteMomentum(), freeParams.absoluteMomentum(),
                  epsMom);
  CHECK_CLOSE_ABS(surfPathLength, freePathLength, epsPos);
}

// consistency tests between two propagators

/// Propagate the initial parameters along their trajectory for the given path
/// length using two different propagators and verify consistent output.
template <
    typename cmp_propagator_t, typename ref_propagator_t, typename charge_t,
    template <typename, typename> class options_t = Acts::PropagatorOptions>
inline void runForwardComparisonTest(
    const cmp_propagator_t& cmpPropagator,
    const ref_propagator_t& refPropagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    double pathLength, double epsPos, double epsDir, double epsMom,
    double tolCov) {
  // propagate twice using the two different propagators
  auto [cmpParams, cmpPath] =
      transportFreely<cmp_propagator_t, charge_t, options_t>(
          cmpPropagator, geoCtx, magCtx, initialParams, pathLength);
  auto [refParams, refPath] =
      transportFreely<ref_propagator_t, charge_t, options_t>(
          refPropagator, geoCtx, magCtx, initialParams, pathLength);
  // check parameter comparison
  checkParametersConsistency(cmpParams, refParams, geoCtx, epsPos, epsDir,
                             epsMom);
  checkCovarianceConsistency(cmpParams, refParams, tolCov);
  CHECK_CLOSE_ABS(cmpPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(refPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(cmpPath, refPath, epsPos);
}

/// Propagate the initial parameters along their trajectory for the given path
/// length using the reference propagator. Use the propagated track parameters
/// to define a target plane. Propagate the initial parameters using two
/// different propagators and verify consistent output.
template <typename cmp_propagator_t, typename ref_propagator_t,
          typename charge_t, typename surface_builder_t,
          template <typename, typename>
          class options_t = Acts::PropagatorOptions>
inline void runToSurfaceComparisonTest(
    const cmp_propagator_t& cmpPropagator,
    const ref_propagator_t& refPropagator, const Acts::GeometryContext& geoCtx,
    const Acts::MagneticFieldContext& magCtx,
    const Acts::SingleCurvilinearTrackParameters<charge_t>& initialParams,
    double pathLength, surface_builder_t&& buildTargetSurface, double epsPos,
    double epsDir, double epsMom, double tolCov) {
  // free propagation with the reference propagator for the given path length
  auto [freeParams, freePathLength] =
      transportFreely<ref_propagator_t, charge_t, options_t>(
          refPropagator, geoCtx, magCtx, initialParams, pathLength);
  CHECK_CLOSE_ABS(freePathLength, pathLength, epsPos);

  // build a target surface at the propagated position
  auto surface = buildTargetSurface(freeParams, geoCtx);
  BOOST_CHECK(surface);

  // propagate twice to the surface using the two different propagators
  // increase path length limit to ensure the surface can be reached
  auto [cmpParams, cmpPath] =
      transportToSurface<cmp_propagator_t, charge_t, options_t>(
          cmpPropagator, geoCtx, magCtx, initialParams, *surface,
          1.5 * pathLength);
  auto [refParams, refPath] =
      transportToSurface<ref_propagator_t, charge_t, options_t>(
          refPropagator, geoCtx, magCtx, initialParams, *surface,
          1.5 * pathLength);
  // check parameter comparison
  checkParametersConsistency(cmpParams, refParams, geoCtx, epsPos, epsDir,
                             epsMom);
  checkCovarianceConsistency(cmpParams, refParams, tolCov);
  CHECK_CLOSE_ABS(cmpPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(refPath, pathLength, epsPos);
  CHECK_CLOSE_ABS(cmpPath, refPath, epsPos);
}
