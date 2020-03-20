// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
namespace {

/// Projection Jacobian from bound parameters on the surface to free parameters.
BoundToFreeMatrix boundToFreeJacobian(const StepperState& state,
                                      const Surface& surface) {
  // convert free to local parameters position
  Vector2D loc = Vector2D::Zero();
  surface.globalToLocal(state.geoContext, state.pos, state.dir, loc);
  BoundVector pars;
  pars[eBoundLoc0] = loc[0];
  pars[eBoundLoc1] = loc[1];
  pars[eBoundTime] = state.t;
  pars[eBoundPhi] = VectorHelpers::phi(state.dir);
  pars[eBoundTheta] = VectorHelpers::theta(state.dir);
  pars[eBoundQOverP] = (state.q != 0) ? (state.q / state.p) : (1 / state.p);
  BoundToFreeMatrix jacToGlobal;
  surface.initJacobianToGlobal(state.geoContext, jacToGlobal, state.pos,
                               state.dir, pars);
  return jacToGlobal;
}

/// Projection Jacobian from curvilinear parameters to free parameters.
///
/// In principle, this should be equivalent to the projection Jacobian for
/// general bound parameters with an appropriate Curvilinear surface. Use an
/// explicit implementation for faster computation.
BoundToFreeMatrix curvilinearToFreeJacobian(const Vector3D& unitDirection) {
  // normalized direction vector components expressed as phi/theta
  const double x = unitDirection[0];  // == cos(phi) * sin(theta)
  const double y = unitDirection[1];  // == sin(phi) * sin(theta)
  const double z = unitDirection[2];  // == cos(theta)
  // extract phi/theta expressions
  const auto cosTheta = z;
  const auto sinTheta = std::hypot(x, y);
  const auto invSinTheta = 1 / sinTheta;
  const auto cosPhi = x * invSinTheta;
  const auto sinPhi = y * invSinTheta;
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  jacToGlobal(eFreePos0, eBoundLoc0) = -sinPhi;
  jacToGlobal(eFreePos0, eBoundLoc1) = -cosPhi * cosTheta;
  jacToGlobal(eFreePos1, eBoundLoc0) = cosPhi;
  jacToGlobal(eFreePos1, eBoundLoc1) = -sinPhi * cosTheta;
  jacToGlobal(eFreePos2, eBoundLoc1) = sinTheta;
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  jacToGlobal(eFreeDir0, eBoundPhi) = -sinTheta * sinPhi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cosTheta * cosPhi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sinTheta * cosPhi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cosTheta * sinPhi;
  jacToGlobal(eFreeDir2, eBoundTheta) = -sinTheta;
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

/// Projection Jacobian from free to curvilinear parameters.
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3D& unitDirection) {
  // normalized direction vector components expressed as phi/theta
  const double x = unitDirection[0];  // == cos(phi) * sin(theta)
  const double y = unitDirection[1];  // == sin(phi) * sin(theta)
  const double z = unitDirection[2];  // == cos(theta)
  // extract phi/theta expressions
  const auto cosTheta = z;
  const auto sinTheta = std::hypot(x, y);
  const auto invSinTheta = 1 / sinTheta;
  const auto cosPhi = x * invSinTheta;
  const auto sinPhi = y * invSinTheta;
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as above
    jacToCurv(eBoundLoc0, eFreePos0) = -sinPhi;
    jacToCurv(eBoundLoc0, eFreePos1) = cosPhi;
    jacToCurv(eBoundLoc1, eFreePos0) = -cosPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos1) = -sinPhi * cosTheta;
    jacToCurv(eBoundLoc1, eFreePos2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const auto c = std::hypot(y, z);
    const auto invC = 1 / c;
    jacToCurv(eBoundLoc0, eFreePos1) = -z * invC;
    jacToCurv(eBoundLoc0, eFreePos2) = y * invC;
    jacToCurv(eBoundLoc1, eFreePos0) = c;
    jacToCurv(eBoundLoc1, eFreePos1) = -x * y * invC;
    jacToCurv(eBoundLoc1, eFreePos2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(eT, eFreeTime) = 1;
  // Directional and momentum parameters for curvilinear
  // TODO this becomes unstable for direction along z
  //      sin(theta) -> 0 -> 1/sin(theta) -> inf
  jacToCurv(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  jacToCurv(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  jacToCurv(eBoundTheta, eFreeDir2) = -invSinTheta;
  jacToCurv(eBoundQOverP, eFreeQOverP) = 1;
  return jacToCurv;
}

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a surface. Since a variation of the start parameters within
/// a given uncertainty would lead to a variation of the end parameters, these
/// need to be propagated onto the target surface. This an approximated approach
/// to treat the (assumed) small change.
///
/// @param [in] state The current state
/// @param [in] surface The surface onto which the projection should be
/// performed
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
///
/// @return The projection jacobian from global end parameters to its local
/// equivalent
const FreeToBoundMatrix surfaceDerivative(StepperState& state,
                                          const Surface* surface = nullptr) {
  // Set the surface projection contributions
  // If no surface is specified it is curvilinear
  if (surface == nullptr) {
    // Transport the covariance
    const ActsRowVectorD<3> normVec(state.dir);
    const BoundRowVector sfactors =
        normVec * state.jacToGlobal.template topLeftCorner<3, BoundParsDim>();
    state.jacToGlobal -= state.derivative * sfactors;
    // Since the jacobian to local needs to calculated for the bound parameters
    // here, it is convenient to do the same here
    return freeToCurvilinearJacobian(state);
  }
  // Else it is bound
  else {
    // Initialize the transport final frame jacobian
    FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
    // Initalize the jacobian to local, returns the transposed ref frame
    auto rframeT = surface->initJacobianToLocal(state.geoContext, jacToLocal,
                                                state.pos, state.dir);
    // Calculate the form factors for the derivatives
    const BoundRowVector sVec = surface->derivativeFactors(
        state.geoContext, state.pos, state.dir, rframeT, state.jacToGlobal);
    state.jacToGlobal -= state.derivative * sVec;
    // Return the jacobian to local
    return jacToLocal;
  }
}

/// @brief This function reinitialises the @p state member @p jacToGlobal.
///
/// @param [in, out] state The state object
/// @param [in] surface The surface the represents the local parametrisation
/// @note The surface is only required for bound parameters since it serves to
/// derive the jacobian from it. In the case of curvilinear parameters this is
/// not needed and can be evaluated without any surface.
void reinitializeJacToGlobal(StepperState& state,
                             const Surface* surface = nullptr) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobian
  state.jacToGlobal = BoundToFreeMatrix::Zero();

  // Fill the jacobian to global for next transport
  // If treating curvilinear parameters
  if (surface == nullptr) {
    // TODO: This was calculated before - can it be reused?
    // Optimized trigonometry on the propagation direction
    const double x = state.dir(0);  // == cos(phi) * sin(theta)
    const double y = state.dir(1);  // == sin(phi) * sin(theta)
    const double z = state.dir(2);  // == cos(theta)
    // can be turned into cosine/sine
    const double cosTheta = z;
    const double sinTheta = sqrt(x * x + y * y);
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = x * invSinTheta;
    const double sinPhi = y * invSinTheta;

    state.jacToGlobal(0, eLOC_0) = -sinPhi;
    state.jacToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
    state.jacToGlobal(1, eLOC_0) = cosPhi;
    state.jacToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
    state.jacToGlobal(2, eLOC_1) = sinTheta;
    state.jacToGlobal(3, eT) = 1;
    state.jacToGlobal(4, ePHI) = -sinTheta * sinPhi;
    state.jacToGlobal(4, eTHETA) = cosTheta * cosPhi;
    state.jacToGlobal(5, ePHI) = sinTheta * cosPhi;
    state.jacToGlobal(5, eTHETA) = cosTheta * sinPhi;
    state.jacToGlobal(6, eTHETA) = -sinTheta;
    state.jacToGlobal(7, eQOP) = 1;
  }
  // If treating bound parameters
  else {
    Vector2D loc{0., 0.};
    surface->globalToLocal(state.geoContext, state.pos, state.dir, loc);
    BoundVector pars;
    pars << loc[eLOC_0], loc[eLOC_1], phi(state.dir), theta(state.dir),
        state.q / state.p, state.t;
    surface->initJacobianToGlobal(state.geoContext, state.jacToGlobal,
                                  state.pos, state.dir, pars);
  }
}

/// @brief Reinitializes the jacobians of @p state and its components
///
/// @param [in, out] state The state of the stepper
/// @param [in] surface Representing surface of the stepper state
void reinitializeJacobians(StepperState& state,
                           const Surface* surface = nullptr) {
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
  reinitializeJacToGlobal(state, surface);
}
}  // namespace

namespace detail {

std::tuple<BoundParameters, BoundMatrix, double> boundState(
    StepperState& state, const Surface& surface) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    // Initialize the transport final frame jacobian
    covarianceTransport(state, &surface);
    cov = state.cov;
  }
  // Create the bound parameters
  BoundParameters parameters(state.geoContext, std::move(cov), state.pos,
                             state.p * state.dir, state.q, state.t,
                             surface.getSharedPtr());
  // Create the bound state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

std::tuple<CurvilinearParameters, BoundMatrix, double> curvilinearState(
    StepperState& state) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state);
    cov = state.cov;
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(std::move(cov), state.pos,
                                   state.p * state.dir, state.q, state.t);
  // Create the bound state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

void covarianceTransport(StepperState& state, const Surface* surface) {
  state.jacToGlobal = state.jacTransport * state.jacToGlobal;

  const FreeToBoundMatrix jacToLocal = surfaceDerivative(state, surface);
  const Jacobian jacFull = jacToLocal * state.jacToGlobal;

  // Apply the actual covariance transport
  state.cov = jacFull * state.cov * jacFull.transpose();

  // Reinitialize
  reinitializeJacobians(state, surface);

  // Store The global and bound jacobian (duplication for the moment)
  state.jacobian = jacFull;
}
}  // namespace detail
}  // namespace Acts
