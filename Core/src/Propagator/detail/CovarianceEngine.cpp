// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;
using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearTrackParameters, Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const Vector3D& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to bound parameters at the final surface
///
// @note Modifications of the jacobian related to the
/// projection onto a surface is considered. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] jacobianLocalToGlobal The projection jacobian from start local
/// to start free parameters
/// @param [in] transportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] derivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] jacFull The full jacobian from start local to bound
/// parameters at the final surface
/// @param [in] surface The final surface onto which the projection should be
/// performed
void localToLocalJacobian(
    std::reference_wrapper<const GeometryContext> geoContext,
    const FreeVector& parameters,
    const BoundToFreeMatrix& jacobianLocalToGlobal,
    const FreeMatrix& transportJacobian, const FreeVector& derivatives,
    Jacobian& jacFull, const Surface& surface) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeRowVector freeToPath =
      surface.freeToPathDerivative(geoContext, parameters);
  // Initalize the jacobian from global to local at the final surface
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Calculate the global to local
  surface.initJacobianToLocal(geoContext, jacToLocal,
                              parameters.segment<3>(eFreePos0),
                              parameters.segment<3>(eFreeDir0));
  // Calculate the full jocobian from the local parameters at the start surface
  // to local parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*transportJac(gloA->gloB) *jac(locA->gloA)
  jacFull = jacToLocal * (FreeMatrix::Identity() + derivatives * freeToPath) *
            transportJacobian * jacobianLocalToGlobal;
}

/// @brief This function calculates the full jacobian from local parameters at
/// the start surface to final curvilinear parameters
///
/// @note Modifications of the jacobian related to the
/// projection onto a curvilinear surface is considered. Since a variation of
/// the start parameters within a given uncertainty would lead to a variation of
/// the end parameters, these need to be propagated onto the target surface.
/// This an approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] jacobianLocalToGlobal The projection jacobian from local start
/// to global final parameters
/// @param [in] transportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in] derivatives Path length derivatives of the final free
/// parameters
/// @param [in, out] jacFull The full jacobian from start local to curvilinear
/// parameters
///
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
void localToLocalJacobian(const Vector3D& direction,
                          const BoundToFreeMatrix& jacobianLocalToGlobal,
                          const FreeMatrix& transportJacobian,
                          const FreeVector& derivatives, Jacobian& jacFull) {
  // Calculate the derivative of path length at the the curvilinear surface
  // w.r.t. free parameters
  FreeRowVector freeToPath = FreeRowVector::Zero();
  freeToPath.head<3>() = -1.0 * direction;
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix jacToLocal = freeToCurvilinearJacobian(direction);
  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*transportJac(gloA->gloB) *jac(locA->gloA)
  jacFull = jacToLocal * (FreeMatrix::Identity() + derivatives * freeToPath) *
            transportJacobian * jacobianLocalToGlobal;
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] transportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] surface The surface the represents the local parametrisation
void reinitializeJacobians(
    std::reference_wrapper<const GeometryContext> geoContext,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacobianLocalToGlobal, const FreeVector& parameters,
    const Surface& surface) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobians
  transportJacobian = FreeMatrix::Identity();
  derivatives = FreeVector::Zero();
  jacobianLocalToGlobal = BoundToFreeMatrix::Zero();

  // Reset the jacobian from local to global
  const Vector3D position = parameters.segment<3>(eFreePos0);
  const Vector3D direction = parameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (not lpResult.ok()) {
    ACTS_LOCAL_LOGGER(
        Acts::getDefaultLogger("CovarianceEngine", Logging::INFO));
    ACTS_FATAL(
        "Inconsistency in global to local transformation during propagation.")
  }
  auto loc = lpResult.value();
  BoundVector pars;
  pars << loc[eBoundLoc0], loc[eBoundLoc1], phi(direction), theta(direction),
      parameters[eFreeQOverP], parameters[eFreeTime];
  surface.initJacobianToGlobal(geoContext, jacobianLocalToGlobal, position,
                               direction, pars);
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] transportJacobian The transport jacobian from start free to
/// final free parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& transportJacobian,
                           FreeVector& derivatives,
                           BoundToFreeMatrix& jacobianLocalToGlobal,
                           const Vector3D& direction) {
  // Reset the jacobians
  transportJacobian = FreeMatrix::Identity();
  derivatives = FreeVector::Zero();
  jacobianLocalToGlobal = BoundToFreeMatrix::Zero();

  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;

  jacobianLocalToGlobal(0, eBoundLoc0) = -sinPhi;
  jacobianLocalToGlobal(0, eBoundLoc1) = -cosPhi * cosTheta;
  jacobianLocalToGlobal(1, eBoundLoc0) = cosPhi;
  jacobianLocalToGlobal(1, eBoundLoc1) = -sinPhi * cosTheta;
  jacobianLocalToGlobal(2, eBoundLoc1) = sinTheta;
  jacobianLocalToGlobal(3, eBoundTime) = 1;
  jacobianLocalToGlobal(4, eBoundPhi) = -sinTheta * sinPhi;
  jacobianLocalToGlobal(4, eBoundTheta) = cosTheta * cosPhi;
  jacobianLocalToGlobal(5, eBoundPhi) = sinTheta * cosPhi;
  jacobianLocalToGlobal(5, eBoundTheta) = cosTheta * sinPhi;
  jacobianLocalToGlobal(6, eBoundTheta) = -sinTheta;
  jacobianLocalToGlobal(7, eBoundQOverP) = 1;
}
}  // namespace

namespace detail {

BoundState boundState(std::reference_wrapper<const GeometryContext> geoContext,
                      Covariance& covarianceMatrix, Jacobian& jacobian,
                      FreeMatrix& transportJacobian, FreeVector& derivatives,
                      BoundToFreeMatrix& jacobianLocalToGlobal,
                      const FreeVector& parameters, bool covTransport,
                      double accumulatedPath, const Surface& surface) {
  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = Jacobian::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // jacobianLocalToGlobal
    covarianceTransport(geoContext, covarianceMatrix, jacobian,
                        transportJacobian, derivatives, jacobianLocalToGlobal,
                        parameters, surface);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the bound parameters
  BoundVector bv =
      detail::transformFreeToBoundParameters(parameters, surface, geoContext);
  // Create the bound state
  return std::make_tuple(
      BoundTrackParameters(surface.getSharedPtr(), bv, std::move(cov)),
      jacobian, accumulatedPath);
}

CurvilinearState curvilinearState(Covariance& covarianceMatrix,
                                  Jacobian& jacobian,
                                  FreeMatrix& transportJacobian,
                                  FreeVector& derivatives,
                                  BoundToFreeMatrix& jacobianLocalToGlobal,
                                  const FreeVector& parameters,
                                  bool covTransport, double accumulatedPath) {
  const Vector3D& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the jacobian from start local to final local
    jacobian = Jacobian::Identity();
    // Calculate the jacobian and transport the covarianceMatrix to final local.
    // Then reinitialize the transportJacobian, derivatives and the
    // jacobianLocalToGlobal
    covarianceTransport(covarianceMatrix, jacobian, transportJacobian,
                        derivatives, jacobianLocalToGlobal, direction);
  }
  if (covarianceMatrix != BoundSymMatrix::Zero()) {
    cov = covarianceMatrix;
  }

  // Create the curvilinear parameters
  Vector4D pos4 = Vector4D::Zero();
  pos4[ePos0] = parameters[eFreePos0];
  pos4[ePos1] = parameters[eFreePos1];
  pos4[ePos2] = parameters[eFreePos2];
  pos4[eTime] = parameters[eFreeTime];
  CurvilinearTrackParameters curvilinearParams(
      pos4, direction, parameters[eFreeQOverP], std::move(cov));
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParams), jacobian,
                         accumulatedPath);
}

void covarianceTransport(Covariance& covarianceMatrix, Jacobian& jacobian,
                         FreeMatrix& transportJacobian, FreeVector& derivatives,
                         BoundToFreeMatrix& jacobianLocalToGlobal,
                         const Vector3D& direction) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current curvilinear parameters
  localToLocalJacobian(direction, jacobianLocalToGlobal, transportJacobian,
                       derivatives, jacobian);

  // Apply the actual covariance transport to get covariance of the current
  // curvilinear parameters
  covarianceMatrix = jacobian * covarianceMatrix * jacobian.transpose();

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The jacobianLocalToGlobal is reinitialized to that at the current
  // curvilinear surface
  reinitializeJacobians(transportJacobian, derivatives, jacobianLocalToGlobal,
                        direction);
}

void covarianceTransport(
    std::reference_wrapper<const GeometryContext> geoContext,
    Covariance& covarianceMatrix, Jacobian& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacobianLocalToGlobal, const FreeVector& parameters,
    const Surface& surface) {
  // Calculate the full jacobian from local parameters at the start surface to
  // current bound parameters
  localToLocalJacobian(geoContext, parameters, jacobianLocalToGlobal,
                       transportJacobian, derivatives, jacobian, surface);

  // Apply the actual covariance transport to get covariance of the current
  // bound parameters
  covarianceMatrix = jacobian * covarianceMatrix * jacobian.transpose();

  // Reinitialize jacobian components:
  // ->The transportJacobian is reinitialized to Identity
  // ->The derivatives is reinitialized to Zero
  // ->The jacobianLocalToGlobal is initialized to that at the current surface
  reinitializeJacobians(geoContext, transportJacobian, derivatives,
                        jacobianLocalToGlobal, parameters, surface);
}

}  // namespace detail
}  // namespace Acts
