// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {
namespace {
/// Some type defs
using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;
using BoundState = std::tuple<BoundParameters, Jacobian, double>;
using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;

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

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a surface. Since a variation of the start parameters within
/// a given uncertainty would lead to a variation of the end parameters, these
/// need to be propagated onto the target surface. This an approximated approach
/// to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] jacobianLocalToGlobal The projection jacobian from local start
/// to global final parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in] surface The surface onto which the projection should be
/// performed
///
/// @return The projection jacobian from global end parameters to its local
/// equivalentconst
FreeToBoundMatrix surfaceDerivative(
    std::reference_wrapper<const GeometryContext> geoContext,
    const FreeVector& parameters, BoundToFreeMatrix& jacobianLocalToGlobal,
    const FreeVector& derivatives, const Surface& surface) {
  // Initialize the transport final frame jacobian
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Initalize the jacobian to local, returns the transposed ref frame
  auto rframeT = surface.initJacobianToLocal(geoContext, jacToLocal,
                                             parameters.segment<3>(eFreePos0),
                                             parameters.segment<3>(eFreeDir0));
  // Calculate the form factors for the derivatives
  const BoundRowVector sVec = surface.derivativeFactors(
      geoContext, parameters.segment<3>(eFreePos0),
      parameters.segment<3>(eFreeDir0), rframeT, jacobianLocalToGlobal);
  jacobianLocalToGlobal -= derivatives * sVec;
  // Return the jacobian to local
  return jacToLocal;
}

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a curvilinear surface. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] jacobianLocalToGlobal The projection jacobian from local start
/// to global final parameters
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
///
/// @return The projection jacobian from global end parameters to its local
/// equivalent
const FreeToBoundMatrix surfaceDerivative(
    const Vector3D& direction, BoundToFreeMatrix& jacobianLocalToGlobal,
    const FreeVector& derivatives) {
  // Transport the covariance
  const ActsRowVectorD<3> normVec(direction);
  const BoundRowVector sfactors =
      normVec *
      jacobianLocalToGlobal.template topLeftCorner<3, eBoundParametersSize>();
  jacobianLocalToGlobal -= derivatives * sfactors;
  // Since the jacobian to local needs to calculated for the bound parameters
  // here, it is convenient to do the same here
  return freeToCurvilinearJacobian(direction);
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in] geoContext The geometry context
/// @param [in, out] jacobian Full jacobian since the last reset
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
  Vector2D loc{0., 0.};
  const Vector3D position = parameters.segment<3>(eFreePos0);
  const Vector3D direction = parameters.segment<3>(eFreeDir0);
  surface.globalToLocal(geoContext, position, direction, loc);
  BoundVector pars;
  pars << loc[eLOC_0], loc[eLOC_1], phi(direction), theta(direction),
      parameters[eFreeQOverP], parameters[eFreeTime];
  surface.initJacobianToGlobal(geoContext, jacobianLocalToGlobal, position,
                               direction, pars);
}

/// @brief This function reinitialises the state members required for the
/// covariance transport
///
/// @param [in, out] jacobian Full jacobian since the last reset
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

  jacobianLocalToGlobal(0, eLOC_0) = -sinPhi;
  jacobianLocalToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
  jacobianLocalToGlobal(1, eLOC_0) = cosPhi;
  jacobianLocalToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
  jacobianLocalToGlobal(2, eLOC_1) = sinTheta;
  jacobianLocalToGlobal(3, eT) = 1;
  jacobianLocalToGlobal(4, ePHI) = -sinTheta * sinPhi;
  jacobianLocalToGlobal(4, eTHETA) = cosTheta * cosPhi;
  jacobianLocalToGlobal(5, ePHI) = sinTheta * cosPhi;
  jacobianLocalToGlobal(5, eTHETA) = cosTheta * sinPhi;
  jacobianLocalToGlobal(6, eTHETA) = -sinTheta;
  jacobianLocalToGlobal(7, eQOP) = 1;
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
    covarianceTransport(geoContext, covarianceMatrix, jacobian,
                        transportJacobian, derivatives, jacobianLocalToGlobal,
                        parameters, surface);
    cov = covarianceMatrix;
  }
  // Create the bound parameters
  const Vector3D& position = parameters.segment<3>(eFreePos0);
  const Vector3D momentum =
      std::abs(1. / parameters[eFreeQOverP]) * parameters.segment<3>(eFreeDir0);
  const double charge = std::copysign(1., parameters[eFreeQOverP]);
  const double time = parameters[eFreeTime];
  BoundParameters boundParameters(geoContext, cov, position, momentum, charge,
                                  time, surface.getSharedPtr());
  // Create the bound state
  return std::make_tuple(std::move(boundParameters), jacobian, accumulatedPath);
  ;
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
    covarianceTransport(covarianceMatrix, jacobian, transportJacobian,
                        derivatives, jacobianLocalToGlobal, direction);
    cov = covarianceMatrix;
  }
  // Create the curvilinear parameters
  const Vector3D& position = parameters.segment<3>(eFreePos0);
  const Vector3D momentum = std::abs(1. / parameters[eFreeQOverP]) * direction;
  const double charge = std::copysign(1., parameters[eFreeQOverP]);
  const double time = parameters[eFreeTime];
  CurvilinearParameters curvilinearParameters(cov, position, momentum, charge,
                                              time);
  // Create the curvilinear state
  return std::make_tuple(std::move(curvilinearParameters), jacobian,
                         accumulatedPath);
}

void covarianceTransport(Covariance& covarianceMatrix, Jacobian& jacobian,
                         FreeMatrix& transportJacobian, FreeVector& derivatives,
                         BoundToFreeMatrix& jacobianLocalToGlobal,
                         const Vector3D& direction) {
  // Build the full jacobian
  jacobianLocalToGlobal = transportJacobian * jacobianLocalToGlobal;
  const FreeToBoundMatrix jacToLocal =
      surfaceDerivative(direction, jacobianLocalToGlobal, derivatives);
  const Jacobian jacFull = jacToLocal * jacobianLocalToGlobal;

  // Apply the actual covariance transport
  covarianceMatrix = jacFull * covarianceMatrix * jacFull.transpose();

  // Reinitialize jacobian components
  reinitializeJacobians(transportJacobian, derivatives, jacobianLocalToGlobal,
                        direction);

  // Store The global and bound jacobian (duplication for the moment)
  jacobian = jacFull;
}

void covarianceTransport(
    std::reference_wrapper<const GeometryContext> geoContext,
    Covariance& covarianceMatrix, Jacobian& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    BoundToFreeMatrix& jacobianLocalToGlobal, const FreeVector& parameters,
    const Surface& surface) {
  // Build the full jacobian
  jacobianLocalToGlobal = transportJacobian * jacobianLocalToGlobal;
  const FreeToBoundMatrix jacToLocal = surfaceDerivative(
      geoContext, parameters, jacobianLocalToGlobal, derivatives, surface);
  const Jacobian jacFull = jacToLocal * jacobianLocalToGlobal;

  // Apply the actual covariance transport
  covarianceMatrix = jacFull * covarianceMatrix * jacFull.transpose();

  // Reinitialize jacobian components
  reinitializeJacobians(geoContext, transportJacobian, derivatives,
                        jacobianLocalToGlobal, parameters, surface);

  // Store The global and bound jacobian (duplication for the moment)
  jacobian = jacFull;
}
}  // namespace detail
}  // namespace Acts