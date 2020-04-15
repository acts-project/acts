// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
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
/// @param [in] state State that will be projected
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
  jacToCurv(3, 6) = -invSinTheta;
  jacToCurv(4, 7) = 1.;

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
const FreeToBoundMatrix surfaceDerivative(std::reference_wrapper<const GeometryContext> geoContext, const Vector3D& position,
										const Vector3D& direction, BoundToFreeMatrix& jacobianLocalToGlobal, 
											const FreeVector& derivatives,
                                          const Surface* surface = nullptr) {
  // Set the surface projection contributions
  // If no surface is specified it is curvilinear
  if (surface == nullptr) {
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
  // Else it is bound
  else {
    // Initialize the transport final frame jacobian
    FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
    // Initalize the jacobian to local, returns the transposed ref frame
    auto rframeT = surface->initJacobianToLocal(geoContext, jacToLocal,
                                                position, direction);
    // Calculate the form factors for the derivatives
    const BoundRowVector sVec = surface->derivativeFactors(
        geoContext, position, direction, rframeT, jacobianLocalToGlobal);
    jacobianLocalToGlobal -= derivatives * sVec;
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
void reinitializeJacToGlobal(std::reference_wrapper<const GeometryContext> geoContext, BoundToFreeMatrix& jacobianLocalToGlobal, const Vector3D& position, const Vector3D& direction,
							double charge, double momentum, double time,
                             const Surface* surface = nullptr) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobian
  jacobianLocalToGlobal = BoundToFreeMatrix::Zero();

  // Fill the jacobian to global for next transport
  // If treating curvilinear parameters
  if (surface == nullptr) {
    // TODO: This was calculated before - can it be reused?
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
  // If treating bound parameters
  else {
    Vector2D loc{0., 0.};
    surface->globalToLocal(geoContext, position, direction, loc);
    BoundVector pars;
    pars << loc[eLOC_0], loc[eLOC_1], phi(direction), theta(direction),
        charge / momentum, time;
    surface->initJacobianToGlobal(geoContext, jacobianLocalToGlobal,
                                  position, direction, pars);
  }
}

/// @brief Reinitializes the jacobians of @p state and its components
///
/// @param [in, out] state The state of the stepper
/// @param [in] surface Representing surface of the stepper state
void reinitializeJacobians(std::reference_wrapper<const GeometryContext> geoContext, 
							FreeMatrix& transportJacobian, FreeVector& derivatives, BoundToFreeMatrix& jacobianLocalToGlobal, 
							const Vector3D& position, const Vector3D& direction,
							double charge, double momentum, double time,
                           const Surface* surface = nullptr) {
  transportJacobian = FreeMatrix::Identity();
  derivatives = FreeVector::Zero();
  reinitializeJacToGlobal(geoContext, jacobianLocalToGlobal, position, direction, charge, momentum, time, surface);
}
}  // namespace

namespace detail {

BoundState boundState(std::reference_wrapper<const GeometryContext> geoContext, Covariance& covarianceMatrix, Jacobian& jacobian,
							FreeMatrix& transportJacobian, FreeVector& derivatives, BoundToFreeMatrix& jacobianLocalToGlobal, 
							const Vector3D& position, const Vector3D& direction,
							double charge, double momentum, double time, bool covTransport, double accumulatedPath,
                         const Surface& surface,  bool reinitialize
                         ) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    // Initialize the transport final frame jacobian
    covarianceTransport(geoContext, covarianceMatrix, jacobian, transportJacobian, derivatives, jacobianLocalToGlobal, position, direction, charge, momentum, time, reinitialize, &surface);
    cov = covarianceMatrix;
  }
  // Create the bound parameters
  BoundParameters parameters(geoContext, cov, position,
                             momentum * direction, charge, time,
                             surface.getSharedPtr());
  // Create the bound state
  BoundState result = std::make_tuple(std::move(parameters), jacobian,
                                      accumulatedPath);
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    jacobian = Jacobian::Identity();
  }
  return result;
}

CurvilinearState curvilinearState(std::reference_wrapper<const GeometryContext> geoContext, Covariance& covarianceMatrix, Jacobian& jacobian,
							FreeMatrix& transportJacobian, FreeVector& derivatives, BoundToFreeMatrix& jacobianLocalToGlobal, 
							const Vector3D& position, const Vector3D& direction,
							double charge, double momentum, double time, bool covTransport, double accumulatedPath,
                         bool reinitialize) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    covarianceTransport(geoContext, covarianceMatrix, jacobian, transportJacobian, derivatives, jacobianLocalToGlobal, position, direction, charge, momentum, time, reinitialize, &surface);
    cov = covarianceMatrix;
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(cov, position, momentum * direction, charge,
                                   time);
  // Create the bound state
  CurvilinearState result = std::make_tuple(
      std::move(parameters), jacobian, accumulatedPath);
  // Reinitialize if asked to do so
  // this is useful for interruption calls
  if (reinitialize) {
    jacobian = Jacobian::Identity();
  }
  return result;
}

void covarianceTransport(
std::reference_wrapper<const GeometryContext> geoContext, Covariance& covarianceMatrix, Jacobian& jacobian,
							FreeMatrix& transportJacobian, FreeVector& derivatives, BoundToFreeMatrix& jacobianLocalToGlobal, 
							const Vector3D& position, const Vector3D& direction,
							double charge, double momentum, double time,                           
                           bool reinitialize,
                         const Surface* surface) {
  jacobianLocalToGlobal = tansportJacobian * jacobianLocalToGlobal;

  const FreeToBoundMatrix jacToLocal = surfaceDerivative(geoContext, position, direction, jacobianLocalToGlobal, derivatives, surface);
  const Jacobian jacFull = jacToLocal * jacobianLocalToGlobal;

  // Apply the actual covariance transport
  covarianceMatrix = jacFull * covarianceMatrix * jacFull.transpose();

  if (reinitialize) {
    reinitializeJacobians(geoContext, transportJacobian, derivatives, jacobianLocalToGlobal, position, direction, charge, momentum, time, surface);
  }

  // Store The global and bound jacobian (duplication for the moment)
  jacobian = jacFull * jacobian;
}
}  // namespace detail
}  // namespace Acts