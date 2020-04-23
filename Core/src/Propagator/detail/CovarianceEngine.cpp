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
using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;
using BoundState = std::tuple<BoundTrackParameters, detail::Jacobian, double>;
using CurvilinearState = std::tuple<CurvilinearTrackParameters, detail::Jacobian, double>;
using FreeState = std::tuple<FreeTrackParameters, detail::Jacobian, double>;

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
  jacToCurv(3, 6) = -invSinTheta * (1. - cosTheta * cosTheta);
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
/// equivalent
FreeToBoundMatrix surfaceDerivative(
    std::reference_wrapper<const GeometryContext> geoContext,
    const FreeVector& parameters, BoundToFreeMatrix& jacobianLocalToGlobal,
    const FreeVector& derivatives, const Surface& surface) {
  // Initialize the transport final frame jacobian
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Initalize the jacobian to local
  surface.initJacobianToLocal(geoContext, jacToLocal,
                              parameters.segment<3>(eFreePos0),
                              parameters.segment<3>(eFreeDir0));
if(jacobianLocalToGlobal.has_value())
{
  // Calculate the form factors for the derivatives
  const FreeRowVector freeToPath =
      surface.freeToPathDerivative(geoContext, parameters);
  jacobianLocalToGlobal += derivatives * freeToPath * jacobianLocalToGlobal;
  // Return the jacobian to local
  return jacToLocal;
}

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a surface. Since a variation of the start parameters within
/// a given uncertainty would lead to a variation of the end parameters, these
/// need to be propagated onto the target surface. This an approximated approach
/// to treat the (assumed) small change.
///
/// @param [in] geoContext The geometry Context
/// @param [in] parameters Free, nominal parametrisation
/// @param [in] jacobianDirToAngle Projection from direction dependent parametrisation to angular
/// @param [in] jacobianAngleToDir Projection from angular dependent parametrisation to direction 
/// @param [in] transportJacobian The transport jacobian 
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in] surface The surface onto which the projection should be
/// performed
///
/// @return The projection jacobian from global end parameters to its local
/// equivalent
FreeToBoundMatrix surfaceDerivative(
    std::reference_wrapper<const GeometryContext> geoContext,
    const FreeVector& parameters, const ActsMatrixD<8, 7>& jacobianDirToAngle, const ActsMatrixD<7, 8>& jacobianAngleToDir, const FreeMatrix& transportJacobian,
    const FreeVector& derivatives, const Surface& surface) {
  // Initialize the transport final frame jacobian
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Initalize the jacobian to local, returns the transposed ref frame
  auto rframeT = surface.initJacobianToLocal(geoContext, jacToLocal,
                                             parameters.segment<3>(eFreePos0),
                                             parameters.segment<3>(eFreeDir0));

	// Calculate the form factors for the derivatives
	const ActsMatrixD<8,7> transport = transportJacobian * jacobianDirToAngle;
	const ActsRowVectorD<7> sVec = surface.derivativeFactors(
		geoContext, parameters.segment<3>(eFreePos0),
        parameters.segment<3>(eFreeDir0), rframeT, transport);
	// Return the jacobian to local
	return jacToLocal * (transport - derivatives * sVec) * jacobianAngleToDir;
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
  const ActsRowVectorD<3> normVec(direction);
  const BoundRowVector sfactors =
      normVec *
      jacobianLocalToGlobal.template topLeftCorner<3, eBoundSize>();
  jacobianLocalToGlobal -= derivatives * sfactors;
  // Since the jacobian to local needs to calculated for the bound parameters
  // here, it is convenient to do the same here
  return freeToCurvilinearJacobian(direction);
}

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a curvilinear surface. Since a variation of the start
/// parameters within a given uncertainty would lead to a variation of the end
/// parameters, these need to be propagated onto the target surface. This an
/// approximated approach to treat the (assumed) small change.
///
/// @param [in] direction Normalised direction vector
/// @param [in] jacobianDirToAngle Projection from direction dependent parametrisation to angular
/// @param [in] jacobianAngleToDir Projection from angular dependent parametrisation to direction 
/// @param [in] transportJacobian The transport jacobian 
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
///
/// @return The projection jacobian from global end parameters to its local
/// equivalent
const FreeToBoundMatrix surfaceDerivative(
    const Vector3D& direction, const ActsMatrixD<8, 7>& jacobianDirToAngle, const ActsMatrixD<7, 8>& jacobianAngleToDir, const FreeMatrix& transportJacobian,
    const FreeVector& derivatives) {
  const ActsRowVectorD<3> normVec(direction);
	const ActsMatrixD<8,7> transport = transportJacobian * jacobianDirToAngle;
	const ActsRowVectorD<7> sfactors =
		normVec * transport.template topLeftCorner<3, 7>();
	// Since the jacobian to local needs to calculated for the bound parameters here, it is convenient to do the same here
	return freeToCurvilinearJacobian(direction) * (transport - derivatives * sfactors) * jacobianAngleToDir;
} 
  
/// @brief This function reinitialises the state members required for the
/// covariance transport in case of a bound target
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
/// covariance transport in case of a curvilinear target
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

/// @brief This function reinitialises the state members required for the
/// covariance transport in case of a free target
///
/// @param [in, out] jacobian Full jacobian since the last reset
/// @param [in, out] derivatives Path length derivatives of the free, nominal
/// parameters
/// @param [in, out] jacobianLocalToGlobal Projection jacobian of the last bound
/// parametrisation to free parameters
/// @param [in] direction Normalised direction vector
void reinitializeJacobians(FreeMatrix& transportJacobian,
                           FreeVector& derivatives,
                           std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal,
                           ActsMatrixD<8, 7>& jacobianDirToAngle, ActsMatrixD<7, 8>& jacobianAngleToDir,
                           const Vector3D& direction) {
	  // Reset the jacobians
  transportJacobian = FreeMatrix::Identity();
  derivatives = FreeVector::Zero();
  jacobianLocalToGlobal = std::nullopt;
  jacobianDirToAngle = detail::jacobianDirectionsToAngles(direction);
  jacobianAngleToDir = detail::jacobianAnglesToDirections(direction);
						}
}  // namespace

namespace detail {

  ActsMatrixD<8, 7>
  jacobianDirectionsToAngles(const Vector3D& dir)
  {
	ActsMatrixD<8, 7> jac = ActsMatrixD<8, 7>::Zero();
	
	const double x = dir(0);  // == cos(phi) * sin(theta)
    const double y = dir(1);  // == sin(phi) * sin(theta)
    const double z = dir(2);  // == cos(theta)
    // can be turned into cosine/sine
    const double cosTheta = z;
    const double sinTheta = sqrt(x * x + y * y);
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = x * invSinTheta;
    const double sinPhi = y * invSinTheta;
    
    jac(0, 0) = 1.;
    jac(1, 1) = 1.;
    jac(2, 2) = 1.;
    jac(3, 3) = 1.;
    jac(7, 6) = 1.;
    
	jac(4, 4) = -sinTheta * sinPhi;
    jac(4, 5) = cosTheta * cosPhi;
    jac(5, 4) = sinTheta * cosPhi;
    jac(5, 5) = cosTheta * sinPhi;
    jac(6, 5) = -sinTheta;
    return jac;
  }
  
ActsMatrixD<7, 8>
jacobianAnglesToDirections(const Vector3D& dir)
{
ActsMatrixD<7, 8> jacobian = ActsMatrixD<7, 8>::Zero();
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  const double z = dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  jacobian(0, 0) = 1.;
  jacobian(1, 1) = 1.;
  jacobian(2, 2) = 1.;
  jacobian(3, 3) = 1.;
  jacobian(6, 7) = 1.;
  
  jacobian(4, 4) = -sinPhi * invSinTheta;
  jacobian(4, 5) = cosPhi * invSinTheta;
  jacobian(5, 4) = cosPhi * cosTheta;
  jacobian(5, 5) = sinPhi * cosTheta;
  jacobian(5, 6) = -invSinTheta * (1. - cosTheta * cosTheta);

  return jacobian;
}

BoundState boundState(std::reference_wrapper<const GeometryContext> geoContext,
                      Covariance& covarianceMatrix, Jacobian& jacobian,
                      FreeMatrix& transportJacobian, FreeVector& derivatives,
                      std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal,
                      const ActsMatrixD<8, 7>& jacobianDirToAngle, const ActsMatrixD<7, 8>& jacobianAngleToDir,
                      const FreeVector& parameters, bool covTransport,
                      double accumulatedPath, const Surface& surface) {
  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    covarianceTransport(geoContext, covarianceMatrix, jacobian,
                        transportJacobian, derivatives, jacobianLocalToGlobal, jacobianDirToAngle, jacobianAngleToDir,
                        parameters, surface);
  }
  if (std::get<BoundSymMatrix>(covarianceMatrix) != BoundSymMatrix::Zero()) {
    cov = std::get<BoundSymMatrix>(covarianceMatrix);
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
                                  std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal,
                                  ActsMatrixD<8, 7>& jacobianDirToAngle, ActsMatrixD<7, 8>& jacobianAngleToDir,
                                  const FreeVector& parameters,
                                  bool covTransport, double accumulatedPath) {
  const Vector3D& direction = parameters.segment<3>(eFreeDir0);

  // Covariance transport
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (covTransport) {
    covarianceTransport(covarianceMatrix, jacobian, transportJacobian,
                        derivatives, jacobianLocalToGlobal, jacobianDirToAngle, jacobianAngleToDir,  direction, true);
  }
  if (std::get<BoundSymMatrix>(covarianceMatrix) != BoundSymMatrix::Zero()) {
    cov = std::get<BoundSymMatrix>(covarianceMatrix);
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

FreeState freeState(Covariance& covarianceMatrix, Jacobian& jacobian, FreeMatrix& transportJacobian,
                                  FreeVector& derivatives,
                                  std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal, ActsMatrixD<8, 7>& jacobianDirToAngle, ActsMatrixD<7, 8>& jacobianAngleToDir, const FreeVector& parameters, bool covTransport, double accumulatedPath) {
  // Transport the covariance to here
  std::optional<FreeSymMatrix> cov = std::nullopt;
  if (covTransport) {
    covarianceTransport(covarianceMatrix, jacobian, transportJacobian, derivatives, jacobianLocalToGlobal, jacobianDirToAngle, jacobianAngleToDir, parameters.segment<3>(4), false);
    cov = std::get<FreeSymMatrix>(covarianceMatrix);
  }
  // Create the free parameters
  FreeParameters freeParameters(std::move(cov), parameters);

  return std::make_tuple(std::move(freeParameters), jacobian,
                         accumulatedPath);
}
  
void covarianceTransport(Covariance& covarianceMatrix, Jacobian& jacobian,
                         FreeMatrix& transportJacobian, FreeVector& derivatives,
                         std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal,
                         ActsMatrixD<8, 7>& jacobianDirToAngle, ActsMatrixD<7, 8>& jacobianAngleToDir,
                         const Vector3D& direction, bool toLocal) {
  // Test if we started on a surface
if(jacobianLocalToGlobal.has_value())
{
	jacobianLocalToGlobal = transportJacobian * (*jacobianLocalToGlobal);
	
	// Test if we went to a surface
	if(toLocal)
	{
		const FreeToBoundMatrix jacToLocal =
			  surfaceDerivative(direction, *jacobianLocalToGlobal, derivatives);
		  const BoundMatrix jacFull = jacToLocal * (*jacobianLocalToGlobal);

		  // Apply the actual covariance transport
		  covarianceMatrix = BoundSymMatrix(jacFull * std::get<BoundSymMatrix>(covarianceMatrix) * jacFull.transpose());
		  
			// Store The global and bound jacobian (duplication for the moment)
		  jacobian = jacFull;
	}  
	else
	{
		covarianceMatrix = FreeSymMatrix((*jacobianLocalToGlobal) * std::get<BoundSymMatrix>(covarianceMatrix) * (*jacobianLocalToGlobal).transpose());
		jacobian = *jacobianLocalToGlobal;
	}
}
else
{
	if(toLocal)
	{
		const FreeToBoundMatrix jacToLocal =
			  surfaceDerivative(direction, jacobianDirToAngle, jacobianAngleToDir, transportJacobian, derivatives);
	  covarianceMatrix = BoundSymMatrix(jacToLocal * std::get<FreeSymMatrix>(covarianceMatrix) * jacToLocal.transpose());
	  jacobian = jacToLocal;
	}
	else
	{
		// Apply the actual covariance transport
		covarianceMatrix = FreeSymMatrix(transportJacobian * std::get<FreeSymMatrix>(covarianceMatrix) * transportJacobian.transpose());
		jacobian = transportJacobian;
	}
}

// Reinitialize jacobian components
if(toLocal)
  reinitializeJacobians(transportJacobian, derivatives, *jacobianLocalToGlobal,
                        direction);
	else
	reinitializeJacobians(transportJacobian, derivatives, jacobianLocalToGlobal, jacobianDirToAngle, jacobianAngleToDir,
                        direction);
}

void covarianceTransport(
    std::reference_wrapper<const GeometryContext> geoContext,
    Covariance& covarianceMatrix, Jacobian& jacobian,
    FreeMatrix& transportJacobian, FreeVector& derivatives,
    std::optional<BoundToFreeMatrix>& jacobianLocalToGlobal, const ActsMatrixD<8, 7>& jacobianDirToAngle, const ActsMatrixD<7, 8>& jacobianAngleToDir, const FreeVector& parameters,
    const Surface& surface) { 
      // Test if we started on a surface
	if(jacobianLocalToGlobal.has_value())
	{
		// Build the full jacobian
  jacobianLocalToGlobal = transportJacobian * (*jacobianLocalToGlobal);
  const FreeToBoundMatrix jacToLocal = surfaceDerivative(
      geoContext, parameters, *jacobianLocalToGlobal, derivatives, surface);
  const BoundMatrix jacFull = jacToLocal * (*jacobianLocalToGlobal);

  // Apply the actual covariance transport
  covarianceMatrix = BoundSymMatrix(jacFull * std::get<BoundSymMatrix>(covarianceMatrix) * jacFull.transpose());
  
  // Store The global and bound jacobian (duplication for the moment)
  jacobian = jacFull;
	}
	else
	{
		const FreeToBoundMatrix jacToLocal = surfaceDerivative(
			geoContext, parameters, jacobianDirToAngle, jacobianAngleToDir, transportJacobian, derivatives, surface);		
		// Apply the actual covariance transport
		covarianceMatrix = BoundSymMatrix(jacToLocal * std::get<FreeSymMatrix>(covarianceMatrix) * jacToLocal.transpose());
		jacobian = jacToLocal;
	}

  // Reinitialize jacobian components
  reinitializeJacobians(geoContext, transportJacobian, derivatives,
                        *jacobianLocalToGlobal, parameters, surface);
}
}  // namespace detail
}  // namespace Acts
