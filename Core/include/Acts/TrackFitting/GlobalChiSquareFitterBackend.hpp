// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Non-template part of the GX2F: equation system, scattering properties and
// the backend functions. Kept apart from GlobalChiSquareFitter.hpp so the
// backend translation units do not include the full fitter.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace Acts::Experimental {

/// @addtogroup track_fitting
/// @{

namespace Gx2fConstants {
constexpr std::string_view gx2fnUpdateColumn = "Gx2fnUpdateColumn";

// Mask for the track states. We don't need Predicted and Filtered
constexpr TrackStatePropMask trackStateMask = TrackStatePropMask::Smoothed |
                                              TrackStatePropMask::Jacobian |
                                              TrackStatePropMask::Calibrated;

// A projector used for scattering. By using Jacobian * phiThetaProjector one
// gets only the derivatives for the variables phi and theta.
const Eigen::Matrix<double, eBoundSize, 2> phiThetaProjector = [] {
  Eigen::Matrix<double, eBoundSize, 2> m =
      Eigen::Matrix<double, eBoundSize, 2>::Zero();
  m(eBoundPhi, 0) = 1.0;
  m(eBoundTheta, 1) = 1.0;
  return m;
}();
}  // namespace Gx2fConstants

/// @brief A container to store scattering properties for each material surface
///
/// This struct holds the scattering angles, the inverse covariance of the
/// material, and a validity flag indicating whether the material is valid for
/// the scattering process.
struct ScatteringProperties {
 public:
  /// @brief Constructor to initialize scattering properties.
  ///
  /// @param scatteringAngles_ The vector of scattering angles.
  /// @param invCovarianceMaterial_ The inverse covariance of the material.
  /// @param materialIsValid_ A boolean flag indicating whether the material is valid.
  ScatteringProperties(const BoundVector& scatteringAngles_,
                       const double invCovarianceMaterial_,
                       const bool materialIsValid_)
      : m_scatteringAngles(scatteringAngles_),
        m_invCovarianceMaterial(invCovarianceMaterial_),
        m_materialIsValid(materialIsValid_) {}

  /// @brief Accessor for the scattering angles (const version)
  /// @return Const reference to the vector of scattering angles
  const BoundVector& scatteringAngles() const { return m_scatteringAngles; }

  /// @brief Accessor for the scattering angles (mutable version)
  /// @return Mutable reference to the vector of scattering angles for modification
  BoundVector& scatteringAngles() { return m_scatteringAngles; }

  /// @brief Accessor for the inverse covariance of the material
  /// @return Inverse covariance value computed from material properties (e.g., Highland formula)
  double invCovarianceMaterial() const { return m_invCovarianceMaterial; }

  /// @brief Accessor for the material validity flag
  /// @return True if material is valid for scattering calculations, false for vacuum or zero thickness
  bool materialIsValid() const { return m_materialIsValid; }

 private:
  /// Vector of scattering angles. The vector is usually all zeros except for
  /// eBoundPhi and eBoundTheta.
  BoundVector m_scatteringAngles;

  /// Inverse covariance of the material. Compute with e.g. the Highland
  /// formula.
  double m_invCovarianceMaterial;

  /// Flag indicating whether the material is valid. Commonly vacuum and zero
  /// thickness material will be ignored.
  bool m_materialIsValid;
};

/// @brief A container to manage all properties of a gx2f system
///
/// This struct manages the mathematical infrastructure for the gx2f. It
/// initializes and maintains the extended aMatrix and extended bVector.
struct Gx2fSystem {
 public:
  /// @brief Constructor to initialize matrices and vectors to zero based on specified dimensions.
  ///
  /// @param nDims Number of dimensions for the extended matrix and vector.
  explicit Gx2fSystem(std::size_t nDims)
      : m_nDims{nDims},
        m_aMatrix{Eigen::MatrixXd::Zero(nDims, nDims)},
        m_bVector{Eigen::VectorXd::Zero(nDims)} {}

  /// @brief Accessor for the number of dimensions of the extended system
  /// @return Number of dimensions for the aMatrix and bVector (bound parameters + scattering angles)
  std::size_t nDims() const { return m_nDims; }

  /// @brief Accessor for the accumulated chi-squared value (const version)
  /// @return Current sum of chi-squared contributions from measurements and material
  double chi2() const { return m_chi2; }

  /// @brief Accessor for the accumulated chi-squared value (mutable version)
  /// @return Mutable reference to chi-squared sum for modification during fitting
  double& chi2() { return m_chi2; }

  /// @brief Accessor for the extended system matrix (const version)
  /// @return Const reference to the aMatrix containing measurement and material contributions
  const Eigen::MatrixXd& aMatrix() const { return m_aMatrix; }

  /// @brief Accessor for the extended system matrix (mutable version)
  /// @return Mutable reference to the aMatrix for adding measurement and material contributions
  Eigen::MatrixXd& aMatrix() { return m_aMatrix; }

  /// @brief Accessor for the extended system vector (const version)
  /// @return Const reference to the bVector containing measurement and material contributions
  const Eigen::VectorXd& bVector() const { return m_bVector; }

  /// @brief Accessor for the extended system vector (mutable version)
  /// @return Mutable reference to the bVector for adding measurement and material contributions
  Eigen::VectorXd& bVector() { return m_bVector; }

  /// @brief Accessor for the number of degrees of freedom (const version)
  /// @return Current number of degrees of freedom from processed measurements
  std::size_t ndf() const { return m_ndf; }

  /// @brief Accessor for the number of degrees of freedom (mutable version)
  /// @return Mutable reference to NDF counter for incrementing during measurement processing
  std::size_t& ndf() { return m_ndf; }

  /// @brief Determines the minimum number of degrees of freedom required for the fit
  ///
  /// Automatically deduces the required NDF based on the system configuration.
  /// We have only 3 cases, because we always have l0, l1, phi, theta:
  /// - 4: no magnetic field -> q/p is empty
  /// - 5: no time measurement -> time is not fittable
  /// - 6: full fit with all parameters
  ///
  /// @return Required NDF based on which parameters can be fitted
  std::size_t findRequiredNdf() {
    std::size_t ndfSystem = 0;
    if (m_aMatrix(4, 4) == 0) {
      ndfSystem = 4;
    } else if (m_aMatrix(5, 5) == 0) {
      ndfSystem = 5;
    } else {
      ndfSystem = 6;
    }

    return ndfSystem;
  }

  /// @brief Checks if the system has sufficient degrees of freedom for fitting
  /// @return True if NDF exceeds the minimum required for the parameter configuration
  bool isWellDefined() { return m_ndf > findRequiredNdf(); }

 private:
  /// Number of dimensions of the (extended) system
  std::size_t m_nDims;

  /// Sum of chi-squared values.
  double m_chi2 = 0.;

  /// Extended matrix for accumulation.
  Eigen::MatrixXd m_aMatrix;

  /// Extended vector for accumulation.
  Eigen::VectorXd m_bVector;

  /// Number of degrees of freedom of the system
  std::size_t m_ndf = 0u;
};

/// @brief Adds a measurement to the GX2F equation system in a modular backend function.
///
/// This function processes measurement data and integrates it into the GX2F
/// system.
///
/// @param extendedSystem All parameters of the current equation system to update.
/// @param jacobianFromStart The Jacobian matrix from the start to the current state.
/// @param covarianceMeasurement The covariance matrix of the measurement.
/// @param predicted The predicted state vector based on the track state.
/// @param measurement The measurement vector.
/// @param projector The projection matrix.
/// @param logger A logger instance.
///
/// @note The dynamic Eigen matrices are suboptimal. We could think of
/// templating again in the future on kMeasDims. We currently use dynamic
/// matrices to reduce the memory during compile time.
void addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger);

/// @brief Solve the gx2f system to get the delta parameters for the update
///
/// This function computes the delta parameters for the GX2F Actor fitting
/// process by solving the linear equation system [a] * delta = b. It uses the
/// column-pivoting Householder QR decomposition for numerical stability.
///
/// @param extendedSystem All parameters of the current equation system
/// @return Delta parameters for the GX2F update
Eigen::VectorXd computeGx2fDeltaParams(const Gx2fSystem& extendedSystem);

/// @brief Update parameters (and scattering angles if applicable)
///
/// @param params Parameters to be updated
/// @param deltaParamsExtended Delta parameters for bound parameter and scattering angles
/// @param nMaterialSurfaces Number of material surfaces in the track
/// @param scatteringMap Map of geometry identifiers to scattering properties,
///        containing all scattering angles and covariances
/// @param geoIdVector Vector of geometry identifiers corresponding to material surfaces
void updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    const std::size_t nMaterialSurfaces,
    std::unordered_map<GeometryIdentifier, ScatteringProperties>& scatteringMap,
    const std::vector<GeometryIdentifier>& geoIdVector);

/// @brief Calculate and update the covariance of the fitted parameters
///
/// This function calculates the covariance of the fitted parameters using
/// cov = inv([a])
/// It then updates the first square block of size ndfSystem. This ensures,
/// that we only update the covariance for fitted parameters. (In case of
/// no qop/time fit)
///
/// @param fullCovariancePredicted The covariance matrix to update
/// @param extendedSystem All parameters of the current equation system
void updateGx2fCovarianceParams(BoundMatrix& fullCovariancePredicted,
                                Gx2fSystem& extendedSystem);

/// @}

}  // namespace Acts::Experimental
