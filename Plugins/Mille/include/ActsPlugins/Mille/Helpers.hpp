// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
namespace ActsPlugins::ActsToMille {

/// @brief: Regularise a covariance matrix for decomposition into Mille.
///         Required especially when running fits with non-timing detectors.
/// @param inputCov: Input covariance matrix, to be regularised
/// @param conditionCutOff: Lowest value (relative to leading EV after clamping)
///                         to clamp the eigenvalues to.
/// @param removeHugeLeading: If set to a positive value, the leading eigenvalue
///                           will be regularised to the value of the
///                           second-leading if the former exceeds the latter by
///                           more than this factor. Use to suppress poorly
///                           constrained directions (e.g. time coordinate).
/// @param stabilisationDiag: If set to a positive value, the main diagonal elements
///                           of the input matrix will be incremented by this
///                           value to regularise the problem.
/// @param return A new matrix, which has been regularised by clamping
///               eigenvalues to the iterval [conditionCutOff.first x max_EV,
///               max_EV], where max_EV is either the leading eigenvalue or, if
///               removeHugeLeading is positive, either the first or the second
///               leading EV (second if largest is more than removeHugeLeading
///               times the second)
Acts::DynamicMatrix regulariseCovariance(const Acts::DynamicMatrix& inputCov,
                                         double conditionCutOff = 1e-10,
                                         double removeHugeLeading = 100.,
                                         double stabilisationDiag = 1.e-10);

/// Calculates the solution X to the matrix equation C = (A + X)^-1,
/// for a poorly conditioned C and known A. Required to decompose the ACTS
/// Kalman covariance matrix into a series of Mille pseudo-measurements.
/// Uses cholesky factorisation of C to solve (1 - CA) = CX
/// to avoid directly inverting C.
/// @param target: The target matrix C to be decomposed - expected to be symmetric
///                and positive (semi)definite
/// @param existing_sol: The partial existing solution A - expected to be symmetric
/// @return the missing piece X solving the above equation.
Acts::DynamicMatrix getInverseComplement(
    const Acts::DynamicMatrix& target, const Acts::DynamicMatrix& existing_sol);

}  // namespace ActsPlugins::ActsToMille
