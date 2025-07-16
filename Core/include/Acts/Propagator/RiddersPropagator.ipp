// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/RiddersPropagator.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

#include <numbers>

namespace Acts {

namespace {

/// @brief This function tests whether the variations on a disc as target
/// surface lead to results on different sides wrt the center of the disc.
/// This would lead to a flip of the phi value on the surface and therewith
/// to a huge variance in that parameter. It can only occur in this
/// algorithm since the ridders algorithm is unaware of the target surface.
///
/// @param [in] derivatives Derivatives of a single parameter
///
/// @return Boolean result whether a phi jump occurred
inline bool inconsistentDerivativesOnDisc(
    const std::vector<BoundVector>& derivatives) {
  // Test each component with each other
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    bool jumpedAngle = true;
    for (unsigned int j = 0; j < derivatives.size(); j++) {
      // If there is at least one with a similar angle then it seems to work
      // properly
      if (i != j && std::abs(derivatives[i](1) - derivatives[j](1)) <
                        std::numbers::pi / 2.) {
        jumpedAngle = false;
        break;
      }
    }
    // Break if a jump was detected
    if (jumpedAngle) {
      return true;
    }
  }
  return false;
}

/// @brief This function fits a linear function through the final state
/// parametrisations
///
/// @param [in] deviations Vector of deviations
/// @param [in] derivatives Vector of resulting derivatives
///
/// @return Vector containing the linear fit
inline BoundVector fitLinear(const std::vector<double>& deviations,
                             const std::vector<BoundVector>& derivatives) {
  assert(!derivatives.empty() && "no fit without data");
  if (derivatives.size() == 1) {
    return derivatives[0];
  }

  BoundVector A = BoundVector::Zero();
  BoundVector C = BoundVector::Zero();
  double B = 0;
  double D = 0;
  const unsigned int N = deviations.size();

  for (unsigned int i = 0; i < N; ++i) {
    A += deviations.at(i) * derivatives.at(i);
    B += deviations.at(i);
    C += derivatives.at(i);
    D += deviations.at(i) * deviations.at(i);
  }

  BoundVector b = (N * A - B * C) / (N * D - B * B);
  BoundVector a = (C - B * b) / N;

  return a;
}

}  // namespace

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
auto RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const propagator_options_t& options) const
    -> Result<actor_list_t_result_t<
        BoundTrackParameters, typename propagator_options_t::actor_list_type>> {
  using ThisResult = Result<actor_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::actor_list_type>>;

  // Remove the covariance from our start parameters in order to skip jacobian
  // transport for the nominal propagation
  BoundTrackParameters startWithoutCov = start;
  startWithoutCov.covariance() = std::nullopt;

  // Propagate the nominal parameters
  auto result = m_propagator.propagate(startWithoutCov, options);
  if (!result.ok()) {
    return ThisResult::failure(result.error());
  }
  // Extract results from the nominal propagation
  auto nominalResult = result.value();
  assert(nominalResult.endParameters);
  const auto& nominalFinalParameters = *nominalResult.endParameters;

  BoundMatrix jacobian =
      wiggleAndCalculateJacobian(startWithoutCov, options, nominalResult);
  nominalResult.transportJacobian = jacobian;

  if (start.covariance()) {
    // use nominal parameters and Ridders covariance
    BoundMatrix cov = jacobian * (*start.covariance()) * jacobian.transpose();
    // replace the covariance of the nominal result w/ the ridders covariance
    nominalResult.endParameters = BoundTrackParameters::createCurvilinear(
        nominalFinalParameters.fourPosition(options.geoContext),
        nominalFinalParameters.direction(), nominalFinalParameters.qOverP(),
        std::move(cov), nominalFinalParameters.particleHypothesis());
  }

  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
auto RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<actor_list_t_result_t<
        BoundTrackParameters, typename propagator_options_t::actor_list_type>> {
  using ThisResult = Result<actor_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::actor_list_type>>;

  // Remove the covariance from our start parameters in order to skip jacobian
  // transport for the nominal propagation
  BoundTrackParameters startWithoutCov = start;
  startWithoutCov.covariance() = std::nullopt;

  // Propagate the nominal parameters
  auto result = m_propagator.propagate(startWithoutCov, target, options);
  if (!result.ok()) {
    return ThisResult::failure(result.error());
  }
  // Extract results from the nominal propagation
  auto nominalResult = result.value();
  assert(nominalResult.endParameters);
  const auto& nominalFinalParameters = *nominalResult.endParameters;

  BoundMatrix jacobian =
      wiggleAndCalculateJacobian(startWithoutCov, options, nominalResult);
  nominalResult.transportJacobian = jacobian;

  if (start.covariance()) {
    // use nominal parameters and Ridders covariance
    BoundMatrix cov = jacobian * (*start.covariance()) * jacobian.transpose();
    // replace the covariance of the nominal result w/ the ridders covariance
    nominalResult.endParameters = BoundTrackParameters(
        nominalFinalParameters.referenceSurface().getSharedPtr(),
        nominalFinalParameters.parameters(), std::move(cov),
        nominalFinalParameters.particleHypothesis());
  }

  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
BoundMatrix RiddersPropagator<propagator_t>::wiggleAndCalculateJacobian(
    const parameters_t& start, const propagator_options_t& options,
    const actor_list_t_result_t<BoundTrackParameters,
                                typename propagator_options_t::actor_list_type>&
        nominalResult) const {
  const auto& nominalFinalParameters = nominalResult.endParameters.value();
  // Use the curvilinear surface of the propagated parameters as target
  const Surface& target = nominalFinalParameters.referenceSurface();

  // Allow larger distances to reach the target surface
  propagator_options_t opts = options;
  opts.pathLimit *= 2;

  // Derivations of each parameter around the nominal parameters
  BoundMatrix jacobian = BoundMatrix::Zero();

  // Wiggle each dimension individually
  for (unsigned int i = 0; i < eBoundSize; i++) {
    // Steps for estimating derivatives
    std::vector<double> deviations = options.deviationFactors;
    for (double& d : deviations) {
      d *= options.deviationScale[i];
    }

    std::vector<BoundVector> derivatives =
        wiggleParameter(opts, start, i, target,
                        nominalFinalParameters.parameters(), deviations);

    // Test if target is disc - this may lead to inconsistent results
    if (target.type() == Surface::Disc) {
      assert(!inconsistentDerivativesOnDisc(derivatives));
    }

    jacobian.col(i) = fitLinear(deviations, derivatives);
  }

  return jacobian;
}

template <typename propagator_t>
template <typename propagator_options_t, typename parameters_t>
std::vector<BoundVector> RiddersPropagator<propagator_t>::wiggleParameter(
    const propagator_options_t& options, const parameters_t& start,
    const unsigned int param, const Surface& target, const BoundVector& nominal,
    const std::vector<double>& deviations) const {
  // Storage of the results
  std::vector<BoundVector> derivatives;
  derivatives.reserve(deviations.size());
  for (double h : deviations) {
    // Treatment for theta
    if (param == eBoundTheta) {
      const double current_theta = start.template get<eBoundTheta>();
      if (current_theta + h > std::numbers::pi) {
        h = std::numbers::pi - current_theta;
      }
      if (current_theta + h < 0) {
        h = -current_theta;
      }
    }

    // Modify start parameter
    BoundVector values = start.parameters();
    values[param] += h;

    // Propagate with updated start parameters
    BoundTrackParameters tp(start.referenceSurface().getSharedPtr(), values,
                            start.covariance(), start.particleHypothesis());
    const auto& r = m_propagator.propagate(tp, target, options).value();
    // Collect the slope
    derivatives.push_back((r.endParameters->parameters() - nominal) / h);

    // Correct for a possible variation of phi around
    if (param == eBoundPhi) {
      double phi0 = nominal(eBoundPhi);
      double phi1 = r.endParameters->parameters()(eBoundPhi);
      if (std::abs(phi1 + 2. * std::numbers::pi - phi0) <
          std::abs(phi1 - phi0)) {
        derivatives.back()[eBoundPhi] =
            (phi1 + 2. * std::numbers::pi - phi0) / h;
      } else if (std::abs(phi1 - 2. * std::numbers::pi - phi0) <
                 std::abs(phi1 - phi0)) {
        derivatives.back()[eBoundPhi] =
            (phi1 - 2. * std::numbers::pi - phi0) / h;
      }
    }
  }

  return derivatives;
}

}  // namespace Acts
