// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
auto Acts::RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const propagator_options_t& options) const
    -> Result<action_list_t_result_t<
        CurvilinearTrackParameters,
        typename propagator_options_t::action_list_type>> {
  using ThisResult = Result<
      action_list_t_result_t<CurvilinearTrackParameters,
                             typename propagator_options_t::action_list_type>>;

  // Propagate the nominal parameters
  auto nominalRet = m_propagator.propagate(start, options);
  if (not nominalRet.ok()) {
    return ThisResult::failure(nominalRet.error());
  }

  // Extract results from the nominal propagation
  auto nominalResult = std::move(nominalRet).value();
  assert(nominalResult.endParameters);
  auto nominalFinalParameters = *nominalResult.endParameters;
  const BoundVector& nominalFinalParameterValues =
      nominalFinalParameters.parameters();
  // Use the curvilinear surface of the propagated parameters as target
  const Surface& surface = nominalFinalParameters.referenceSurface();

  // Steps for estimating derivatives
  std::vector<double> deviations = {-4e-4, -2e-4, 2e-4, 4e-4};

  // Allow larger distances for the oscillation
  propagator_options_t opts = options;
  opts.pathLimit *= 2.;

  // Derivations of each parameter around the nominal parameters
  std::array<std::vector<BoundVector>, eBoundSize> derivatives;

  // Wiggle each dimension individually
  for (unsigned int i = 0; i < eBoundSize; i++) {
    derivatives[i] = wiggleParameter(opts, start, i, surface,
                                     nominalFinalParameterValues, deviations);
  }

  auto jacobian = calculateJacobian(deviations, derivatives);
  nominalResult.transportJacobian = jacobian;

  if (start.covariance()) {
    auto cov = jacobian * (*start.covariance()) * jacobian.transpose();
    // replace the covariance of the nominal result w/ the ridders covariance
    nominalResult.endParameters = CurvilinearTrackParameters(
        nominalFinalParameters.fourPosition(options.geoContext),
        nominalFinalParameters.unitDirection(),
        nominalFinalParameters.absoluteMomentum(),
        nominalFinalParameters.charge(), std::move(cov));
  }

  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
template <typename parameters_t, typename propagator_options_t>
auto Acts::RiddersPropagator<propagator_t>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<action_list_t_result_t<
        BoundTrackParameters,
        typename propagator_options_t::action_list_type>> {
  using ThisResult = Result<action_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::action_list_type>>;

  // Propagate the nominal parameters
  auto nominalRet = m_propagator.propagate(start, target, options);
  if (not nominalRet.ok()) {
    return ThisResult::failure(nominalRet.error());
  }

  // Extract results from the nominal propagation
  auto nominalResult = std::move(nominalRet).value();
  assert(nominalResult.endParameters);
  auto nominalFinalParameters = *nominalResult.endParameters;
  const BoundVector& nominalFinalParameterValues =
      nominalFinalParameters.parameters();
  // Use the curvilinear surface of the propagated parameters as target
  const Surface& surface = nominalFinalParameters.referenceSurface();

  // Steps for estimating derivatives
  std::vector<double> deviations = {-4e-4, -2e-4, 2e-4, 4e-4};
  if (target.type() == Surface::Disc) {
    deviations = {{-3e-5, -1e-5, 1e-5, 3e-5}};
  }

  // - for planar surfaces the dest surface is a perfect destination
  // surface for the numerical propagation, as reference frame
  // aligns with the referenceSurface.transform().rotation() at
  // at any given time
  //
  // - for straw & cylinder, where the error is given
  // in the reference frame that re-aligns with a slightly different
  // intersection solution

  // Allow larger distances for the oscillation
  propagator_options_t opts = options;
  opts.pathLimit *= 2.;

  // Derivations of each parameter around the nominal parameters
  std::array<std::vector<BoundVector>, eBoundSize> derivatives;

  // Wiggle each dimension individually
  for (unsigned int i = 0; i < eBoundSize; i++) {
    derivatives[i] = wiggleParameter(opts, start, i, surface,
                                     nominalFinalParameterValues, deviations);
  }

  auto jacobian = calculateJacobian(deviations, derivatives);
  nominalResult.transportJacobian = jacobian;

  if (start.covariance()) {
    // Test if target is disc - this may lead to inconsistent results
    if (target.type() == Surface::Disc) {
      for (const std::vector<BoundVector>& deriv : derivatives) {
        if (inconsistentDerivativesOnDisc(deriv)) {
          // Set covariance to zero and return
          // TODO: This should be changed to indicate that something went
          // wrong
          return ThisResult::success(std::move(nominalResult));
        }
      }
    }
    // use nominal parameters and Ridders covariance
    auto cov = jacobian * (*start.covariance()) * jacobian.transpose();
    // replace the covariance of the nominal result w/ the ridders covariance
    nominalResult.endParameters = CurvilinearTrackParameters(
        nominalFinalParameters.fourPosition(options.geoContext),
        nominalFinalParameters.unitDirection(),
        nominalFinalParameters.absoluteMomentum(),
        nominalFinalParameters.charge(), std::move(cov));
  }

  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
bool Acts::RiddersPropagator<propagator_t>::inconsistentDerivativesOnDisc(
    const std::vector<Acts::BoundVector>& derivatives) {
  // Test each component with each other
  for (unsigned int i = 0; i < derivatives.size(); i++) {
    bool jumpedAngle = true;
    for (unsigned int j = 0; j < derivatives.size(); j++) {
      // If there is at least one with a similar angle then it seems to work
      // properly
      if (i != j &&
          std::abs(derivatives[i](1) - derivatives[j](1)) < 0.5 * M_PI) {
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

template <typename propagator_t>
template <typename options_t, typename parameters_t>
std::vector<Acts::BoundVector>
Acts::RiddersPropagator<propagator_t>::wiggleParameter(
    const options_t& options, const parameters_t& start,
    const unsigned int param, const Surface& target,
    const Acts::BoundVector& nominal,
    const std::vector<double>& deviations) const {
  // Storage of the results
  std::vector<BoundVector> derivatives;
  derivatives.reserve(deviations.size());
  for (double h : deviations) {
    // Treatment for theta
    if (param == eBoundTheta) {
      const double current_theta = start.template get<eBoundTheta>();
      if (current_theta + h > M_PI) {
        h = M_PI - current_theta;
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
                            start.covariance());
    const auto& r = m_propagator.propagate(tp, target, options).value();
    // Collect the slope
    derivatives.push_back((r.endParameters->parameters() - nominal) / h);

    // Correct for a possible variation of phi around
    if (param == eBoundPhi) {
      double phi0 = nominal(Acts::eBoundPhi);
      double phi1 = r.endParameters->parameters()(Acts::eBoundPhi);
      if (std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0)) {
        derivatives.back()[Acts::eBoundPhi] = (phi1 + 2. * M_PI - phi0) / h;
      } else if (std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0)) {
        derivatives.back()[Acts::eBoundPhi] = (phi1 - 2. * M_PI - phi0) / h;
      }
    }
  }

  return derivatives;
}

template <typename propagator_t>
typename Acts::RiddersPropagator<propagator_t>::Jacobian
Acts::RiddersPropagator<propagator_t>::calculateJacobian(
    const std::vector<double>& deviations,
    const std::array<std::vector<Acts::BoundVector>, Acts::eBoundSize>&
        derivatives) {
  Jacobian jacobian = Jacobian::Identity();
  jacobian.col(eBoundLoc0) = fitLinear(deviations, derivatives[eBoundLoc0]);
  jacobian.col(eBoundLoc1) = fitLinear(deviations, derivatives[eBoundLoc1]);
  jacobian.col(eBoundPhi) = fitLinear(deviations, derivatives[eBoundPhi]);
  jacobian.col(eBoundTheta) = fitLinear(deviations, derivatives[eBoundTheta]);
  jacobian.col(eBoundQOverP) = fitLinear(deviations, derivatives[eBoundQOverP]);
  jacobian.col(eBoundTime) = fitLinear(deviations, derivatives[eBoundTime]);
  return jacobian;
}

template <typename propagator_t>
Acts::BoundVector Acts::RiddersPropagator<propagator_t>::fitLinear(
    const std::vector<double>& deviations,
    const std::vector<Acts::BoundVector>& derivatives) {
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
