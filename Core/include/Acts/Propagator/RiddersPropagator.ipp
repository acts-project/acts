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
  const BoundVector& nominalParameters =
      nominalResult.endParameters->parameters();
  // Use the curvilinear surface of the propagated parameters as target
  const Surface& surface = nominalResult.endParameters->referenceSurface();

  // Steps for estimating derivatives
  std::vector<double> deviations = {-4e-4, -2e-4, 2e-4, 4e-4};

  // Allow larger distances for the oscillation
  propagator_options_t opts = options;
  opts.pathLimit *= 2.;

  // Derivations of each parameter around the nominal parameters
  std::array<std::vector<BoundVector>, eBoundSize> derivatives;

  // Wiggle each dimension individually
  for (unsigned int i = 0; i < eBoundSize; i++) {
    derivatives[i] =
        wiggleDimension(opts, start, i, surface, nominalParameters, deviations);
  }
  if (start.covariance()) {
    auto cov =
        calculateCovariance(derivatives, *start.covariance(), deviations);
    // replace the covariance of the nominal result w/ the ridders covariance
    auto& nom = *nominalResult.endParameters;
    nom = CurvilinearTrackParameters(
        nom.fourPosition(options.geoContext), nom.unitDirection(),
        nom.absoluteMomentum(), nom.charge(), std::move(cov));
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
  const BoundVector& nominalParameters =
      nominalResult.endParameters->parameters();

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
    derivatives[i] =
        wiggleDimension(opts, start, i, target, nominalParameters, deviations);
  }
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
    auto cov =
        calculateCovariance(derivatives, *start.covariance(), deviations);
    auto& nom = *nominalResult.endParameters;
    nom = BoundTrackParameters(nom.referenceSurface().getSharedPtr(),
                               nom.parameters(), nom.charge(), std::move(cov));
  }
  return ThisResult::success(std::move(nominalResult));
}

template <typename propagator_t>
bool Acts::RiddersPropagator<propagator_t>::inconsistentDerivativesOnDisc(
    const std::vector<Acts::BoundVector>& derivatives) const {
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
Acts::RiddersPropagator<propagator_t>::wiggleDimension(
    const options_t& options, const parameters_t& startPars,
    const unsigned int param, const Surface& target,
    const Acts::BoundVector& nominal,
    const std::vector<double>& deviations) const {
  // Storage of the results
  std::vector<BoundVector> derivatives;
  derivatives.reserve(deviations.size());
  for (double h : deviations) {
    // Treatment for theta
    if (param == eBoundTheta) {
      const double current_theta = startPars.template get<eBoundTheta>();
      if (current_theta + h > M_PI) {
        h = M_PI - current_theta;
      }
      if (current_theta + h < 0) {
        h = -current_theta;
      }
    }

    // Modify start parameter
    BoundVector values = startPars.parameters();
    values[param] += h;

    // Propagate with updated start parameters
    BoundTrackParameters tp(startPars.referenceSurface().getSharedPtr(), values,
                            startPars.covariance());
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
auto Acts::RiddersPropagator<propagator_t>::calculateCovariance(
    const std::array<std::vector<Acts::BoundVector>, Acts::eBoundSize>&
        derivatives,
    const Acts::BoundSymMatrix& startCov,
    const std::vector<double>& deviations) const -> Covariance {
  Jacobian jacobian;
  jacobian.setIdentity();
  jacobian.col(eBoundLoc0) = fitLinear(derivatives[eBoundLoc0], deviations);
  jacobian.col(eBoundLoc1) = fitLinear(derivatives[eBoundLoc1], deviations);
  jacobian.col(eBoundPhi) = fitLinear(derivatives[eBoundPhi], deviations);
  jacobian.col(eBoundTheta) = fitLinear(derivatives[eBoundTheta], deviations);
  jacobian.col(eBoundQOverP) = fitLinear(derivatives[eBoundQOverP], deviations);
  jacobian.col(eBoundTime) = fitLinear(derivatives[eBoundTime], deviations);
  return jacobian * startCov * jacobian.transpose();
}

template <typename propagator_t>
Acts::BoundVector Acts::RiddersPropagator<propagator_t>::fitLinear(
    const std::vector<Acts::BoundVector>& values,
    const std::vector<double>& deviations) const {
  BoundVector A;
  BoundVector C;
  A.setZero();
  C.setZero();
  double B = 0;
  double D = 0;
  const unsigned int N = deviations.size();

  for (unsigned int i = 0; i < N; ++i) {
    A += deviations.at(i) * values.at(i);
    B += deviations.at(i);
    C += values.at(i);
    D += deviations.at(i) * deviations.at(i);
  }

  BoundVector b = (N * A - B * C) / (N * D - B * B);
  BoundVector a = (C - B * b) / N;

  return a;
}
