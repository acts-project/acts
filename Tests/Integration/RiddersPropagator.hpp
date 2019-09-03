// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <optional>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

namespace IntegrationTest {

using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;

/// @brief This class performs the Ridders algorithm to estimate the propagation
/// of the covariance to a certain point in space.
///
/// The algorithm is based on the small deviations of the start parameters based
/// on their uncertainty at the beginning of the propgation. This deviation is
/// represented here by a vector of relative deviations of these parameters and
/// fix for all parameters. So, a common choice has to be found that is able to
/// actually fit into the order of magnitude of the uncertainty of each
/// parameter. Using these deviations, the propagation is repeated multiple
/// times and the final covariance matrix at a given target surface is
/// afterwards evaluated by first order derivatives of the final state
/// parameters wrt. the inital parameters. Therefore this evaluation represents
/// a first order approximation of the transport jacobian. Since performing
/// multiple propagations and a numerical evaluation of the covariance requires
/// more time than a single propagation towards a target + a common propagation
/// of the covariance, this class just serves to verify the results of the
/// latter classes.
template <typename propagator_t>
class RiddersPropagator {
 private:
  ///
  /// @note The result_type_helper struct and the action_list_t_result_t are
  /// here to allow a look'n'feel of this class like the Propagator itself
  ///

  /// @copydoc Propagator::result_type_helper
  template <typename parameters_t, typename action_list_t>
  struct result_type_helper {
    /// @copydoc Propagator::result_type_helper::this_result_type
    template <typename... args>
    using this_result_type = PropagatorResult<parameters_t, args...>;

    /// @copydoc Propagator::result_type_helper::type
    using type = typename action_list_t::template result_type<this_result_type>;
  };

  /// @copydoc Propagator::action_list_t_result_t
  template <typename parameters_t, typename action_list_t>
  using action_list_t_result_t =
      typename result_type_helper<parameters_t, action_list_t>::type;

  /// @brief This function tests whether the variations on a disc as target
  /// surface lead to results on different sides wrt the center of the disc.
  /// This would lead to a flip of the phi value on the surface and therewith to
  /// a huge variance in that parameter. It can only occur in this algorithm
  /// since the ridders algorithm is unaware of the target surface.
  ///
  /// @param [in] derivatives Derivatives of a single parameter
  ///
  /// @return Boolean result whether a phi jump occured
  bool inconsistentDerivativesOnDisc(
      const std::vector<BoundVector>& derivatives) const {
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

 public:
  /// @brief Constructor using a propagator
  ///
  /// @param [in] propagator Underlying propagator that will be used
  RiddersPropagator(propagator_t& propagator) : m_propagator(propagator) {}

  /// @brief Constructor building a propagator
  ///
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in] stepper Stepper that will be used
  /// @param [in] navigator Navigator that will be used
  template <typename stepper_t, typename navigator_t = detail::VoidNavigator>
  RiddersPropagator(stepper_t stepper, navigator_t navigator = navigator_t())
      : m_propagator(Propagator(stepper, navigator)) {}

  /// @brief Propagation method targeting curvilinear parameters
  ///
  /// @tparam parameters_t Type of the start parameters
  /// @tparam action_list_t Type of the action list
  /// @tparam aborter_list_t Type of the aborter list
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  template <typename parameters_t, typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t>
  Result<
      action_list_t_result_t<typename propagator_t::Stepper::
                                 template return_parameter_type<parameters_t>,
                             action_list_t>>
  propagate(const parameters_t& start,
            const propagator_options_t<action_list_t, aborter_list_t>& options)
      const {
    // Launch nominal propagation and collect results
    auto nominalResult = m_propagator.propagate(start, options).value();
    const BoundVector& nominalParameters =
        nominalResult.endParameters->parameters();
    // Pick the surface of the propagation as target
    const Surface& surface = nominalResult.endParameters->referenceSurface();

    // Steps for estimating derivatives
    std::vector<double> deviations = {-2e-4, -1e-4, 1e-4, 2e-4};

    // Allow larger distances for the oscillation
    propagator_options_t<action_list_t, aborter_list_t> opts = options;
    opts.pathLimit *= 2.;

    // Derivations of each parameter around the nominal parameters
    std::array<std::vector<BoundVector>, BoundParsDim> derivatives;

    // Wiggle each dimension individually
    for (unsigned int i = 0; i < BoundParsDim; i++) {
      derivatives[i] = wiggleDimension(opts, start, i, surface,
                                       nominalParameters, deviations);
    }
    // Exchange the result by Ridders Covariance
    const FullParameterSet& parSet =
        nominalResult.endParameters->getParameterSet();
    FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
    if (start.covariance()) {
      mParSet->setCovariance(
          calculateCovariance(derivatives, *start.covariance(), deviations));
    }

    return std::move(nominalResult);
  }

  /// @brief Propagation method targeting bound parameters
  ///
  /// @tparam parameters_t Type of the start parameters
  /// @tparam surface_t Type of target surface
  /// @tparam action_list_t Type of the action list
  /// @tparam aborter_list_t Type of the aborter list
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  /// @note If the target surface is a disc, the resulting covariance may be
  /// inconsistent. In this case a zero matrix is returned.
  template <typename parameters_t, typename surface_t, typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t>
  Result<action_list_t_result_t<
      typename propagator_t::Stepper::template return_parameter_type<
          parameters_t, surface_t>,
      action_list_t>>
  propagate(const parameters_t& start, const surface_t& target,
            const propagator_options_t<action_list_t, aborter_list_t>& options)
      const {
    // Launch nominal propagation and collect results
    //~ Result<action_list_t_result_t<typename propagator_t::Stepper::template
    //return_parameter_type<parameters_t>, action_list_t>> nominalResult =
    //m_propagator.propagate(start, target, options);
    auto nominalResult = m_propagator.propagate(start, target, options).value();
    const BoundVector& nominalParameters =
        nominalResult.endParameters->parameters();

    // Steps for estimating derivatives
    std::vector<double> deviations = {-2e-4, -1e-4, 1e-4, 2e-4};
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
    propagator_options_t<action_list_t, aborter_list_t> opts = options;
    opts.pathLimit *= 2.;

    // Derivations of each parameter around the nominal parameters
    std::array<std::vector<BoundVector>, BoundParsDim> derivatives;

    // Wiggle each dimension individually
    for (unsigned int i = 0; i < BoundParsDim; i++) {
      derivatives[i] = wiggleDimension(opts, start, i, target,
                                       nominalParameters, deviations);
    }
    // Exchange the result by Ridders Covariance
    const FullParameterSet& parSet =
        nominalResult.endParameters->getParameterSet();
    FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
    if (start.covariance()) {
      // Test if target is disc - this may lead to inconsistent results
      if (target.type() == Surface::Disc) {
        for (const std::vector<BoundVector>& deriv : derivatives) {
          if (inconsistentDerivativesOnDisc(deriv)) {
            // Set covariance to zero and return
            // TODO: This should be changed to indicate that something went
            // wrong
            mParSet->setCovariance(Covariance::Zero());
            return std::move(nominalResult);
          }
        }
      }
      mParSet->setCovariance(
          calculateCovariance(derivatives, *start.covariance(), deviations));
    }
    return std::move(nominalResult);
  }

 private:
  /// @brief This function wiggles one dimension of the starting parameters,
  /// performs the propagation to a surface and collects for each change of the
  /// start parameters the slope
  ///
  /// @tparam options_t PropagatorOptions object
  /// @tparam parameters+t Type of the parameters to start the propagation with
  ///
  /// @param [in] options Options do define how to wiggle
  /// @param [in] startPart Start parameters that are modified
  /// @param [in] param Index to get the parameter that will be modified
  /// @param [in] target Target surface
  /// @param [in] nominal Nominal end parameters
  ///
  /// @return Vector containing each slope
  template <typename options_t, typename parameters_t>
  std::vector<BoundVector> wiggleDimension(
      const options_t& options, const parameters_t& startPars,
      const unsigned int param, const Surface& target,
      const BoundVector& nominal, const std::vector<double>& deviations) const {
    // Storage of the results
    std::vector<BoundVector> derivatives;
    derivatives.reserve(deviations.size());
    for (double h : deviations) {
      parameters_t tp = startPars;

      // Treatment for theta
      if (param == eTHETA) {
        const double current_theta = tp.template get<eTHETA>();
        if (current_theta + h > M_PI) {
          h = M_PI - current_theta;
        }
        if (current_theta + h < 0) {
          h = -current_theta;
        }
      }

      // Modify start parameter and propagate
      switch (param) {
        case 0: {
          tp.template set<eLOC_0>(options.geoContext,
                                  tp.template get<eLOC_0>() + h);
          break;
        }
        case 1: {
          tp.template set<eLOC_1>(options.geoContext,
                                  tp.template get<eLOC_1>() + h);
          break;
        }
        case 2: {
          tp.template set<ePHI>(options.geoContext,
                                tp.template get<ePHI>() + h);
          break;
        }
        case 3: {
          tp.template set<eTHETA>(options.geoContext,
                                  tp.template get<eTHETA>() + h);
          break;
        }
        case 4: {
          tp.template set<eQOP>(options.geoContext,
                                tp.template get<eQOP>() + h);
          break;
        }
        case 5: {
          tp.template set<eT>(options.geoContext, tp.template get<eT>() + h);
          break;
        }
        default:
          return {};
      }
      const auto& r = m_propagator.propagate(tp, target, options).value();
      // Collect the slope
      derivatives.push_back((r.endParameters->parameters() - nominal) / h);

      // Correct for a possible variation of phi around
      if (param == 2) {
        double phi0 = nominal(Acts::ePHI);
        double phi1 = r.endParameters->parameters()(Acts::ePHI);
        if (std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 + 2. * M_PI - phi0) / h;
        else if (std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0))
          derivatives.back()[Acts::ePHI] = (phi1 - 2. * M_PI - phi0) / h;
      }
    }
    return derivatives;
  }

  /// @brief This function propagates the covariance matrix
  ///
  /// @param [in] derivatives Slopes of each modification of the parameters
  /// @param [in] startCov Starting covariance
  ///
  /// @return Propagated covariance matrix
  const Covariance calculateCovariance(
      const std::array<std::vector<BoundVector>, BoundParsDim>& derivatives,
      const Covariance& startCov, const std::vector<double>& deviations) const {
    Jacobian jacobian;
    jacobian.setIdentity();
    jacobian.col(eLOC_0) = fitLinear(derivatives[eLOC_0], deviations);
    jacobian.col(eLOC_1) = fitLinear(derivatives[eLOC_1], deviations);
    jacobian.col(ePHI) = fitLinear(derivatives[ePHI], deviations);
    jacobian.col(eTHETA) = fitLinear(derivatives[eTHETA], deviations);
    jacobian.col(eQOP) = fitLinear(derivatives[eQOP], deviations);
    jacobian.col(eT) = fitLinear(derivatives[eT], deviations);
    return jacobian * startCov * jacobian.transpose();
  }

  /// @brief This function fits a linear function through the final state
  /// parametrisations
  ///
  /// @param [in] values Vector containing the final state parametrisations
  ///
  /// @return Vector containing the linear fit
  BoundVector fitLinear(const std::vector<BoundVector>& values,
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

  /// Propagator
  propagator_t m_propagator;
};
}  // namespace IntegrationTest

}  // namespace Acts