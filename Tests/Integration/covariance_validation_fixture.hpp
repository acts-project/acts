// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Propagator/Propagator.hpp"

namespace Acts {

namespace IntegrationTest {

using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;

template <typename propagator_t>
class RiddersPropagator
{

private:
  ///
  /// @note The result_type_helper struct and the action_list_t_result_t are here to allow a look'n'feel of this class like the Propagator itself
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
      

public:

	const unsigned int n = 4;
	std::array<double, n> deviations = {-2e-4, -1e-4, 1e-4, 2e-4}; // TODO: get rid of fixed dimension

	RiddersPropagator(propagator_t& propagator) : m_propagator(propagator){}
	
	template<typename stepper_t, typename navigator_t = detail::VoidNavigator>
	RiddersPropagator(stepper_t stepper, navigator_t navigator  = navigator_t()) : m_propagator(Propagator(stepper, navigator)) {}
	
template <typename parameters_t, typename action_list_t,
		typename aborter_list_t,
		template <typename, typename> class propagator_options_t>
Result<action_list_t_result_t<
  typename propagator_t::Stepper::template return_parameter_type<parameters_t>,
  action_list_t>>
propagate(
  const parameters_t& start,
  const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {
	// Launch nominal propagation and collect results
	auto nominalResult = m_propagator.propagate(start, options).value();
	const BoundVector& nominalParameters = nominalResult.endParameters->parameters();
	// Pick the surface of the propagation as target
	const Surface& surface = nominalResult.endParameters->referenceSurface();

	// Allow larger distances for the oscillation
	propagator_options_t<action_list_t, aborter_list_t> opts = options;
	opts.pathLimit *= 2.;
	
	// Derivations of each parameter around the nominal parameters
	std::array<std::array<BoundVector, n>, BoundParsDim> derivatives;
	
	// Wiggle each dimension individually
	for(unsigned int i = 0; i < BoundParsDim; i++)
	{
		derivatives[i] = wiggleDimension(opts, start, i, surface, nominalParameters);
	}
	
	// Exchange the result by Ridders Covariance
	if(nominalResult.endParameters->charge() == 0.)
		nominalResult.endParameters = std::make_unique<const CurvilinearParameters>(calculateCovariance(derivatives, *start.covariance()), nominalParameters.head(3), nominalParameters.template segment<3>(3), nominalParameters.tail(1));
	else	
		nominalResult.endParameters = std::make_unique<const CurvilinearParameters>(calculateCovariance(derivatives, *start.covariance()), nominalParameters.head(3), nominalParameters.template segment<3>(3), nominalResult.endParameters->charge(), nominalParameters.tail(1));
	return nominalResult;
  }
  
  template <typename parameters_t, typename surface_t, typename action_list_t,
            typename aborter_list_t,
            template <typename, typename> class propagator_options_t>
  Result<
      action_list_t_result_t<typename propagator_t::Stepper::template return_parameter_type<
                                 parameters_t, surface_t>,
                             action_list_t>>
  propagate(
      const parameters_t& start, const surface_t& target,
      const propagator_options_t<action_list_t, aborter_list_t>& options) const
  {
	// Launch nominal propagation and collect results
	auto nominalResult = m_propagator.propagate(start, target, options).value();
	const BoundVector& nominalParameters = nominalResult.endParameters->parameters();
	
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
	std::array<std::array<BoundVector, n>, BoundParsDim> derivatives;

	// Wiggle each dimension individually
	for(unsigned int i = 0; i < BoundParsDim; i++)
	{
		derivatives[i] = wiggleDimension(opts, start, i, target, nominalParameters);
	}
	
	// Exchange the result by Ridders Covariance
	if(nominalResult.endParameters->charge() == 0.)
		nominalResult.endParameters = std::make_unique<const BoundParameters>(opts.geoContext, calculateCovariance(derivatives, *start.covariance()), nominalParameters.head(3), nominalParameters.template segment<3>(3), nominalParameters.tail(1), target.getSharedPtr());
	else
		nominalResult.endParameters = std::make_unique<const BoundParameters>(opts.geoContext, calculateCovariance(derivatives, *start.covariance()), nominalParameters.head(3), nominalParameters.template segment<3>(3), nominalResult.endParameters->charge(), nominalParameters.tail(1), target.getSharedPtr());
	return nominalResult;
  }

private:  
     
    /// @brief This function wiggles one dimension of the starting parameters, performs the propagation to a surface and collects for each change of the start parameters the slope
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
	std::vector<BoundVector>
	wiggleDimension(const options_t& options, const parameters_t& startPars, const unsigned int param, const Surface& target, const BoundVector& nominal) const
	{
		// Storage of the results
		std::array<BoundVector, n> derivatives;
		derivatives.reserve(deviations.size());
		//~ for (double h : deviations) {
		for (unsigned int i = 0; i < n; i++) {
			double h = deviations[i];
		  parameters_t tp = startPars;
		  
		  // Treatment for theta
		  if(param == eTHETA)
		  {
			  const double current_theta = tp.template get<Acts::eTHETA>();
			  if (current_theta + h > M_PI) {
				h = M_PI - current_theta;
			  }
			  if (current_theta + h < 0) {
				h = -current_theta;
			  }
		  }
		  
		  // Modify start parameter and propagate
		  tp.template set<param>(options.geoContext,
										tp.template get<param>() + h);
		  const auto& r = m_propagator.propagate(tp, target, options).value();
		  
		  // Collect the slope
		  derivatives[i] = (r.endParameters->parameters() - nominal) / h;
		}
		return derivatives;
	}

	//~ /// @brief This function propagates the covariance matrix
	//~ ///
	//~ /// @tparam options_t PropagatorOptions object
	//~ ///
	//~ /// @param [in] options Options that store the variations
	//~ /// @param [in] derivatives Slopes of each modification of the parameters
	//~ /// @param [in] startCov Starting covariance
	//~ ///
	//~ /// @return Propagated covariance matrix
	std::unique_ptr<const Covariance>
	calculateCovariance(const std::array<std::array<BoundVector, n>, Acts::BoundParsDim>& derivatives, const Covariance& startCov) const
	{
		Jacobian jacobian;
		jacobian.setIdentity();
		jacobian.col(eLOC_0) = fitLinear(derivatives[eLOC_0], deviations);
		jacobian.col(eLOC_1) = fitLinear(derivatives[eLOC_1], deviations);
		jacobian.col(ePHI) = fitLinear(derivatives[ePHI], deviations);
		jacobian.col(eTHETA) = fitLinear(derivatives[eTHETA], deviations);
		jacobian.col(eQOP) = fitLinear(derivatives[eQOP], deviations);
		jacobian.col(eT) = fitLinear(derivatives[eT], deviations);
		return std::make_unique<const Covariance>(jacobian * startCov * jacobian.transpose());
    }

  template <unsigned long int N>
  BoundVector fitLinear(const std::vector<BoundVector>& values,
                               const std::array<double, N>& h) const {
    BoundVector A;
    BoundVector C;
    A.setZero();
    C.setZero();
    double B = 0;
    double D = 0;

    for (unsigned int i = 0; i < N; ++i) {
      A += h.at(i) * values.at(i);
      B += h.at(i);
      C += values.at(i);
      D += h.at(i) * h.at(i);
    }

    BoundVector b = (N * A - B * C) / (N * D - B * B);
    BoundVector a = (C - B * b) / N;

    return a;
  }

	/// Propagator
	propagator_t m_propagator;

};

template <typename T>
struct covariance_validation_fixture {
 public:
  covariance_validation_fixture(T propagator)
      : m_propagator(std::move(propagator)) {}

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
      std::vector<BoundVector>& derivatives) const {
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

  /// Numerical transport of covariance using the ridder's algorithm
  /// this is for covariance propagation validation
  /// it can either be used for curvilinear transport
  template <typename StartParameters, typename EndParameters, typename U>
  Covariance calculateCovariance(const StartParameters& startPars,
                                 const Covariance& startCov,
                                 const EndParameters& endPars,
                                 const U& options) const {
    // nominal propagation
    const auto& nominal = endPars.parameters();
    const Surface& dest = endPars.referenceSurface();

    // steps for estimating derivatives
    std::array<double, 4> h_steps = {{-4e-4, -2e-4, 2e-4, 4e-4}};
    if (dest.type() == Surface::Disc) {
      h_steps = {{-3e-5, -1e-5, 1e-5, 3e-5}};
    }

    // - for planar surfaces the dest surface is a perfect destination
    // surface for the numerical propagation, as reference frame
    // aligns with the referenceSurface.transform().rotation() at
    // at any given time
    //
    // - for straw & cylinder, where the error is given
    // in the reference frame that re-aligns with a slightly different
    // intersection solution

    // avoid stopping before the surface because of path length reached
    U var_options = options;
    var_options.pathLimit *= 2;

    // variation in x
    std::vector<BoundVector> x_derivatives;
    x_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eLOC_0>(options.geoContext,
                                    tp.template get<Acts::eLOC_0>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      x_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in y
    std::vector<BoundVector> y_derivatives;
    y_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eLOC_1>(options.geoContext,
                                    tp.template get<Acts::eLOC_1>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      y_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in phi
    std::vector<BoundVector> phi_derivatives;
    phi_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::ePHI>(options.geoContext,
                                  tp.template get<Acts::ePHI>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      phi_derivatives.push_back((r.endParameters->parameters() - nominal) / h);

      double phi0 = nominal(Acts::ePHI);
      double phi1 = r.endParameters->parameters()(Acts::ePHI);
      if (std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0))
        phi_derivatives.back()[Acts::ePHI] = (phi1 + 2. * M_PI - phi0) / h;
      else if (std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0))
        phi_derivatives.back()[Acts::ePHI] = (phi1 - 2. * M_PI - phi0) / h;
    }

    // variation in theta
    std::vector<BoundVector> theta_derivatives;
    theta_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      const double current_theta = tp.template get<Acts::eTHETA>();
      if (current_theta + h > M_PI) {
        h = M_PI - current_theta;
      }
      if (current_theta + h < 0) {
        h = -current_theta;
      }
      tp.template set<Acts::eTHETA>(options.geoContext,
                                    tp.template get<Acts::eTHETA>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      theta_derivatives.push_back((r.endParameters->parameters() - nominal) /
                                  h);
    }

    // variation in q/p
    std::vector<BoundVector> qop_derivatives;
    qop_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eQOP>(options.geoContext,
                                  tp.template get<Acts::eQOP>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      qop_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in t
    std::vector<BoundVector> t_derivatives;
    t_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eT>(options.geoContext,
                                tp.template get<Acts::eT>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      t_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    if (dest.type() == Surface::Disc &&
        (inconsistentDerivativesOnDisc(x_derivatives) ||
         inconsistentDerivativesOnDisc(y_derivatives) ||
         inconsistentDerivativesOnDisc(phi_derivatives) ||
         inconsistentDerivativesOnDisc(theta_derivatives) ||
         inconsistentDerivativesOnDisc(qop_derivatives) ||
         inconsistentDerivativesOnDisc(t_derivatives))) {
      return startCov;
    }

    Jacobian jacobian;
    jacobian.setIdentity();
    jacobian.col(Acts::eLOC_0) = fitLinear(x_derivatives, h_steps);
    jacobian.col(Acts::eLOC_1) = fitLinear(y_derivatives, h_steps);
    jacobian.col(Acts::ePHI) = fitLinear(phi_derivatives, h_steps);
    jacobian.col(Acts::eTHETA) = fitLinear(theta_derivatives, h_steps);
    jacobian.col(Acts::eQOP) = fitLinear(qop_derivatives, h_steps);
    jacobian.col(Acts::eT) = fitLinear(t_derivatives, h_steps);
    return jacobian * startCov * jacobian.transpose();
  }

 private:
  template <unsigned long int N>
  static BoundVector fitLinear(const std::vector<BoundVector>& values,
                               const std::array<double, N>& h) {
    BoundVector A;
    BoundVector C;
    A.setZero();
    C.setZero();
    double B = 0;
    double D = 0;

    for (unsigned int i = 0; i < N; ++i) {
      A += h.at(i) * values.at(i);
      B += h.at(i);
      C += values.at(i);
      D += h.at(i) * h.at(i);
    }

    BoundVector b = (N * A - B * C) / (N * D - B * B);
    BoundVector a = (C - B * b) / N;

    return a;
  }

  T m_propagator;
};

}  // namespace IntegrationTest

}  // namespace Acts