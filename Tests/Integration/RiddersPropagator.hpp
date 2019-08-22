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

/// @brief This class performs the Ridders algorithm to estimate the propagation of the covariance to a certain point in space.
///
/// The algorithm is based on the small deviations of the start parameters based on their uncertainty at the beginning of the propgation. This deviation is represented here by a vector of relative deviations of these parameters and fix for all parameters. So, a common choice has to be found that is able to actually fit into the order of magnitude of the uncertainty of each parameter. Using these deviations, the propagation is repeated multiple times and the final covariance matrix at a given target surface is afterwards evaluated by first order derivatives of the final state parameters wrt. the inital parameters. Therefore this evaluation represents a first order approximation of the transport jacobian. Since performing multiple propagations and a numerical evaluation of the covariance requires more time than a single propagation towards a target + a common propagation of the covariance, this class just serves to verify the results of the latter classes.
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

	std::vector<double> deviations = {-2e-4, -1e-4, 1e-4, 2e-4};

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
	std::array<std::vector<BoundVector>, BoundParsDim> derivatives;

	// Wiggle each dimension individually
	for(unsigned int i = 0; i < BoundParsDim; i++)
	{
		derivatives[i] = wiggleDimension(opts, start, i, surface, nominalParameters);
	}
	// Exchange the result by Ridders Covariance
	const FullParameterSet& parSet = nominalResult.endParameters->getParameterSet();
	FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
	mParSet->setCovariance((start.covariance() != nullptr) ? calculateCovariance(derivatives, *start.covariance()) : nullptr);

	return std::move(nominalResult);
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
	std::array<std::vector<BoundVector>, BoundParsDim> derivatives;

	// Wiggle each dimension individually
	for(unsigned int i = 0; i < BoundParsDim; i++)
	{
		derivatives[i] = wiggleDimension(opts, start, i, target, nominalParameters);
	}
	// Exchange the result by Ridders Covariance
	const FullParameterSet& parSet = nominalResult.endParameters->getParameterSet();
	FullParameterSet* mParSet = const_cast<FullParameterSet*>(&parSet);
	mParSet->setCovariance((start.covariance() != nullptr) ? calculateCovariance(derivatives, *start.covariance()) : nullptr);

	return std::move(nominalResult);
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
		std::vector<BoundVector> derivatives;
		derivatives.reserve(deviations.size());
		for (double h : deviations) {
		  parameters_t tp = startPars;
		  
		  // Treatment for theta
		  if(param == eTHETA)
		  {
			  const double current_theta = tp.template get<eTHETA>();
			  if (current_theta + h > M_PI) {
				h = M_PI - current_theta;
			  }
			  if (current_theta + h < 0) {
				h = -current_theta;
			  }
		  }
		  
		  // Modify start parameter and propagate
		  switch(param)
		  {
			  case 0:
			  {
				tp.template set<eLOC_0>(options.geoContext,
											tp.template get<eLOC_0>() + h);
				break;
			  }
			  case 1:
			  {
				tp.template set<eLOC_1>(options.geoContext,
											tp.template get<eLOC_1>() + h);
				break;
			  }
			  case 2:
			  {
				tp.template set<ePHI>(options.geoContext,
											tp.template get<ePHI>() + h);
				break;
			  }
			  case 3:
			  {
				tp.template set<eTHETA>(options.geoContext,
											tp.template get<eTHETA>() + h);
				break;
			  }
			  case 4:
			  {
				tp.template set<eQOP>(options.geoContext,
											tp.template get<eQOP>() + h);
				break;
			  }
			  case 5:
			  {
				tp.template set<eT>(options.geoContext,
											tp.template get<eT>() + h);
				break;
			  }
			  default:
				return {};
		  }
		  const auto& r = m_propagator.propagate(tp, target, options).value();
		  // Collect the slope
		  derivatives.push_back((r.endParameters->parameters() - nominal) / h);
		}
		return derivatives;
	}

	/// @brief This function propagates the covariance matrix
	///
	/// @param [in] derivatives Slopes of each modification of the parameters
	/// @param [in] startCov Starting covariance
	///
	/// @return Propagated covariance matrix
	std::unique_ptr<const Covariance>
	calculateCovariance(const std::array<std::vector<BoundVector>, BoundParsDim>& derivatives, const Covariance& startCov) const
	{
		Jacobian jacobian;
		jacobian.setIdentity();
		jacobian.col(eLOC_0) = fitLinear(derivatives[eLOC_0]);
		jacobian.col(eLOC_1) = fitLinear(derivatives[eLOC_1]);
		jacobian.col(ePHI) = fitLinear(derivatives[ePHI]);
		jacobian.col(eTHETA) = fitLinear(derivatives[eTHETA]);
		jacobian.col(eQOP) = fitLinear(derivatives[eQOP]);
		jacobian.col(eT) = fitLinear(derivatives[eT]);
		return std::make_unique<const Covariance>(jacobian * startCov * jacobian.transpose());
    }
    
  /// @brief This function fits a linear function through the final state parametrisations
  ///
  /// @param [in] values Vector containing the final state parametrisations
  ///
  /// @return Vector containing the linear fit
  BoundVector fitLinear(const std::vector<BoundVector>& values) const {
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