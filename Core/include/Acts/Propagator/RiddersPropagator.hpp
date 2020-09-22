// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

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
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;

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
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  template <typename parameters_t, typename propagator_options_t>
  Result<
      action_list_t_result_t<CurvilinearTrackParameters,
                             typename propagator_options_t::action_list_type>>
  propagate(const parameters_t& start,
            const propagator_options_t& options) const;

  /// @brief Propagation method targeting bound parameters
  ///
  /// @tparam parameters_t Type of the start parameters
  /// @tparam propagator_options_t Type of the propagator options
  ///
  /// @param [in] start Start parameters
  /// @param [in] options Options of the propagations
  ///
  /// @return Result of the propagation
  /// @note If the target surface is a disc, the resulting covariance may be
  /// inconsistent. In this case a zero matrix is returned.
  template <typename parameters_t, typename propagator_options_t>
  Result<action_list_t_result_t<
      BoundTrackParameters, typename propagator_options_t::action_list_type>>
  propagate(const parameters_t& start, const Surface& target,
            const propagator_options_t& options) const;

 private:
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
      const std::vector<BoundVector>& derivatives) const;

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
      const BoundVector& nominal, const std::vector<double>& deviations) const;

  /// @brief This function propagates the covariance matrix
  ///
  /// @param [in] derivatives Slopes of each modification of the parameters
  /// @param [in] startCov Starting covariance
  ///
  /// @return Propagated covariance matrix
  const Covariance calculateCovariance(
      const std::array<std::vector<BoundVector>, eBoundSize>& derivatives,
      const Covariance& startCov, const std::vector<double>& deviations) const;

  /// @brief This function fits a linear function through the final state
  /// parametrisations
  ///
  /// @param [in] values Vector containing the final state parametrisations
  ///
  /// @return Vector containing the linear fit
  BoundVector fitLinear(const std::vector<BoundVector>& values,
                        const std::vector<double>& deviations) const;

  /// Propagator
  propagator_t m_propagator;
};
}  // namespace Acts

#include "Acts/Propagator/RiddersPropagator.ipp"
