// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @brief Evaluater of the k_i's and elements of the transport matrix
/// D of the RKN4 stepping. This implementation involves energy loss due to
/// ioninisation, bremsstrahlung, pair production and photonuclear interaction
/// in the propagation and the jacobian. These effects will only occur if the
/// propagation is in a TrackingVolume with attached material.
struct DenseEnvironmentExtension
{

  /// Momentum at a certain point
  double currentMomentum = 0.;
  /// Particles momentum at k1
  double initialMomentum = 0.;
  /// Particles mass in SI units
  double massSI = 0.;
  /// Material that will be passed
  Material const* material = nullptr;
  /// Derivatives dLambda''dlambda at each sub-step point
  std::array<double, 4> dLdl;
  /// q/p at each sub-step in SI units
  std::array<double, 4> qop;
  /// Derivatives dPds at each sub-step
  std::array<double, 4> dPds;
  /// Derivative d(dEds)d(q/p) evaluated at the initial point
  double dgdqopValue = 0.;
  /// Derivative dEds at the initial point
  double g = 0.;

  /// Local store for conversion of momentum from SI to natural units
  const double conv = units::SI2Nat<units::MOMENTUM>(1);

  /// Energy loss calculator
  static const detail::IonisationLoss ionisationLoss;
  static const detail::RadiationLoss  radiationLoss;

  /// @brief Default constructor
  DenseEnvironmentExtension() = default;

  /// @brief Control function if the step evaluation would be valid
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @return Boolean flag if the step would be valid
  template <typename propagator_state_t, typename stepper_t>
  int
  bid(const propagator_state_t& state, const stepper_t& stepper) const
  {
    // Check for valid particle properties
    if (stepper.charge(state.stepping) == 0. || state.options.mass == 0.
        || stepper.momentum(state.stepping)
            < conv * state.options.momentumCutOff) {
      return 0;
    }

    // Check existence of a volume with material
    if (!state.navigation.currentVolume
        || !state.navigation.currentVolume->material()) {
      return 0;
    }
    return 2;
  }

  /// @brief Evaluater of the k_i's of the RKN4. For the case of i = 0 this
  /// step sets up member parameters, too.
  ///
  /// @tparam stepper_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] knew Next k_i that is evaluated
  /// @param [in] bField B-Field at the evaluation position
  /// @param [in] i Index of the k_i, i = [0, 3]
  /// @param [in] h Step size (= 0. ^ 0.5 * StepSize ^ StepSize)
  /// @param [in] kprev Evaluated k_{i - 1}
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool
  k(const propagator_state_t& state,
    const stepper_t&          stepper,
    Vector3D&                 knew,
    const Vector3D&           bField,
    const int                 i     = 0,
    const double              h     = 0.,
    const Vector3D&           kprev = Vector3D())
  {
    // i = 0 is used for setup and evaluation of k
    if (i == 0) {
      // Set up container for energy loss
      massSI   = units::Nat2SI<units::MASS>(state.options.mass);
      material = state.navigation.currentVolume->material().get();
      initialMomentum
          = units::Nat2SI<units::MOMENTUM>(stepper.momentum(state.stepping));
      currentMomentum = initialMomentum;
      qop[0]          = stepper.charge(state.stepping) / initialMomentum;
      initializeEnergyLoss(state);
      // Evaluate k
      knew = qop[0] * stepper.direction(state.stepping).cross(bField);
    } else {
      // Update parameters and check for momentum condition
      updateEnergyLoss(h, state.stepping, stepper, i);
      if (currentMomentum < state.options.momentumCutOff) {
        return false;
      }
      // Evaluate k
      knew = qop[i]
          * (stepper.direction(state.stepping) + h * kprev).cross(bField);
    }
    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in evaluation
  /// of the energy loss and the therewith constrained to keep the momentum
  /// after the step in reasonable values.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool
  finalize(propagator_state_t& state,
           const stepper_t&    stepper,
           const double        h) const
  {
    // Evaluate the new momentum
    double newMomentum = stepper.momentum(state.stepping)
        + conv * (h / 6.) * (dPds[0] + 2. * (dPds[1] + dPds[2]) + dPds[3]);

    // Break propagation if momentum becomes below cut-off
    if (units::Nat2SI<units::MOMENTUM>(newMomentum)
        < state.options.momentumCutOff) {
      return false;
    }

    // Add derivative dlambda/ds = Lambda''
    state.stepping.derivative(6)
        = -std::sqrt(state.options.mass * state.options.mass
                     + newMomentum * newMomentum)
        * units::SI2Nat<units::ENERGY>(g)
        / (newMomentum * newMomentum * newMomentum);

    // Update momentum
    state.stepping.p = newMomentum;

    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in the
  /// evaluation
  /// of the energy loss, the therewith constrained to keep the momentum
  /// after the step in reasonable values and the evaluation of the transport
  /// matrix.
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  /// @return Boolean flag if the calculation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool
  finalize(propagator_state_t& state,
           const stepper_t&    stepper,
           const double        h,
           ActsMatrixD<7, 7>& D) const
  {
    return finalize(state, stepper, h) && transportMatrix(state, stepper, h, D);
  }

private:
  /// @brief Evaluates the transport matrix D for the jacobian
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] h Step size
  /// @param [out] D Transport matrix
  /// @return Boolean flag if evaluation is valid
  template <typename propagator_state_t, typename stepper_t>
  bool
  transportMatrix(propagator_state_t& state,
                  const stepper_t&    stepper,
                  const double        h,
                  ActsMatrixD<7, 7>& D) const
  {
    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by dFdT and dGdT. The second is given by the top 3
    /// lines of the rightmost column. This is calculated by dFdL and dGdL.
    /// The remaining non-zero term is calculated directly. The naming of the
    /// variables is explained in eq. 11 and are directly related to the
    /// initial problem in eq. 7.
    /// The evaluation is based on propagating the parameters T and lambda
    /// (including g(lambda) and E(lambda)) as given in eq. 16 and evaluating
    /// the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dFdu' (due to the
    /// missing Lambda part) and only exists for dGdu' in dlambda/dlambda.

    auto& sd  = state.stepping.stepData;
    auto  dir = stepper.direction(state.stepping);

    D                   = ActsMatrixD<7, 7>::Identity();
    const double half_h = h * 0.5;

    // This sets the reference to the sub matrices
    // dFdx is already initialised as (3x3) zero
    auto dFdT = D.block<3, 3>(0, 3);
    auto dFdL = D.block<3, 1>(0, 6);
    // dGdx is already initialised as (3x3) identity
    auto dGdT = D.block<3, 3>(3, 3);
    auto dGdL = D.block<3, 1>(3, 6);

    ActsMatrixD<3, 3> dk1dT = ActsMatrixD<3, 3>::Zero();
    ActsMatrixD<3, 3> dk2dT = ActsMatrixD<3, 3>::Identity();
    ActsMatrixD<3, 3> dk3dT = ActsMatrixD<3, 3>::Identity();
    ActsMatrixD<3, 3> dk4dT = ActsMatrixD<3, 3>::Identity();

    ActsVectorD<3> dk1dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk2dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk3dL = ActsVectorD<3>::Zero();
    ActsVectorD<3> dk4dL = ActsVectorD<3>::Zero();

    /// Propagation of derivatives of dLambda''dlambda at each sub-step
    std::array<double, 4> jdL;

    // Evaluation of the rightmost column without the last term.
    jdL[0] = dLdl[0];
    dk1dL  = dir.cross(sd.B_first);
    jdL[1] = dLdl[1] * (1. + half_h * jdL[0]);
    dk2dL  = (1. + half_h * jdL[0]) * (dir + half_h * sd.k1).cross(sd.B_middle)
        + qop[1] * half_h * dk1dL.cross(sd.B_middle);
    jdL[2] = dLdl[2] * (1. + half_h * jdL[1]);
    dk3dL  = (1. + half_h * jdL[1]) * (dir + half_h * sd.k2).cross(sd.B_middle)
        + qop[2] * half_h * dk2dL.cross(sd.B_middle);
    jdL[3] = dLdl[3] * (1. + h * jdL[2]);
    dk4dL  = (1. + h * jdL[2]) * (dir + h * sd.k3).cross(sd.B_last)
        + qop[3] * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= qop[0];

    dk2dT += half_h * dk1dT;
    dk2dT = qop[1] * VectorHelpers::cross(dk2dT, sd.B_middle);

    dk3dT += half_h * dk2dT;
    dk3dT = qop[2] * VectorHelpers::cross(dk3dT, sd.B_middle);

    dk4dT += h * dk3dT;
    dk4dT = qop[3] * VectorHelpers::cross(dk4dT, sd.B_last);

    dFdT.setIdentity();
    dFdT += h / 6. * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = conv * h * h / 6. * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6. * (dk1dT + 2. * (dk2dT + dk3dT) + dk4dT);

    dGdL = conv * h / 6. * (dk1dL + 2. * (dk2dL + dk3dL) + dk4dL);

    // Evaluation of the dLambda''/dlambda term
    D(6, 6) += (h / 6.) * (jdL[0] + 2. * (jdL[1] + jdL[2]) + jdL[3]);

    return true;
  }

  /// @brief This function calculates the energy loss dE per path length ds of
  /// a particle through material. The energy loss consists of ionisation and
  /// radiation.
  /// @note The calculations use SI units and it is assumed that the arguments
  /// are given in SI units.
  ///
  /// @param [in] momentum Initial momentum of the particle
  /// @param [in] energy Initial energy of the particle
  /// @param [in] pdg PDG code of the particle
  /// @param [in] meanEnergyLoss Boolean indicator if mean or mode of the energy
  /// loss will be evaluated
  /// @return Infinitesimal energy loss
  double
  dEds(const double momentum,
       const double energy,
       const int    pdg,
       const bool   meanEnergyLoss) const
  {
    // Easy exit if material is invalid
    if (material->X0() == 0 || material->Z() == 0) {
      return 0.;
    }

    // Calculate energy loss by
    // a) ionisation
    double ionisationEnergyLoss = ionisationLoss
                                      .dEds(massSI,
                                            momentum * units::_c / energy,
                                            energy / (massSI * units::_c2),
                                            *(material),
                                            1.,
                                            meanEnergyLoss,
                                            true)
                                      .first;
    // b) radiation
    double radiationEnergyLoss
        = radiationLoss.dEds(energy, massSI, *(material), pdg, 1., true);

    // Rescaling for mode evaluation.
    // C.f. ATL-SOFT-PUB-2008-003 section 3. The mode evaluation for the energy
    // loss by ionisation can be directly evaluated.
    if (!meanEnergyLoss) {
      radiationEnergyLoss *= 0.15;
    }

    // Return sum of contributions
    return ionisationEnergyLoss + radiationEnergyLoss;
  }

  /// @brief This function calculates the derivation of g=dE/dx by d(q/p)
  ///
  /// @param [in] energy Initial energy of the particle
  /// @param [in] pdg PDG code of the particle
  /// @param [in] meanEnergyLoss Boolean indicator if mean or mode of the energy
  /// loss will be evaluated
  /// @return Derivative evaluated at the point defined by the
  /// function parameters
  double
  dgdqop(const double energy, const int pdg, const bool meanEnergyLoss) const
  {
    // Fast exit if material is invalid
    if (material->X0() == 0. || material->Z() == 0.
        || material->zOverAtimesRho() == 0.) {
      return 0.;
    }

    // Bethe-Bloch
    const double betheBlochDerivative
        = ionisationLoss.dqop(energy, qop[0], massSI, *(material), true, true);

    // Bethe-Heitler (+ pair production & photonuclear interaction for muons)
    const double radiationDerivative
        = radiationLoss.dqop(massSI, *(material), qop[0], energy, pdg, true);

    // Return the total derivative
    if (meanEnergyLoss) {
      return betheBlochDerivative + radiationDerivative;
    } else {
      // C.f. ATL-SOFT-PUB-2008-003 section 3
      return 0.9 * betheBlochDerivative + 0.15 * radiationDerivative;
    }
  }

  /// @brief Initializer of all parameters related to a RKN4 step with energy
  /// loss of a particle in material
  ///
  /// @tparam propagator_state_t Type of the state of the propagator
  /// @param [in] state Deliverer of configurations
  template <typename propagator_state_t>
  void
  initializeEnergyLoss(const propagator_state_t& state)
  {
    double E = std::sqrt(initialMomentum * initialMomentum * units::_c2
                         + massSI * massSI * units::_c4);
    // Use the same energy loss throughout the step.
    g = dEds(initialMomentum,
             E,
             state.options.absPdgCode,
             state.options.meanEnergyLoss);
    // Change of the momentum per path length
    // dPds = dPdE * dEds
    dPds[0] = g * E / (initialMomentum * units::_c2);
    if (state.stepping.covTransport) {
      // Calculate the change of the energy loss per path length and
      // inverse momentum
      if (state.options.includeGgradient) {
        dgdqopValue = dgdqop(
            E,
            state.options.absPdgCode,
            state.options
                .meanEnergyLoss);  // Use this value throughout the step.
      }
      // Calculate term for later error propagation
      dLdl[0] = (-qop[0] * qop[0] * g * E
                     * (3.
                        - (initialMomentum * initialMomentum * units::_c2)
                            / (E * E))
                 - qop[0] * qop[0] * qop[0] * E * dgdqopValue)
          / units::_c2;
    }
  }

  /// @brief Update of the kinematic parameters of the RKN4 sub-steps after
  /// initialization with energy loss of a particle in material
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @tparam stepper_t Type of the stepper
  /// @param [in] h Stepped distance of the sub-step (1-3)
  /// @param [in] state State of the stepper
  /// @param [in] i Index of the sub-step (1-3)
  template <typename stepper_state_t, typename stepper_t>
  void
  updateEnergyLoss(const double           h,
                   const stepper_state_t& state,
                   const stepper_t&       stepper,
                   const int              i)
  {
    // Update parameters related to a changed momentum
    currentMomentum = initialMomentum + h * dPds[i - 1];

    double E = std::sqrt(currentMomentum * currentMomentum * units::_c2
                         + massSI * massSI * units::_c4);
    dPds[i] = g * E / (currentMomentum * units::_c2);
    qop[i]  = stepper.charge(state) / currentMomentum;
    // Calculate term for later error propagation
    if (state.covTransport) {
      dLdl[i] = (-qop[i] * qop[i] * g * E
                     * (3.
                        - (currentMomentum * currentMomentum * units::_c2)
                            / (E * E))
                 - qop[i] * qop[i] * qop[i] * E * dgdqopValue)
          / units::_c2;
    }
  }
};

template <typename action_list_t  = ActionList<>,
          typename aborter_list_t = AbortList<>>
struct DenseStepperPropagatorOptions
    : public PropagatorOptions<action_list_t, aborter_list_t>
{
  /// Toggle between mean and mode evaluation of energy loss
  bool meanEnergyLoss = true;

  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool includeGgradient = true;

  /// Cut-off value for the momentum in SI units
  double momentumCutOff = 0.;

  /// @brief Expand the Options with extended aborters
  ///
  /// @tparam extended_aborter_list_t Type of the new aborter list
  ///
  /// @param aborters The new aborter list to be used (internally)
  template <typename extended_aborter_list_t>
  DenseStepperPropagatorOptions<action_list_t, extended_aborter_list_t>
  extend(extended_aborter_list_t aborters) const
  {
    DenseStepperPropagatorOptions<action_list_t, extended_aborter_list_t>
        eoptions;
    // Copy the options over
    eoptions.direction       = this->direction;
    eoptions.absPdgCode      = this->absPdgCode;
    eoptions.mass            = this->mass;
    eoptions.maxSteps        = this->maxSteps;
    eoptions.maxStepSize     = this->maxStepSize;
    eoptions.targetTolerance = this->targetTolerance;
    eoptions.pathLimit       = this->pathLimit;
    eoptions.loopProtection  = this->loopProtection;
    eoptions.loopFraction    = this->loopFraction;
    // Output option
    eoptions.debug         = this->debug;
    eoptions.debugString   = this->debugString;
    eoptions.debugPfxWidth = this->debugPfxWidth;
    eoptions.debugMsgWidth = this->debugMsgWidth;
    // Action / abort list
    eoptions.actionList = this->actionList;
    eoptions.abortList  = std::move(aborters);
    // Copy dense environment specific parameters
    eoptions.meanEnergyLoss   = meanEnergyLoss;
    eoptions.includeGgradient = includeGgradient;
    eoptions.momentumCutOff   = momentumCutOff;
    // And return the options
    return eoptions;
  }
};

}  // namespace Acts
