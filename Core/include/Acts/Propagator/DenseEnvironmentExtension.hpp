// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @brief Evaluater of the k_i's and elements of the transport matrix
/// D of the RKN4 stepping. This implementation involves energy loss due to
/// ioninisation, bremsstrahlung, pair production and photonuclear interaction
/// in the propagation and the jacobian. These effects will only occur if the
/// propagation is in a TrackingVolume with attached material.
struct DenseEnvironmentExtension
{
  /// @brief This struct serves as data container to keep track of all
  /// parameters that are related to an energy loss of a particle in matter.
  struct EnergyLossData
  {
    /// Momentum at a certain point
    double currentMomentum = 0.;
    /// Particles momentum at k1
    double initialMomentum = 0.;
    /// Particles mass in SI units
    double massSI = 0.;
    /// Material that will be passed
    std::shared_ptr<const Material> material = nullptr;
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
  };

  /// Toggle between mean and mode evaluation of energy loss
  bool meanEnergyLoss = true;

  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool includeGgradient = true;

  /// Cut-off value for the momentum in SI units
  double momentumCutOff = 0.;

  /// Data container for the energy loss
  EnergyLossData eld;

  /// Local store for conversion of momentum from SI to natural units
  const double conv = units::SI2Nat<units::MOMENTUM>(1);

  /// @brief Default constructor
  DenseEnvironmentExtension() = default;

  /// @brief Control function if the step evaluation would be valid
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  /// @return Boolean flag if the step would be valid
  template <typename stepper_state_t>
  bool
  validExtensionForStep(const stepper_state_t& state) const
  {
    // Check for valid particle properties
    if (state.q == 0. || state.p < conv * momentumCutOff) return false;

    // Check existence of a volume with material
    if (!state.volume || !(*state.volume) || !(*state.volume)->material())
      return false;
    return true;
  }

  /// @brief Evaluater of the k_i's of the RKN4. For the case of i = 0 this
  /// step sets up member parameters, too.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  /// @param [out] knew Next k_i that is evaluated
  /// @param [in] bField B-Field at the evaluation position
  /// @param [in] i Index of the k_i
  /// @param [in] h Step size (= 0. ^ 0.5 * StepSize ^ StepSize)
  /// @param [in] kprev Evaluated k_{i - 1}
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t>
  bool
  k(const stepper_state_t& state,
    Vector3D&              knew,
    const Vector3D&        bField,
    const int              i     = 0,
    const double           h     = 0.,
    const Vector3D&        kprev = Vector3D())
  {
    // i = 0 is used for setup and evaluation of k
    if (i == 0) {
      // Set up container for energy loss
      eld.massSI          = units::Nat2SI<units::MASS>(state.mass);
      eld.material        = (*state.volume)->material();
      eld.initialMomentum = units::Nat2SI<units::MOMENTUM>(state.p);
      eld.currentMomentum = eld.initialMomentum;
      eld.qop[0]          = state.q / eld.initialMomentum;
      initializeEnergyLoss(state);
      // Evaluate k
      knew = eld.qop[0] * state.dir.cross(bField);
    } else {
      // Update parameters and check for momentum condition
      updateEnergyLoss(h, state, i);
      if (eld.currentMomentum < momentumCutOff) return false;
      // Evaluate k
      knew = eld.qop[i] * (state.dir + h * kprev).cross(bField);
    }
    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in evaluation
  /// of the energy loss and the therewith constrained to keep the momentum
  /// after the step in reasonable values.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state State of the stepper
  /// @param [in] h Step size
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t>
  bool
  finalize(stepper_state_t& state, const double h) const
  {
    // Evaluate the new momentum
    double newMomentum = state.p
        + conv * (h / 6.)
            * (eld.dPds[0] + 2. * (eld.dPds[1] + eld.dPds[2]) + eld.dPds[3]);

    // Break propagation if momentum becomes below cut-off
    if (units::Nat2SI<units::MOMENTUM>(newMomentum) < momentumCutOff)
      return false;

    // Update momentum
    state.p = newMomentum;
    return true;
  }

  /// @brief After a RKN4 step was accepted by the stepper this method has an
  /// additional veto on the quality of the step. The veto lies in the
  /// evaluation
  /// of the energy loss, the therewith constrained to keep the momentum
  /// after the step in reasonable values and the evaluation of the transport
  /// matrix.
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @tparam stepper_data_t Type of the data collected in the step
  /// @param [in] state State of the stepper
  /// @param [in] h Step size
  /// @param [in] data Data of B-field and k_i's
  /// @param [out] D Transport matrix
  /// @return Boolean flag if the calculation is valid
  template <typename stepper_state_t, typename stepper_data_t>
  bool
  finalize(stepper_state_t&      state,
           const double          h,
           const stepper_data_t& data,
           ActsMatrixD<7, 7>& D) const
  {
    return finalize(state, h) && transportMatrix(state.dir, h, data, D);
  }

protected:
  /// @brief Evaluates the transport matrix D for the jacobian
  ///
  /// @param [in] dir Direction of the particle
  /// @param [in] h Step size
  /// @param [in] sd Data of B-field and k_i's
  /// @param [out] D Transport matrix
  /// @return Boolean flag if evaluation is valid
  template <typename stepper_data_t>
  bool
  transportMatrix(const Vector3D&       dir,
                  const double          h,
                  const stepper_data_t& sd,
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
    /// u_n nor in u_n'. The second and forth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dGdu' (due to the
    /// missing Lambda part) and only exists for dFdu' in dlambda/dlambda.

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
    jdL[0] = eld.dLdl[0];
    dk1dL  = dir.cross(sd.B_first);
    jdL[1] = eld.dLdl[1] * (1. + half_h * jdL[0]);
    dk2dL  = (1. + half_h * jdL[0]) * (dir + half_h * sd.k1).cross(sd.B_middle)
        + eld.qop[1] * half_h * dk1dL.cross(sd.B_middle);
    jdL[2] = eld.dLdl[2] * (1. + half_h * jdL[1]);
    dk3dL  = (1. + half_h * jdL[1]) * (dir + half_h * sd.k2).cross(sd.B_middle)
        + eld.qop[2] * half_h * dk2dL.cross(sd.B_middle);
    jdL[3] = eld.dLdl[3] * (1. + h * jdL[2]);
    dk4dL  = (1. + h * jdL[2]) * (dir + h * sd.k3).cross(sd.B_last)
        + eld.qop[3] * h * dk3dL.cross(sd.B_last);

    dk1dT(0, 1) = sd.B_first.z();
    dk1dT(0, 2) = -sd.B_first.y();
    dk1dT(1, 0) = -sd.B_first.z();
    dk1dT(1, 2) = sd.B_first.x();
    dk1dT(2, 0) = sd.B_first.y();
    dk1dT(2, 1) = -sd.B_first.x();
    dk1dT *= eld.qop[0];

    dk2dT += half_h * dk1dT;
    dk2dT *= VectorHelpers::cross(dk2dT, sd.B_middle);
    dk2dT *= eld.qop[1];

    dk3dT += half_h * dk2dT;
    dk3dT *= VectorHelpers::cross(dk3dT, sd.B_middle);
    dk3dT *= eld.qop[2];

    dk4dT += h * dk3dT;
    dk4dT *= VectorHelpers::cross(dk4dT, sd.B_last);
    dk4dT *= eld.qop[3];

    dFdT.setIdentity();
    dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
    dFdT *= h;

    dFdL = conv * h * h / 6 * (dk1dL + dk2dL + dk3dL);

    dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

    dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

    // Evaluation of the dLambda''/dlambda term
    D(6, 6) += conv * (h / 6.) * (jdL[0] + 2. * (jdL[1] + jdL[2]) + jdL[3]);
    return true;
  }

  /// Energy loss calculator
  detail::IonisationLoss ionisationLoss;
  detail::RadiationLoss  radiationLoss;

  /// @brief This function calculates the energy loss dE per path length ds of
  /// a particle through material. The energy loss consists of ionisation and
  /// radiation.
  /// @note The calculations use SI units and it is assumed that the arguments
  /// are given in SI units.
  ///
  /// @tparam material_t Type of the material
  /// @param [in] momentum Initial momentum of the particle
  /// @param [in] energy Initial energy of the particle
  /// @param [in] material Penetrated material
  /// @param [in] pdg PDG code of the particle
  /// @return Infinitesimal energy loss
  template <typename material_t>
  double
  dEds(const double      momentum,
       const double      energy,
       const material_t& material,
       const int         pdg) const
  {
    // Easy exit if material is invalid
    if (material.X0() == 0 || material.Z() == 0) return 0.;

    // Calculate energy loss by
    // a) ionisation
    double ionisationEnergyLoss
        = ionisationLoss(eld.massSI,
                         momentum * units::_c / energy,
                         energy / (eld.massSI * units::_c2),
                         material,
                         1.,
                         meanEnergyLoss,
                         true)
              .first;
    // b) radiation
    double radiationEnergyLoss
        = radiationLoss(energy, eld.massSI, material, pdg, 1., true);

    // Rescaling for mode evaluation.
    // C.f. ATL-SOFT-PUB-2008-003 section 3. The mode evaluation for the energy
    // loss by ionisation can be directly evaluated.
    if (!meanEnergyLoss) radiationEnergyLoss *= 0.15;

    // Return sum of contributions
    return ionisationEnergyLoss + radiationEnergyLoss;
  }

  /// @brief This function calculates the derivation of g=dE/dx by d(q/p)
  ///
  /// @tparam material_t Type of the material
  /// @param [in] energy Initial energy of the particle
  /// @param [in] qop Initial value of q/p of the particle
  /// @param [in] material Penetrated material
  /// @param [in] pdg PDG code of the particle
  /// @return Derivative evaluated at the point defined by the
  /// function parameters
  template <typename material_t>
  double
  dgdqop(const double      energy,
         const double      qop,
         const material_t& material,
         const int         pdg) const
  {
    // Fast exit if material is invalid
    if (material.X0() == 0. || material.Z() == 0.
        || material.zOverAtimesRho() == 0.)
      return 0.;

    // Bethe-Bloch
    const double betheBlochDerivative
        = ionisationLoss.dqop(energy, qop, eld.massSI, material, true, true);

    // Bethe-Heitler (+ pair production & photonuclear interaction for muons)
    const double radiationDerivative
        = radiationLoss.dqop(eld.massSI, material, qop, energy, pdg, true);

    // Return the total derivative
    if (meanEnergyLoss)
      return betheBlochDerivative + radiationDerivative;
    else
      // C.f. ATL-SOFT-PUB-2008-003 section 3
      return 0.9 * betheBlochDerivative + 0.15 * radiationDerivative;
  }

  /// @brief Initializer of all parameters related to a RKN4 step with energy
  /// loss of a particle in material
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] state Deliverer of configurations
  template <typename stepper_state_t>
  void
  initializeEnergyLoss(const stepper_state_t& state)
  {
    double E = std::sqrt(eld.initialMomentum * eld.initialMomentum * units::_c2
                         + eld.massSI * eld.massSI * units::_c4);
    // Use the same energy loss throughout the step.
    eld.g = dEds(eld.initialMomentum, E, *(eld.material), state.pdg);
    // Change of the momentum per path length
    // dPds = dPdE * dEds
    eld.dPds[0] = eld.g * E / (eld.initialMomentum * units::_c2);
    if (state.covTransport) {
      // Calculate the change of the energy loss per path length and
      // inverse momentum
      if (includeGgradient) {
        eld.dgdqopValue
            = dgdqop(E,
                     eld.qop[0],
                     *(eld.material),
                     state.pdg);  // Use this value throughout the step.
      }
      // Calculate term for later error propagation
      eld.dLdl[0]
          = (-eld.qop[0] * eld.qop[0] * eld.g * E
                 * (3.
                    - (eld.initialMomentum * eld.initialMomentum * units::_c2)
                        / (E * E))
             - eld.qop[0] * eld.qop[0] * eld.qop[0] * E * eld.dgdqopValue)
          / units::_c3;
    }
  }

  /// @brief Update of the kinematic parameters of the RKN4 sub-steps after
  /// initialization with energy loss of a particle in material
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @param [in] h Stepped distance of the sub-step (1-3)
  /// @param [in] state State of the stepper
  /// @param [in] i Index of the sub-step (1-3)
  template <typename stepper_state_t>
  void
  updateEnergyLoss(const double h, const stepper_state_t& state, const int i)
  {
    // Update parameters related to a changed momentum
    eld.currentMomentum = eld.initialMomentum + h * eld.dPds[i - 1];

    double E = std::sqrt(eld.currentMomentum * eld.currentMomentum * units::_c2
                         + eld.massSI * eld.massSI * units::_c4);
    eld.dPds[i] = eld.g * E / (eld.currentMomentum * units::_c2);
    eld.qop[i]  = state.q / eld.currentMomentum;
    // Calculate term for later error propagation
    if (state.covTransport) {
      eld.dLdl[i]
          = (-eld.qop[i] * eld.qop[i] * eld.g * E
                 * (3.
                    - (eld.currentMomentum * eld.currentMomentum * units::_c2)
                        / (E * E))
             - eld.qop[i] * eld.qop[i] * eld.qop[i] * E * eld.dgdqopValue)
          / units::_c3;
    }
  }
};

/// @brief Actor as configurator of the Stepper for working with the
/// DenseEnvironmentExtension. It sets up steering properties by the user.
struct DenseEnvironmentExtensionActor
{
  // Configurations for Stepper
  /// Toggle between mean and mode evaluation of energy loss
  bool m_meanEnergyLoss = true;
  /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
  bool m_includeGgradient = true;
  /// Cut-off value for the momentum in SI units
  double m_momentumCutOff = 0.;

  /// @brief Main call operator for setting up stepper properties
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param [in, out] state State of the propagator
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {
    // Initialize all parameters
    if (state.stepping.pathAccumulated == 0.) {
      // Initialize user defined parameters
      DenseEnvironmentExtension& ext
          = state.stepping.extension
                .template get<Acts::DenseEnvironmentExtension>();
      ext.momentumCutOff   = m_momentumCutOff;
      ext.meanEnergyLoss   = m_meanEnergyLoss;
      ext.includeGgradient = m_includeGgradient;
    }
  }
};

}  // namespace Acts
