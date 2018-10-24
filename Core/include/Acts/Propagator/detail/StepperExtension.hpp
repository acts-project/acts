// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"

namespace Acts {
namespace detail {

  /// @brief Default evaluater of the k_i's and elements of the transport matrix
  /// D of the RKN4 stepping. This is a pure implementation by textbook.
  // TODO: note used unit systems
  struct DefaultExtension
  {
    /// @brief Default constructor
    DefaultExtension() = default;

    /// Local store for q/p
    double qop = 0.;

    /// Local store for conversion of momentum from SI to natural units
    const double conv = units::SI2Nat<units::MOMENTUM>(1);

    /// @brief Evaluates the k1 of the RKN4 by textbook
    ///
    /// @tparam stepper_state_t Type of the state of the stepper
    /// @param [in] state State of the stepper
    /// @param [in] bField B-Field at the first position
    /// @param [out] k1 Vector of k1
    template <typename stepper_state_t>
    void
    evaluatek1(const stepper_state_t& state,
               const Vector3D&        bField,
               Vector3D&              k1)
    {
      // Store qop, it is always used if valid
      qop = state.q / units::Nat2SI<units::MOMENTUM>(state.p);

      k1 = qop * state.dir.cross(bField);
    }

    /// @brief Evaluates the k2 of the RKN4 by textbook
    ///
    /// @tparam stepper_state_t Type of the state of the stepper
    /// @param [in] dir Direction of the particle
    /// @param [in] half_h Half the step size
    /// @param [in] k1 Evaluated value of k1
    /// @param [in] bField B-Field at the second position
    /// @param [out] k2 Vector of k2
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek2(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k1,
               const Vector3D&        bField,
               Vector3D&              k2) const
    {
      k2 = qop * (state.dir + half_h * k1).cross(bField);
      return true;
    }

    /// @brief Evaluates the k3 of the RKN4 by textbook
    ///
    /// @tparam stepper_state_t Type of the state of the stepper
    /// @param [in] dir Direction of the particle
    /// @param [in] half_h Half the step size
    /// @param [in] k2 Evaluated value of k2
    /// @param [in] bField B-Field at the second position
    /// @param [out] k3 Vector of k3
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek3(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k2,
               const Vector3D&        bField,
               Vector3D&              k3) const
    {
      k3 = qop * (state.dir + half_h * k2).cross(bField);
      return true;
    }

    /// @brief Evaluates the k4 of the RKN4 by textbook
    ///
    /// @tparam stepper_state_t Type of the state of the steppers
    /// @param [in] dir Direction of the particle
    /// @param [in] h Step size
    /// @param [in] k3 Evaluated value of k3
    /// @param [in] bField B-Field at the last position
    /// @param [out] k4 Vector of k4
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek4(const stepper_state_t& state,
               const double           h,
               const Vector3D&        k3,
               const Vector3D&        bField,
               Vector3D&              k4) const
    {
      k4 = qop * (state.dir + h * k3).cross(bField);
      return true;
    }

	/// @brief Veto function after a RKN4 step was accepted by judging on the error of the step. Since the textbook does not deliver further vetos, this is a dummy function.
	///
	template<typename stepper_state_t>
    bool
    finalizeStep(stepper_state_t&, const double)
    {
		return true;
	}
	
    /// @brief Evaluates the transport matrix D for the jacobian
    ///
    /// @param [in] dir Direction of the particle
    /// @param [in] bFiedl1 B-field at the first position
    /// @param [in] bField2 B-field at the middle position
    /// @param [in] bField3 B-field at the last position
    /// @param [in] h Step size
    /// @param [in] k1 Vector of k1
    /// @param [in] k2 Vector of k2
    /// @param [in] k3 Vector of k3
    /// @param [out] D Transport matrix
    /// @return Boolean flag is step evaluation is valid
    bool
    evaluateD(const Vector3D& dir,
              const Vector3D& bField1,
              const Vector3D& bField2,
              const Vector3D& bField3,
              const double    h,
              const Vector3D& k1,
              const Vector3D& k2,
              const Vector3D& k3,
              ActsMatrixD<7, 7>& D) const
    {
      double half_h = h * 0.5;
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

      // For the case without energy loss
      dk1dL = dir.cross(bField1);
      dk2dL = (dir + half_h * k1).cross(bField2)
          + qop * half_h * dk1dL.cross(bField2);
      dk3dL = (dir + half_h * k2).cross(bField2)
          + qop * half_h * dk2dL.cross(bField2);
      dk4dL = (dir + h * k3).cross(bField3) + qop * h * dk3dL.cross(bField3);

      dk1dT(0, 1) = bField1.z();
      dk1dT(0, 2) = -bField1.y();
      dk1dT(1, 0) = -bField1.z();
      dk1dT(1, 2) = bField1.x();
      dk1dT(2, 0) = bField1.y();
      dk1dT(2, 1) = -bField1.x();
      dk1dT *= qop;

      dk2dT += half_h * dk1dT;
      dk2dT *= cross(dk2dT, bField2);
      dk2dT *= qop;

      dk3dT += half_h * dk2dT;
      dk3dT *= cross(dk3dT, bField2);
      dk3dT *= qop;

      dk4dT += h * dk3dT;
      dk4dT *= cross(dk4dT, bField3);
      dk4dT *= qop;

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h * h / 6 * (dk1dL + dk2dL + dk3dL);

      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      return true;
    }
  };

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
      std::shared_ptr<const Material> material;
      /// Derivatives dLambda''dlambda at each sub-step point
      std::array<double, 4> dLdl;
      /// q/p at each sub-step
      std::array<double, 4> qop;
      /// Derivatives dPds at each sub-step
      std::array<double, 4> dPds;
      /// Propagation of derivatives of dLambda''dlambda at each sub-step
      std::array<double, 4> jdL;
      /// Derivative d(dEds)d(q/p) evaluated at the initial point
      double dgdqopValue = 0.;
      /// Derivative dEds at the initial point
      double g = 0.;
    };

    /// Mass
    double mass = 0.;

    /// PDG code
    int pdg = 0;

    /// Volume with material that is passed
    TrackingVolume const* const* volume = nullptr;

    /// Toggle between mean and mode evaluation of energy loss
    bool meanEnergyLoss = true;

    /// Tolerance for the error of the integration
    double tolerance = 5e-5;

    /// Boolean flag for inclusion of d(dEds)d(q/p) into energy loss
    bool includeGgradient = true;

    /// Cut-off value for the momentum in SI units
    double momentumCutOff = 0.;

    /// Cut-off value for the step size
    double stepSizeCutOff = 0.;

    /// Data container for the energy loss
    EnergyLossData elData;

    /// Local store for conversion of momentum from SI to natural units
    const double conv = units::SI2Nat<units::MOMENTUM>(1);

    /// @brief Default constructor
    DenseEnvironmentExtension() = default;

    /// @brief Evaluates the k1 of the RKN4 in a dense environment
    ///
    /// @tparam steppter_state_t Type of the state of the stepper
    /// @param [in] state State of the stepper
    /// @param [in] bField B-Field at the first position
    /// @param [out] k1 Vector of k1
    template <typename stepper_state_t>
    void
    evaluatek1(const stepper_state_t& state,
               const Vector3D&        bField,
               Vector3D&              k1)
    {
      if (!volume || !(*volume) || !(*volume)->material()) return;

      // Set up container for energy loss
      elData.massSI          = units::Nat2SI<units::MASS>(mass);
      elData.material        = (*volume)->material();
      elData.initialMomentum = units::Nat2SI<units::MOMENTUM>(state.p);
      elData.currentMomentum = elData.initialMomentum;
      elData.qop[0]          = state.q / elData.initialMomentum;
      k1                     = elData.qop[0] * state.dir.cross(bField);
      initializeEnergyLoss(state);
    }

    /// @brief Evaluates the k2 of the RKN4 in a dense environment
    ///
    /// @tparam steppter_state_t Type of the state of the stepper
    /// @param [in] state State of the stepper
    /// @param [in] half_h Half the step size
    /// @param [in] k1 Evaluated value of k1
    /// @param [in] bField B-Field at the second position
    /// @param [out] k2 Vector of k2
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek2(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k1,
               const Vector3D&        bField,
               Vector3D&              k2)
    {
      if (!volume || !(*volume) || !(*volume)->material()) return true;

      // Update parameters and check for momentum condition
      updateEnergyLoss(half_h, state, 1);
      if (elData.currentMomentum < momentumCutOff) return false;
      k2 = elData.qop[1] * (state.dir + half_h * k1).cross(bField);
      return true;
    }

    /// @brief Evaluates the k3 of the RKN4 in a dense environment
    ///
    /// @tparam steppter_state_t Type of the state of the stepper
    /// @param [in] state State of the stepper
    /// @param [in] half_h Half the step size
    /// @param [in] k2 Evaluated value of k2
    /// @param [in] bField B-Field at the second position
    /// @param [out] k3 Vector of k3
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek3(const stepper_state_t& state,
               const double           half_h,
               const Vector3D&        k2,
               const Vector3D&        bField,
               Vector3D&              k3)
    {
      if (!volume || !(*volume) || !(*volume)->material()) return true;

      // Update parameters and check for momentum condition
      updateEnergyLoss(half_h, state, 2);
      if (elData.currentMomentum < momentumCutOff) return false;
      k3 = elData.qop[2] * (state.dir + half_h * k2).cross(bField);
      return true;
    }

    /// @brief Evaluates the k4 of the RKN4 in a dense environment
    ///
    /// @tparam steppter_state_t Type of the state of the stepper
    /// @param [in] state State of the stepper
    /// @param [in] h Step size
    /// @param [in] k3 Evaluated value of k3
    /// @param [in] bField B-Field at the last position
    /// @param [out] k4 Vector of k4
    /// @return Boolean flag if step evaluation is valid
    template <typename stepper_state_t>
    bool
    evaluatek4(const stepper_state_t& state,
               const double           h,
               const Vector3D&        k3,
               const Vector3D&        bField,
               Vector3D&              k4)
    {
      if (!volume || !(*volume) || !(*volume)->material()) return true;

      // Update parameters and check for momentum condition
      updateEnergyLoss(h, state, 3);
      if (elData.currentMomentum < momentumCutOff)
        return false;
      k4 = elData.qop[3] * (state.dir + h * k3).cross(bField);
      return true;
    }

	/// @brief After a RKN4 step was accepted by the stepper this method has an additional veto on the quality of the step. The veto lies in evaluation of the energy loss and the therewith constrained to keep the momentum after the step in reasonable values.
	///
	/// @tparam stepper_state_t Type of the state of the stepper
	/// @param [in, out] state State of the stepper
	/// @param [in] h Step size
	/// @return Boolean flag if step evaluation is valid
	template<typename stepper_state_t>
    bool
    finalizeStep(stepper_state_t& state, const double h)
    {
		// Evaluate the new momentum
		double newMomentum = state.p + conv * (h / 6.) * (elData.dPds[0] + 2. * (elData.dPds[1] + elData.dPds[2]) + elData.dPds[3]);
	
		// Break propagation if momentum becomes below cut-off
		if (units::Nat2SI<units::MOMENTUM>(newMomentum) < momentumCutOff)
			return false;
		else
		{
			// Update momentum
			state.p = newMomentum;
			return true;
		}
	}

    /// @brief Evaluates the transport matrix D for the jacobian
    ///
    /// @param [in] dir Direction of the particle
    /// @param [in] bFiedl1 B-field at the first position
    /// @param [in] bField2 B-field at the middle position
    /// @param [in] bField3 B-field at the last position
    /// @param [in] h Step size
    /// @param [in] k1 Vector of k1
    /// @param [in] k2 Vector of k2
    /// @param [in] k3 Vector of k3
    /// @param [out] D Transport matrix
    /// @return Boolean flag is step evaluation is valid
    bool
    evaluateD(const Vector3D& dir,
              const Vector3D& bField1,
              const Vector3D& bField2,
              const Vector3D& bField3,
              const double    h,
              const Vector3D& k1,
              const Vector3D& k2,
              const Vector3D& k3,
              ActsMatrixD<7, 7>& D)
    {
      if (!volume || !(*volume) || !(*volume)->material()) return true;

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

      // Evaluation of the rightmost column without the last term.
      elData.jdL[0] = elData.dLdl[0];
      dk1dL         = dir.cross(bField1);
      elData.jdL[1] = elData.dLdl[1] * (1. + half_h * elData.jdL[0]);
      dk2dL = (1. + half_h * elData.jdL[0]) * (dir + half_h * k1).cross(bField2)
          + elData.qop[1] * half_h * dk1dL.cross(bField2);
      elData.jdL[2] = elData.dLdl[2] * (1. + half_h * elData.jdL[1]);
      dk3dL = (1. + half_h * elData.jdL[1]) * (dir + half_h * k2).cross(bField2)
          + elData.qop[2] * half_h * dk2dL.cross(bField2);
      elData.jdL[3] = elData.dLdl[3] * (1. + h * elData.jdL[2]);
      dk4dL         = (1. + h * elData.jdL[2]) * (dir + h * k3).cross(bField3)
          + elData.qop[3] * h * dk3dL.cross(bField3);

      dk1dT(0, 1) = bField1.z();
      dk1dT(0, 2) = -bField1.y();
      dk1dT(1, 0) = -bField1.z();
      dk1dT(1, 2) = bField1.x();
      dk1dT(2, 0) = bField1.y();
      dk1dT(2, 1) = -bField1.x();
      dk1dT *= elData.qop[0];

      dk2dT += half_h * dk1dT;
      dk2dT *= cross(dk2dT, bField2);
      dk2dT *= elData.qop[1];

      dk3dT += half_h * dk2dT;
      dk3dT *= cross(dk3dT, bField2);
      dk3dT *= elData.qop[2];

      dk4dT += h * dk3dT;
      dk4dT *= cross(dk4dT, bField3);
      dk4dT *= elData.qop[3];

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h * h / 6 * (dk1dL + dk2dL + dk3dL);

      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      // Evaluation of the dLambda''/dlambda term
      D(6, 6) += conv * (h / 6.)
          * (elData.jdL[0] + 2. * (elData.jdL[1] + elData.jdL[2])
             + elData.jdL[3]);
      return true;
    }

  private:
    /// Energy loss calculator
    detail::IonisationLoss ionisationLoss;
    detail::RadiationLoss  radiationLoss;

    /// @brief This function calculates the energy loss dE per path length ds of
    /// a
    /// particle through material. The energy loss consists of ionisation and
    /// radiation.
    /// @note The calculations use SI units and it is assumed that the arguments
    /// are given in SI units.
    ///
    /// @tparam material_t Type of the material
    /// @param [in] momentum Initial momentum of the particle
    /// @param [in] energy Initial energy of the particle
    /// @param [in] mass Mass of the particle
    /// @param [in] material Penetrated material
    /// @param [in] pdg PDG code of the particle
    /// @param [in] meanEnergyLoss Boolean flag if mean or mode should be
    /// evaluated for the energy loss
    /// @return Infinitesimal energy loss
    template <typename material_t>
    double
    dEds(const double      momentum,
         const double      energy,
         const material_t& material) const
    {
      // Easy exit if material is invalid
      if (material.X0() == 0 || material.Z() == 0) return 0.;

      // Calculate energy loss by
      // a) ionisation
      double ionisationEnergyLoss
          = ionisationLoss(elData.massSI,
                           momentum * units::_c / energy,
                           energy / (elData.massSI * units::_c2),
                           material,
                           1.,
                           meanEnergyLoss,
                           true)
                .first;
      // b) radiation
      double radiationEnergyLoss
          = radiationLoss(energy, elData.massSI, material, pdg, 1., true);

      // Rescaling for mode evaluation.
      // TODO: Factor just copied from Athena but not tested for correctness
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
    /// @return Derivative evaluated at the point defined by the
    /// function parameters
    template <typename material_t>
    double
    dgdqop(const double      energy,
           const double      qop,
           const material_t& material) const
    {
      // Fast exit if material is invalid
      if (material.X0() == 0. || material.Z() == 0.
          || material.zOverAtimesRho() == 0.)
        return 0.;

      // Bethe-Bloch
      const double betheBlochDerivative
          = ionisationLoss.dqop(energy, qop, elData.massSI, material, true);

      // Bethe-Heitler (+ pair production & photonuclear interaction for muons)
      const double radiationDerivative
          = radiationLoss.dqop(elData.massSI, material, qop, energy, pdg, true);

      // Return the total derivative
      if (meanEnergyLoss)
        return betheBlochDerivative + radiationDerivative;
      else
        // TODO: The scaling factors are just copied from Athena without any
        // test
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
      double E = std::sqrt(elData.initialMomentum * elData.initialMomentum
                               * units::_c2
                           + elData.massSI * elData.massSI * units::_c4);
      // Use the same energy loss throughout the step.
      elData.g = dEds(elData.initialMomentum, E, *(elData.material));
      // Change of the momentum per path length
      // dPds = dPdE * dEds
      elData.dPds[0] = elData.g * E / (elData.initialMomentum * units::_c2);
      if (state.covTransport) {
        // Calculate the change of the energy loss per path length and
        // inverse momentum
        if (includeGgradient) {
          elData.dgdqopValue = dgdqop(
              E,
              elData.qop[0],
              *(elData.material));  // Use this value throughout the step.
        }
        // Calculate term for later error propagation
        elData.dLdl[0]
            = (-elData.qop[0] * elData.qop[0] * elData.g * E
                   * (3.
                      - (elData.initialMomentum * elData.initialMomentum
                         * units::_c2)
                          / (E * E))
               - elData.qop[0] * elData.qop[0] * elData.qop[0] * E
                   * elData.dgdqopValue)
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
      elData.currentMomentum = elData.initialMomentum + h * elData.dPds[i - 1];
      // if (momentum <= momentumCutOff) return false; //Abort propagation
      double E = std::sqrt(elData.currentMomentum * elData.currentMomentum
                               * units::_c2
                           + elData.massSI * elData.massSI * units::_c4);
      elData.dPds[i] = elData.g * E / (elData.currentMomentum * units::_c2);
      elData.qop[i]  = state.q / elData.currentMomentum;
      // Calculate term for later error propagation
      if (state.covTransport) {
        elData.dLdl[i]
            = (-elData.qop[i] * elData.qop[i] * elData.g * E
                   * (3.
                      - (elData.currentMomentum * elData.currentMomentum
                         * units::_c2)
                          / (E * E))
               - elData.qop[i] * elData.qop[i] * elData.qop[i] * E
                   * elData.dgdqopValue)
            / units::_c3;
      }
    }
  };

}  // namespace detail
}  // namespace Acts
