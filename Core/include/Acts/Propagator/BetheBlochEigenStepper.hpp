// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/detail/Constants.hpp"
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Detector/TrackingVolume.hpp"

namespace Acts {

/// @brief Runge-Kutta-Nystroem stepper based on Eigen implementation
/// for the following ODE:
///
/// r = (x,y,z)    ... global position
/// T = (Ax,Ay,Az) ... momentum direction (normalized)
///
/// dr/ds = T
/// dT/ds = q/p * (T x B)
///
/// with s being the arc length of the track, q the charge of the particle,
/// p its momentum and B the magnetic field
///
template <typename BField, typename corrector_t = VoidCorrector>
class BetheBlochEigenStepper : public EigenStepper<BField, corrector_t>
{

public:
  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State : public EigenStepper<BField, corrector_t>::State
  {
    /// Constructor from the initial track parameters
    /// @param[in] par The track parameters at start
    /// @param[in] ndir The navigation direciton w.r.t momentum
    /// @param[in] sszice is the maximum step size
    ///
    /// @note the covariance matrix is copied when needed
    template <typename T>
    explicit State(const T&            par,
                   NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max())
      : EigenStepper<BField, corrector_t>::State(par, ndir, ssize)
    {
    }

    /// Mass
    double* mass = nullptr;
    
    /// PDG code
    int* pdg = nullptr;

    /// Volume with material that is passed
    TrackingVolume const** volume = nullptr;

    /// Boolean flag for energy loss while stepping
    bool energyLossFlag = true;
    
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
  };

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using state_type = BetheBlochEigenStepper::State;

  /// Constructor requires knowledge of the detector's magnetic field
  BetheBlochEigenStepper(BField bField = BField())
    : m_bField(std::move(bField)){};

private:
  /// @brief This function calculates the energy loss dE per path length ds of a
  /// particle through material. The energy loss consists of ionisation and
  /// radiation.
  /// @note The calculations use SI units and it is assumed that the arguments are given in SI units.
  ///
  /// @tparam material_t Type of the material
  /// @param [in] momentum Initial momentum of the particle
  /// @param [in] energy Initial energy of the particle
  /// @param [in] mass Mass of the particle
  /// @param [in] material Penetrated material
  /// @param [in] pdg PDG code of the particle
  /// @param [in] meanEnergyLoss Boolean flag if mean or mode should be evaluated for the energy loss due to ionisation
  /// @return Infinitesimal energy loss
  template <typename material_t>
  double
  dEds(const double      momentum,
       const double      energy,
       const double      mass,
       const material_t& material,
       const int pdg,
       const bool meanEnergyLoss = true) const
  {
    // Easy exit if material is invalid
    if (material.X0() == 0 || material.Z() == 0) return 0.;

    // Calculate energy loss by
    // a) ionisation
    // TODO: make return as pure double
    double ionisationEnergyLoss = ionisationLoss(mass,
                                             momentum * units::_c / energy,
                                             energy / (mass * units::_c2),
                                             material,
                                             1.,
                                             meanEnergyLoss,
                                             true)
                                      .first;
    // b) radiation
    double radiationEnergyLoss = radiationLoss(energy, mass, material, pdg, 1., true);
               
    // Rescaling for mode evaluation. 
    // TODO: Factor just copied from Athena but not tested for correctness
    if(!meanEnergyLoss)
		radiationEnergyLoss *= 0.15;

    // Return sum of contributions
    return ionisationEnergyLoss + radiationEnergyLoss;
  }

  /// @brief This function calculates the derivation of g=dE/dx by d(q/p)
  ///
  /// @tparam material_t Type of the material
  /// @param [in] energy Initial energy of the particle
  /// @param [in] qop Initial value of q/p of the particle
  /// @param [in] mass Mass of the particle
  /// @param [in] material Penetrated material
  /// @return Derivative evaEnergyLossDataluated at the point defined by the
  /// function
  /// parameters
  template <typename material_t>
  double
  dgdqop(const double      energy,
         const double      qop,
         const double      mass,
         const material_t& material) const
  {
	  // TODO: all unit conversions from nat to SI
	  // TODO: camelcase naming
	  // TODO: Move derivatives to InteractionFormulas.hpp
    // Fast exit if material is invalid
    if (material.X0() == 0 || material.Z() == 0
        || material.zOverAtimesRho() == 0)
      return 0.;

    // Bethe-Bloch
    // Constants for readability
    const double me    = constants::me;
    const double me2   = me * me;
    const double qopNU = qop / units::SI2Nat<units::MOMENTUM>(1.);
    const double qop3  = qopNU * qopNU * qopNU;
    const double qop4  = qop3 * qopNU;
    const double m     = units::SI2Nat<units::MASS>(mass);
    const double m2    = m * m;
    const double m4    = m2 * m2;
    const double I     = constants::eionisation * std::pow(material.Z(), 0.9);
    const double I2    = I * I;
    const double E     = units::SI2Nat<units::ENERGY>(energy);
    const double gamma = E / m;
    const double beta  = std::abs(1 / (E * qopNU));
    const double beta2 = beta * beta;
    const double kaz
        = 0.5 * units::Nat2SI<units::ENERGY>(30.7075 * units::_MeV)
            * units::_mm2 * units::_e2 / units::_g * material.zOverAtimesRho();
    // Parts of the derivative
    double lnCore
        = 4. * me2 / (m4 * I2 * qop4) / (1. + 2. * gamma * me / m + me2 / m2);
    double lnCore_deriv = -4. * me2 / (m4 * I2)
        * std::pow(qop4 + 2. * gamma * qop4 * me / m + qop4 * me2 / m2, -2.)
        * (4. * qop3 + 8. * me * qop3 * gamma / m
           - 2. * me * qopNU / (m2 * m * gamma)
           + 4. * qop3 * me2 / m2);

    // Combine parts
    double ln_deriv
        = 2. * qopNU * m2 * std::log(lnCore) + lnCore_deriv / (lnCore * beta2);
    double Bethe_Bloch_deriv = -kaz * ln_deriv * units::Nat2SI<units::MOMENTUM>(1.);
    
    // Density effect, only valid for high energies (gamma > 10 -> p > 1GeV for
    // muons)
    if (gamma > 10.) {
      double delta
          = 2. * std::log(constants::eplasma
                          * std::sqrt(1000. * material.zOverAtimesRho())
                          / I)
          + 2. * std::log(beta * gamma) - 1.;
      double delta_deriv = -2. / (qopNU * beta2) + 2. * delta * qopNU * m2;
      Bethe_Bloch_deriv += kaz * delta_deriv * units::Nat2SI<units::MOMENTUM>(1.);
    }

    // Bethe-Heitler
    double Bethe_Heitler_deriv = me2 / (m2 * material.X0() * qop * qop * qop * energy);

    // Radiative corrections (e+e- pair production + photonuclear) for muons at
    // energies above 8 GeV and below 1 TeV
    // TODO: no dgdqop for radiation if there is no radiation
    double radiative_deriv = 0.;

    // Return the total derivative
    return Bethe_Bloch_deriv + Bethe_Heitler_deriv + radiative_deriv;
  }

	/// @brief This struct serves as data container to keep track of all parameters that are related to an energy loss of a particle in matter.
  struct EnergyLossData
  {
	  /// Particles momentum at k1
    double                          initialMomentum;
    /// Particles mass in SI units
    double                          massSI;
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
    double g;
  };

  /// @brief Initializer of all parameters related to a RKN4 step with energy
  /// loss of a particle in material
  /// @note This function serves for reducing the length of the actual
  /// BetheBlochEigenStepper::step() function.
  ///
  /// @param [in, out] eld Data container related to material interaction of
  /// the particle
  /// @param [in] state Deliverer of configurations
  /// @return Modified EnergyLossData
  std::unique_ptr<EnergyLossData>
  initializeEnergyLoss(std::unique_ptr<EnergyLossData> eld,
						const State& state) const
  {
    double E = std::sqrt(eld->initialMomentum * eld->initialMomentum * units::_c2
                         + eld->massSI * eld->massSI * units::_c4);
    // Use the same energy loss throughout the step.
    eld->g
        = dEds(eld->initialMomentum, E, eld->massSI, *(eld->material), *(state.pdg), state.meanEnergyLoss);
    // Change of the momentum per path length
    // dPds = dPdE * dEds
    eld->dPds[0] = eld->g * E / (eld->initialMomentum * units::_c2);
    if (state.covTransport) {
      // Calculate the change of the the energy loss per path length and
      // inverse momentum
      if (state.includeGgradient) {
        eld->dgdqopValue = dgdqop(
            E,
            eld->qop[0],
            eld->massSI,
            *(eld->material));  // Use this value throughout the step.
      }
      // Calculate term for later error propagation
      eld->dLdl[0] = (-eld->qop[0] * eld->qop[0] * eld->g * E
              * (3.
                 - (eld->initialMomentum * eld->initialMomentum * units::_c2)
                     / (E * E))
          - eld->qop[0] * eld->qop[0] * eld->qop[0] * E
              * eld->dgdqopValue) / units::_c3;
    }
    return std::move(eld);
  }

  /// @brief Update of the kinematic parameters of the RKN4 sub-steps after
  /// initialization with energy loss of a particle in material
  /// @note This function serves for reducing the length of the actual
  /// BetheBlochEigenStepper::step() function.
  ///
  /// @param [in] eld Data container related to material interaction of
  /// the particle
  /// @param [out] momentum Updated momentum
  /// @param [in] h Stepped distance of the sub-step (1-3)
  /// @param [in] q Charge of the particle
  /// @param [in] covTransport Boolean flag if covariance should be transported
  /// @param [in] i Index of the sub-step (1-3)
  /// @return Modified EnergyLossData
  std::unique_ptr<EnergyLossData>
  updateEnergyLoss(std::unique_ptr<EnergyLossData> eld,
                   double&         momentum,
                   const double    h,
                   const int       q,
                   const bool      covTransport,
                   const int       i) const
  {
    // Update parameters related to a changed momentum
    momentum = eld->initialMomentum + h * eld->dPds[i - 1];
    // if (momentum <= momentumCutOff) return false; //Abort propagation
    double E = std::sqrt(momentum * momentum * units::_c2
                         + eld->massSI * eld->massSI * units::_c4);
    eld->dPds[i] = eld->g * E / (momentum * units::_c2);
    eld->qop[i]  = q / momentum;
    // Calculate term for later error propagation
    if (covTransport) {
      eld->dLdl[i] = (-eld->qop[i] * eld->qop[i] * eld->g
              * E * (3. - (momentum * momentum * units::_c2) / (E * E))
          - eld->qop[i] * eld->qop[i] * eld->qop[i] * E
              * eld->dgdqopValue) / units::_c3;
    }
    return std::move(eld);
  }

// TODO: using B field gradient

public:
  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                      parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  double
  step(State& state) const
  {

    std::unique_ptr<EnergyLossData> elData;
    double momentum, qop0;

	// Set up initial data
    if (state.energyLossFlag && (*state.volume) && (*state.volume)->material()) {
		// Set up container for energy loss
      elData           = std::make_unique<EnergyLossData>();
      elData->massSI   = units::Nat2SI<units::MASS>(*state.mass);
      elData->material = (*state.volume)->material();
      elData->initialMomentum = units::Nat2SI<units::MOMENTUM>(state.p);
      momentum                = elData->initialMomentum;
      elData->qop[0]          = state.q / momentum;
    } else {
		// Set up data for the case without energy loss
      momentum = units::Nat2SI<units::MOMENTUM>(state.p);
      qop0     = state.q / momentum;
    }
    // Runge-Kutta integrator state
    double   h2, half_h;
    Vector3D B_middle, B_last, k2, k3, k4;

    // First Runge-Kutta point (at current position)
    const Vector3D B_first = this->getField(state, state.pos);
    const Vector3D k1      = (elData ? elData->qop[0] : qop0)
        * state.dir.cross(
              B_first);

    // Calculate the energy loss
    if (elData) {
      elData = initializeEnergyLoss(std::move(elData), state);
    }

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the local
    // integration error. The results are stated in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory
    const auto tryRungeKuttaStep = [&](const double h) -> double {

      // State the square and half of the step size
      h2     = h * h;
      half_h = h * 0.5;

      // Second Runge-Kutta point
      if (elData) {
		  // Update parameters and check for momentum condition
        elData = updateEnergyLoss(
            std::move(elData), momentum, h * 0.5, state.q, state.covTransport, 1);
        if(momentum < state.momentumCutOff)
			return std::numeric_limits<double>::max();
      }

      const Vector3D pos1 = state.pos + half_h * state.dir + h2 * 0.125 * k1;
      B_middle            = this->getField(state, pos1);
      k2                  = (elData ? elData->qop[1] : qop0)
          * (state.dir + half_h * k1).cross(B_middle);

      // Third Runge-Kutta point
      if (elData) {
		  // Update parameters and check for momentum condition
        elData = updateEnergyLoss(
            std::move(elData), momentum, h * 0.5, state.q, state.covTransport, 2);
       if(momentum < state.momentumCutOff)
			return std::numeric_limits<double>::max();
      }

      k3 = (elData ? elData->qop[2] : qop0)
          * (state.dir + half_h * k2).cross(B_middle);

      // Last Runge-Kutta point
      if (elData) {
		  // Update parameters and check for momentum condition
        elData = updateEnergyLoss(std::move(elData), momentum, h, state.q, state.covTransport, 3);
        if(momentum < state.momentumCutOff)
			return std::numeric_limits<double>::max();
      }

      const Vector3D pos2 = state.pos + h * state.dir + h2 * 0.5 * k3;
      B_last              = this->getField(state, pos2);
      k4                  = (elData ? elData->qop[3] : qop0)
          * (state.dir + h * k3).cross(B_last);
      // Return an estimate of the local integration error
      return h2 * (k1 - k2 - k3 + k4).template lpNorm<1>();
    };

    // Select and adjust the appropriate Runge-Kutta step size in ATLAS style
    double error_estimate = std::max(tryRungeKuttaStep(state.stepSize), 1e-20);
    while (error_estimate > 4. * state.tolerance) {
      state.stepSize = state.stepSize
          * std::min(std::max(
                         0.25,
                         std::pow((state.tolerance / error_estimate), 0.25)),
                     4.);
      // If step size becomes too small the particle remains at the initial place
      if (state.stepSize < state.stepSizeCutOff) {
        return 0.; // Not moving due to too low momentum needs an aborter
      }
      error_estimate = std::max(tryRungeKuttaStep(state.stepSize), 1e-20);
    }

    // use the adjusted step size
    const double h = state.stepSize;
    
    // Break propagation if momentum becomes below cut-off
    if (elData) {
		double newMomentum = state.p + units::SI2Nat<units::MOMENTUM>(
          (h / 6.) * (elData->dPds[0] + 2. * (elData->dPds[1] + elData->dPds[2])
                      + elData->dPds[3]));
        if(units::Nat2SI<units::MOMENTUM>(newMomentum) < state.momentumCutOff)
			return 0.;
		else
			// Update momentum
			state.p = newMomentum;
	}

    // When doing error propagation, update the associated Jacobian matrix
    if (state.covTransport) {
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
      /// The evaluation is based by propagating the parameters T and lambda
      /// (including g(lambda) and E(lambda)) as given in eq. 16 and evaluating
      /// the derivations for matrix A.

      // The step transport matrix in global coordinates
      ActsMatrixD<7, 7> D = ActsMatrixD<7, 7>::Identity();
      const double conv = units::SI2Nat<units::MOMENTUM>(1);

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
      if (elData) {
		  // For the case of energy loss
        elData->jdL[0] = elData->dLdl[0];
        dk1dL          = state.dir.cross(B_first);
        elData->jdL[1] = elData->dLdl[1] * (1. + half_h * elData->jdL[0]);
        dk2dL          = (1. + half_h * elData->jdL[0])
                * (state.dir + half_h * k1).cross(B_middle)
            + elData->qop[1] * half_h * dk1dL.cross(B_middle);
        elData->jdL[2] = elData->dLdl[2] * (1. + half_h * elData->jdL[1]);
        dk3dL          = (1. + half_h * elData->jdL[1])
                * (state.dir + half_h * k2).cross(B_middle)
            + elData->qop[2] * half_h * dk2dL.cross(B_middle);
        elData->jdL[3] = elData->dLdl[3] * (1. + h * elData->jdL[2]);
        dk4dL = (1. + h * elData->jdL[2]) * (state.dir + h * k3).cross(B_last)
            + elData->qop[3] * h * dk3dL.cross(B_last);
      } else {
		  // For the case without energy loss
        dk1dL = state.dir.cross(B_first);
        dk2dL = (state.dir + half_h * k1).cross(B_middle)
            + qop0 * half_h * dk1dL.cross(B_middle);
        dk3dL = (state.dir + half_h * k2).cross(B_middle)
            + qop0 * half_h * dk2dL.cross(B_middle);
        dk4dL = (state.dir + h * k3).cross(B_last)
            + qop0 * h * dk3dL.cross(B_last);
      }

      dk1dT(0, 1) = B_first.z();
      dk1dT(0, 2) = -B_first.y();
      dk1dT(1, 0) = -B_first.z();
      dk1dT(1, 2) = B_first.x();
      dk1dT(2, 0) = B_first.y();
      dk1dT(2, 1) = -B_first.x();
      dk1dT *= (elData ? elData->qop[0] : qop0);

      dk2dT += h / 2 * dk1dT;
      dk2dT *= cross(dk2dT, B_middle);
      dk2dT *= (elData ? elData->qop[1] : qop0);

      dk3dT += h / 2 * dk2dT;
      dk3dT *= cross(dk3dT, B_middle);
      dk3dT *= (elData ? elData->qop[2] : qop0);

      dk4dT += h * dk3dT;
      dk4dT *= cross(dk4dT, B_last);
      dk4dT *= (elData ? elData->qop[3] : qop0);

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h2 / 6 * (dk1dL + dk2dL + dk3dL);

      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      // Evaluation of the dLambda''/dlambda term
      if (elData)
        D(6, 6) += conv * (h / 6.)
            * (elData->jdL[0] + 2. * (elData->jdL[1] + elData->jdL[2])
               + elData->jdL[3]);
      std::cout << "D:\n" << D << std::endl;
      std::cout << "jac:\n" << state.jacTransport << std::endl;
      // for moment, only update the transport part
      state.jacTransport = D * state.jacTransport;
    }

    // Update the track parameters according to the equations of motion
    state.pos += h * state.dir + h2 / 6. * (k1 + k2 + k3);
    state.dir += h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
    state.dir /= state.dir.norm();
    state.derivative.template head<3>()     = state.dir;
    state.derivative.template segment<3>(3) = k4;

    std::cout << "result pos: " << state.pos << std::endl;
    std::cout << "result dir: " << state.dir << std::endl;
    std::cout << "result p: " << state.p << std::endl;
    std::cout << "result cov:\n" << state.jacTransport << std::endl;
    state.pathAccumulated += h;
    std::exit(1);
    return h;
  }

private:
  /// Magnetic field inside of the detector
  BField m_bField;
	/// Energy loss calculator
  detail::IonisationLoss ionisationLoss;
  detail::RadiationLoss radiationLoss;
};

}  // namespace Acts
