// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/detail/Constants.hpp"
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"

// tmp includes
#include <cstdlib>
#include "Acts/Material/Material.hpp"

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
class BetheBlochEigenStepper
{

private:
  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    using type = BoundParameters;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    using type = CurvilinearParameters;
  };

public:
  using cstep = detail::ConstrainedStep;

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State
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
      : pos(par.position())
      , dir(par.momentum().normalized())
      , p(par.momentum().norm())
      , q(par.charge())
      , navDir(ndir)
      , stepSize(ndir * std::abs(ssize))
    {
      // remember the start parameters
      startPos = pos;
      startDir = dir;
      // Init the jacobian matrix if needed
      if (par.covariance()) {
        // Get the reference surface for navigation
        const auto& surface = par.referenceSurface();
        // set the covariance transport flag to true and copy
        covTransport = true;
        cov          = ActsSymMatrixD<5>(*par.covariance());
        surface.initJacobianToGlobal(jacToGlobal, pos, dir, par.parameters());
      }
    }

    /// Global particle position accessor
    Vector3D
    position() const
    {
      return pos;
    }

    /// Momentum direction accessor
    Vector3D
    direction() const
    {
      return dir;
    }

    /// Actual momentum accessor
    Vector3D
    momentum() const
    {
      return p * dir;
    }

    /// Charge access
    double
    charge() const
    {
      return q;
    }

    /// Return a corrector
    corrector_t
    corrector() const
    {
      return corrector_t(startPos, startDir, pathAccumulated);
    }

    /// Method to update momentum, direction and p
    ///
    /// @param uposition the updated position
    /// @param udirection the updated direction
    /// @param p the updated momentum value
    void
    update(const Vector3D& uposition, const Vector3D& udirection, double up)
    {
      pos = uposition;
      dir = udirection;
      p   = up;
    }

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the state
    ///
    /// @param reinitialize is a flag to steer whether the
    ///        state should be reinitialized at the new
    ///        position
    ///
    /// @return the full transport jacobian
    const ActsMatrixD<5, 5>
    covarianceTransport(bool reinitialize = false)
    {
      // Optimized trigonometry on the propagation direction
      const double x = dir(0);  // == cos(phi) * sin(theta)
      const double y = dir(1);  // == sin(phi) * sin(theta)
      const double z = dir(2);  // == cos(theta)
      // can be turned into cosine/sine
      const double cosTheta    = z;
      const double sinTheta    = sqrt(x * x + y * y);
      const double invSinTheta = 1. / sinTheta;
      const double cosPhi      = x * invSinTheta;
      const double sinPhi      = y * invSinTheta;
      // prepare the jacobian to curvilinear
      ActsMatrixD<5, 7> jacToCurv = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
        // We normally operate in curvilinear coordinates defined as follows
        jacToCurv(0, 0) = -sinPhi;
        jacToCurv(0, 1) = cosPhi;
        jacToCurv(1, 0) = -cosPhi * cosTheta;
        jacToCurv(1, 1) = -sinPhi * cosTheta;
        jacToCurv(1, 2) = sinTheta;
      } else {
        // Under grazing incidence to z, the above coordinate system definition
        // becomes numerically unstable, and we need to switch to another one
        const double c    = sqrt(y * y + z * z);
        const double invC = 1. / c;
        jacToCurv(0, 1) = -z * invC;
        jacToCurv(0, 2) = y * invC;
        jacToCurv(1, 0) = c;
        jacToCurv(1, 1) = -x * y * invC;
        jacToCurv(1, 2) = -x * z * invC;
      }
      // Directional and momentum parameters for curvilinear
      jacToCurv(2, 3) = -sinPhi * invSinTheta;
      jacToCurv(2, 4) = cosPhi * invSinTheta;
      jacToCurv(3, 5) = -invSinTheta;
      jacToCurv(4, 6) = 1;
      // Apply the transport from the steps on the jacobian
      jacToGlobal = jacTransport * jacToGlobal;
      // Transport the covariance
      ActsRowVectorD<3>       normVec(dir);
      const ActsRowVectorD<5> sfactors
          = normVec * jacToGlobal.topLeftCorner<3, 5>();
      // The full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jacFull
          = jacToCurv * (jacToGlobal - derivative * sfactors);
      // Apply the actual covariance transport
      cov = (jacFull * cov * jacFull.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jacToGlobal  = ActsMatrixD<7, 5>::Zero();
        jacTransport = ActsMatrixD<7, 7>::Identity();
        // fill the jacobian to global for next transport
        jacToGlobal(0, eLOC_0) = -sinPhi;
        jacToGlobal(0, eLOC_1) = -cosPhi * cosTheta;
        jacToGlobal(1, eLOC_0) = cosPhi;
        jacToGlobal(1, eLOC_1) = -sinPhi * cosTheta;
        jacToGlobal(2, eLOC_1) = sinTheta;
        jacToGlobal(3, ePHI)   = -sinTheta * sinPhi;
        jacToGlobal(3, eTHETA) = cosTheta * cosPhi;
        jacToGlobal(4, ePHI)   = sinTheta * cosPhi;
        jacToGlobal(4, eTHETA) = cosTheta * sinPhi;
        jacToGlobal(5, eTHETA) = -sinTheta;
        jacToGlobal(6, eQOP)   = 1;
      }
      // return the full transport jacobian
      return jacFull;
    }

    /// Method for on-demand transport of the covariance
    /// to a new curvilinear frame at current  position,
    /// or direction of the state
    ///
    /// @tparam S the Surfac type
    ///
    /// @param surface is the surface to which the covariance is
    ///        forwarded to
    /// @param reinitialize is a flag to steer whether the
    ///        state should be reinitialized at the new
    ///        position
    /// @note no check is done if the position is actually on the surface
    ///
    /// @return the full transport jacobian
    template <typename S>
    const ActsMatrixD<5, 5>
    covarianceTransport(const S& surface, bool reinitialize = false)
    {
      using VectorHelpers::phi;
      using VectorHelpers::theta;
      // Initialize the transport final frame jacobian
      ActsMatrixD<5, 7> jacToLocal = ActsMatrixD<5, 7>::Zero();
      // initalize the jacobian to local, returns the transposed ref frame
      auto rframeT = surface.initJacobianToLocal(jacToLocal, pos, dir);
      // Update the jacobian with the transport from the steps
      jacToGlobal = jacTransport * jacToGlobal;
      // calculate the form factors for the derivatives
      const ActsRowVectorD<5> sVec
          = surface.derivativeFactors(pos, dir, rframeT, jacToGlobal);
      // the full jacobian is ([to local] jacobian) * ([transport] jacobian)
      const ActsMatrixD<5, 5> jacFull
          = jacToLocal * (jacToGlobal - derivative * sVec);
      // Apply the actual covariance transport
      cov = (jacFull * cov * jacFull.transpose());
      // Reinitialize if asked to do so
      // this is useful for interruption calls
      if (reinitialize) {
        // reset the jacobians
        jacToGlobal  = ActsMatrixD<7, 5>::Zero();
        jacTransport = ActsMatrixD<7, 7>::Identity();
        // reset the derivative
        derivative = ActsVectorD<7>::Zero();
        // fill the jacobian to global for next transport
        Vector2D loc{0., 0.};
        surface.globalToLocal(pos, dir, loc);
        ActsVectorD<5> pars;
        pars << loc[eLOC_0], loc[eLOC_1], phi(dir), theta(dir), q / p;
        surface.initJacobianToGlobal(jacToGlobal, pos, dir, pars);
      }
      // store in the global jacobian
      jacobian = jacFull * jacobian;
      // return the full transport jacobian
      return jacFull;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);
    /// Global start particle position
    Vector3D startPos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);
    /// Momentum start direction (normalized)
    Vector3D startDir = Vector3D(1, 0, 0);

    /// Momentum
    double p = 0.;

    /// The charge
    double q = 1.;

    /// Mass
    double mass = 0.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// The full jacobian of the transport
    ActsMatrixD<5, 5> jacobian = ActsMatrixD<5, 5>::Identity();

    /// Jacobian from local to the global frame
    ActsMatrixD<7, 5> jacToGlobal = ActsMatrixD<7, 5>::Zero();

    /// Pure transport jacobian part from runge kutta integration
    ActsMatrixD<7, 7> jacTransport = ActsMatrixD<7, 7>::Identity();

    /// The propagation derivative
    ActsVectorD<7> derivative = ActsVectorD<7>::Zero();

    /// Covariance matrix (and indicator)
    //// associated with the initial error on track parameters
    bool              covTransport = false;
    ActsSymMatrixD<5> cov          = ActsSymMatrixD<5>::Zero();

    /// Lazily initialized state of the field Cache
    bool fieldCacheReady = false;

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    concept::AnyFieldCell<> fieldCache;

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
    cstep stepSize = std::numeric_limits<double>::max();
  };

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using state_type = State;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename T>
  using step_parameter_type = CurvilinearParameters;
  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we usually return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor requires knowledge of the detector's magnetic field
  BetheBlochEigenStepper(BField bField = BField())
    : m_bField(std::move(bField)){};

  /// Convert the propagation state (global) to curvilinear parameters
  /// @param state The stepper state
  /// @param reinitialize is a flag to (optionally) reinitialse the state
  /// @return curvilinear parameters
  static CurvilinearParameters
  convert(State& state, bool reinitialize = false)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    // only do the covariance transport if needed
    if (state.covTransport) {
      // transport the covariance forward
      state.covarianceTransport(reinitialize);
      covPtr = std::make_unique<const ActsMatrixD<5, 5>>(state.cov);
    }
    // return the parameters
    return CurvilinearParameters(
        std::move(covPtr), state.pos, state.p * state.dir, state.q);
  }

  /// Convert the propagation state to track parameters at a certain surface
  ///
  /// @tparam S The surface type
  ///
  /// @param [in] state Propagation state used
  /// @param [in] surface Destination surface to which the conversion is done
  template <typename S>
  static BoundParameters
  convert(State& state, const S& surface, bool reinitialize = false)
  {
    std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
    // Perform error propagation if an initial covariance matrix was provided
    if (state.covTransport) {
      // transport the covariance forward
      state.covarianceTransport(surface, reinitialize);
      covPtr = std::make_unique<const ActsSymMatrixD<5>>(state.cov);
    }
    // return the bound parameters
    return BoundParameters(
        std::move(covPtr), state.pos, state.p * state.dir, state.q, surface);
  }

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(State& state, const Vector3D& pos) const
  {
    if (!state.fieldCacheReady || !state.fieldCache.isInside(pos)) {
      state.fieldCacheReady = true;
      state.fieldCache      = m_bField.getFieldCell(pos);
    }
    // get the field from the cell
    return std::move(state.fieldCache.getField(pos));
  }

  /// @brief This function calculates the energy loss dE per path length ds of a
  /// particle through material. The energy loss consists of ionisation and
  /// radiation.
  ///
  /// @tparam material_t Type of the material
  /// @param [in] momentum Initial momentum of the particle
  /// @param [in] energy Initial energy of the particle
  /// @param [in] mass Mass of the particle
  /// @param [in] material Penetrated material
  /// @return Infinitesimal energy loss
  template <typename material_t>
  double
  dEds(const double      momentum,
       const double      energy,
       const double      mass,
       const material_t& material) const
  {  
    // Easy exit if material is invalid
    if (material.X0() == 0 || material.Z() == 0) return 0.;

    // Calculate energy loss by
    // a) ionisation
    // TODO: Allow change between mean and mode & make return as pure double
    double ionisationEnergyLoss
        = energyLoss(mass, momentum * units::_c / energy, energy / (mass * c2), material, 1., true, true).first;
std::cout << "ion loss: " << energy << "\t" << ionisationEnergyLoss << "\t" << units::SI2Nat<units::ENERGY>(ionisationEnergyLoss) << std::endl;
    // b) radiation
    // TODO: There doesn't radiate anything
    double radiationEnergyLoss = 0.;

    // return sum of contributions
    return ionisationEnergyLoss + radiationEnergyLoss;
  }

  /// @brief This function calculates the derivation of g=dE/dx by d(q/p)
  ///
  /// @tparam material_t Type of the material
  /// @param [in] energy Initial energy of the particle
  /// @param [in] qop Initial value of q/p of the particle
  /// @param [in] mass Mass of the particle
  /// @param [in] material Penetrated material
  /// @return Derivative evaluated at the point defined by the function
  /// parameters
  template <typename material_t>
  double
  dgdqop(const double      energy,
         const double      qop,
         const double      mass,
         const material_t& material) const
  {
    // Fast exit if material is invalid
    if (material.X0() == 0 || material.Z() == 0
        || material.zOverAtimesRho() == 0)
      return 0.;

    // Bethe-Bloch
    // Constants for readability
    const double me    = constants::me;
    const double me2   = me * me;
    const double qop3  = qop * qop * qop;
    const double qop4  = qop3 * qop;
    const double m     = mass;
    const double m2    = m * m;
    const double m4    = m2 * m2;
    const double I     = constants::eionisation * std::pow(material.Z(), 0.9);
    const double I2    = I * I;
    const double E     = energy;
    const double gamma = E / m;
    const double beta  = std::abs(1 / (E * qop));
    const double beta2 = beta * beta;
    const double kaz
        = 0.5 * constants::ka_BetheBloch * material.zOverAtimesRho();

    // Parts of the derivative
    double lnCore
        = 4. * me2 / (m4 * I2 * qop4) / (1. + 2. * gamma * me / m + me2 / m2);
    double lnCore_deriv = -4. * me2 / (m4 * I2)
        * std::pow(qop4 + 2. * gamma * qop4 * me / m + qop4 * me2 / m2, -2.)
        * (4. * qop3 + 8. * me * qop3 * gamma / m
           - 2. * me * qop / (m2 * m * gamma)
           + 4. * qop3 * me2 / m2);

    // Combine parts
    double ln_deriv
        = 2. * qop * m2 * std::log(lnCore) + lnCore_deriv / (lnCore * beta2);
    double Bethe_Bloch_deriv = -kaz * ln_deriv;

    // Density effect, only valid for high energies (gamma > 10 -> p > 1GeV for
    // muons)
    if (gamma > 10.) {
      double delta
          = 2. * std::log(28.816e-6
                          * std::sqrt(1000. * material.zOverAtimesRho())
                          / I)
          + 2. * std::log(beta * gamma) - 1.;
      double delta_deriv = -2. / (qop * beta2) + 2. * delta * qop * m2;
      Bethe_Bloch_deriv += kaz * delta_deriv;
    }

    // Bethe-Heitler
    double Bethe_Heitler_deriv = me2 / (m2 * material.X0() * qop3 * E);

    // Radiative corrections (e+e- pair production + photonuclear) for muons at
    // energies above 8 GeV and below 1 TeV
    // TODO: no dgdqop for radiation if there is no radiation
    double radiative_deriv = 0.;
    // if ((m_particle == Trk::muon) && (E > 8000.)) {
    // if (E < 1.e6) {
    // radiative_deriv = 6.803e-5/(m_material->x0()*l*l*l*E) +
    //~ //2.*2.278e-11/(m_material->x0()*l*l*l) -
    // 3.*9.899e-18*E/(m_material->x0()*l*l*l);
    //} else {
    // radiative_deriv = 9.253e-5/(m_material->x0()*l*l*l*E);
    //}
    //}

    // Return the total derivative
    return Bethe_Bloch_deriv + Bethe_Heitler_deriv + radiative_deriv;
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param[in,out] state is the propagation state associated with the track
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
    bool     energyLossFlag = true;
    Material material(
        352.8, 394.133, 9.012, 4., 1.848e-3);         // TODO: how to pass this?
        //~ 95.7, 465.2, 28.03, 14., 2.32e-3);         // TODO: how to pass this?
    state.mass            = 139.57018 * units::_MeV;  // TODO: hardcoded
    double massSI            = units::Nat2SI<units::MASS>(state.mass);  // TODO: hardcoded
    bool includeGgradient = false;
    std::array<double, 4> dL, qop, dP;
    double dgdqopValue    = 0.;
    //~ double momentumCutOff = 0.;
    double g;

    // Charge-momentum ratio, in SI units
    double initialMomentum = units::Nat2SI<units::MOMENTUM>(state.p);
    std::cout << "step " << state.p << "\t" << initialMomentum << std::endl;
    double momentum = initialMomentum;
    qop[0]          = state.q / momentum;

    // Runge-Kutta integrator state
    double   h2, half_h;
    Vector3D B_middle, B_last, k2, k3, k4;

    // First Runge-Kutta point (at current position)
    const Vector3D B_first = getField(state, state.pos);
    const Vector3D k1      = qop[0]
        * state.dir.cross(
              B_first);  // TODO: athena multiplies c to that expression

    // Calculate the energy loss
    if (energyLossFlag && material) {
      double E = std::sqrt(momentum * momentum * c2 + massSI * massSI * c4);
      // Use the same energy loss throughout the step.
      g = dEds(momentum, E, massSI, material);
      // Change of the momentum per path length
      dP[0] = g * E / (momentum * c2);
      if (state.covTransport) {
        // Calculate the change of the the energy loss per path length and
        // inverse momentum
        if (includeGgradient) {
          dgdqopValue
              = dgdqop(E,
                       qop[0],
                       state.mass,
                       material);  // Use this value throughout the step.
        }
        // Calculate term for later error propagation
        dL[0]
            = -qop[0] * qop[0] * g * E * (3. - (momentum * momentum * c2) / (E * E)) / c2
            - qop[0] * qop[0] * qop[0] * E * dgdqopValue;
      }
      std::cout << "dP[0]: " << dP[0] << "\t" << dL[0] << "\t" << (momentum * momentum * c2) / (E * E) << std::endl;
    }


    // TODO: abort condition for too many steps, too low momentum (at any point
    // in the propagation)

    // The following functor starts to perform a Runge-Kutta step of a certain
    // size, going up to the point where it can return an estimate of the local
    // integration error. The results are stated in the local variables above,
    // allowing integration to continue once the error is deemed satisfactory
    const auto tryRungeKuttaStep = [&](const double h) -> double {
		
      // State the square and half of the step size
      h2     = h * h;
      half_h = h * 0.5;

      // Second Runge-Kutta point
      if (energyLossFlag && material) {
        // Update parameters related to a changed momentum
        momentum = initialMomentum + h * 0.5 * dP[0];
        // if (momentum <= momentumCutOff) return false; //Abort propagation
        double E = std::sqrt(momentum * momentum * c2 + massSI * massSI * c4);
        dP[1]    = g * E / (momentum * c2);
        qop[1]   = state.q / momentum;
        // Calculate term for later error propagation
        if (state.covTransport) {
          dL[1] = -qop[1] * qop[1] * g * E
                  * (3. - (momentum * momentum * c2) / (E * E)) / c2
              - qop[1] * qop[1] * qop[1] * E * dgdqopValue;
        }
      }

      const Vector3D pos1 = state.pos + half_h * state.dir + h2 * 0.125 * k1;
      B_middle            = getField(state, pos1);
      k2                  = qop[1] * (state.dir + half_h * k1).cross(B_middle);

      // Third Runge-Kutta point
      if (energyLossFlag && material) {
        // Update parameters related to a changed momentum
        momentum = initialMomentum + h * 0.5 * dP[1];
        // if (momentum <= momentumCutOff) return false; //Abort propagation
        double E = std::sqrt(momentum * momentum * c2 + massSI * massSI * c4);
        dP[2]    = g * E / (momentum * c2);
        qop[2]   = state.q / momentum;
        // Calculate term for later error propagation
        if (state.covTransport) {
          dL[2] = -qop[2] * qop[2] * g * E
                  * (3. - (momentum * momentum * c2) / (E * E)) / c2
              - qop[2] * qop[2] * qop[2] * E * dgdqopValue;
        }
      }

      k3 = qop[2] * (state.dir + half_h * k2).cross(B_middle);

      // Last Runge-Kutta point
      if (energyLossFlag && material) {
        // Update parameters related to a changed momentum
        momentum = initialMomentum + h * dP[2];
        // if (momentum <= momentumCutOff) return false; //Abort propagation
        double E = std::sqrt(momentum * momentum * c2 + massSI * massSI * c4);
        dP[3]    = g * E / (momentum * c2);
        qop[3]   = state.q / momentum;
        // Calculate term for later error propagation
        if (state.covTransport) {
          dL[3] = -qop[3] * qop[3] * g * E
                  * (3. - (momentum * momentum * c2) / (E * E)) / c2
              - qop[3] * qop[3] * qop[3] * E * dgdqopValue;
        }
      }

      const Vector3D pos2 = state.pos + h * state.dir + h2 * 0.5 * k3;
      B_last              = getField(state, pos2);
      k4                  = qop[3] * (state.dir + h * k3).cross(B_last);
std::cout << "h2: " << h2 << std::endl;
      // Return an estimate of the local integration error
      return h2 * (k1 - k2 - k3 + k4).template lpNorm<1>();
    };

    double tolerance = 5e-5;  // TODO: hardcode

    // Select and adjust the appropriate Runge-Kutta step size
    // @todo remove magic numbers and implement better step estimation
    double error_estimate = std::max(tryRungeKuttaStep(state.stepSize), 1e-20);
    std::cout << "call " << error_estimate << std::endl;
    while (error_estimate > 4. * tolerance) {
      state.stepSize = state.stepSize
          * std::min(std::max(0.25,
                              std::pow((tolerance / error_estimate), 0.25)),
                     4.);
      error_estimate = std::max(tryRungeKuttaStep(state.stepSize), 1e-20);
      std::cout << "call" << error_estimate << std::endl;
    }

    // use the adjusted step size
    const double h = state.stepSize;

    // When doing error propagation, update the associated Jacobian matrix
    if (state.covTransport) {
std::cout << "dL: " << dL[0] << "\t" << dL[1] << "\t" << dL[2] << "\t" << dL[3] << std::endl;
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
      const double jdL1 = dL[0];
      dk1dL             = state.dir.cross(B_first);
      const double jdL2 = dL[1] * (1. + half_h * jdL1);
      dk2dL = (1. + half_h * jdL1) * (state.dir + half_h * k1).cross(B_middle)
          + qop[1] * half_h * dk1dL.cross(B_middle);
      const double jdL3 = dL[2] * (1. + half_h * jdL2);
      dk3dL = (1. + half_h * jdL2) * (state.dir + half_h * k2).cross(B_middle)
          + qop[2] * half_h * dk2dL.cross(B_middle);
      const double jdL4 = dL[3] * (1. + h * jdL3);
      dk4dL             = (1. + h * jdL3) * (state.dir + h * k3).cross(B_last)
          + qop[3] * h * dk3dL.cross(B_last);

      dk1dT(0, 1) = B_first.z();
      dk1dT(0, 2) = -B_first.y();
      dk1dT(1, 0) = -B_first.z();
      dk1dT(1, 2) = B_first.x();
      dk1dT(2, 0) = B_first.y();
      dk1dT(2, 1) = -B_first.x();
      dk1dT *= qop[0];

      dk2dT += h / 2 * dk1dT;
      dk2dT *= cross(dk2dT, B_middle);
      dk2dT *= qop[1];

      dk3dT += h / 2 * dk2dT;
      dk3dT *= cross(dk3dT, B_middle);
      dk3dT *= qop[2];

      dk4dT += h * dk3dT;
      dk4dT *= cross(dk4dT, B_last);
      dk4dT *= qop[3];

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      Vector3D onesVec(1., 1., 1.);
      dFdL = onesVec * h + conv * h2 / 6 * (dk1dL + dk2dL + dk3dL);

      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = onesVec + conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      // Evaluation of the dLambda''/dlambda term
      D(6, 6) += conv * (h / 6.) * (jdL1 + 2. * (jdL2 + jdL3) + jdL4);
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

    // Update inverse momentum if energyLossFlag is switched on
    if (energyLossFlag && material) {
std::cout << "p abzug: " << units::SI2Nat<units::MOMENTUM>((h / 6.) * (dP[0] + 2. * (dP[1] + dP[2]) + dP[3])) << std::endl;
      state.p += units::SI2Nat<units::MOMENTUM>((h / 6.) * (dP[0] + 2. * (dP[1] + dP[2]) + dP[3]));
      // if (momentum <= m_momentumCutOff) return false; //Abort propagation
    }
    std::cout << "result pos: " << state.pos << std::endl;
    std::cout << "result dir: " << state.dir << std::endl;
    std::cout << "result p: " << state.p << std::endl;
    std::cout << "result cov:\n" << state.jacTransport << std::endl;
    state.pathAccumulated += h;
    //~ std::exit(1);
    return h;
  }

private:
  /// Magnetic field inside of the detector
  BField m_bField;

const double c2 = units::_c * units::_c;
const double c4 = units::_c * units::_c * units::_c * units::_c;

  detail::IonisationLoss energyLoss;

  ActsMatrixD<3, 3>
  cross(const ActsMatrixD<3, 3>& m, const Vector3D& v) const
  {
    ActsMatrixD<3, 3> r;
    r.col(0) = m.col(0).cross(v);
    r.col(1) = m.col(1).cross(v);
    r.col(2) = m.col(2).cross(v);

    return r;
  }
};

}  // namespace Acts
