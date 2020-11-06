// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/ConstrainedStepControl.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/StepperExtensionList.hpp"
#include "Acts/Propagator/StepperState.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Units.hpp"

#include <cmath>
#include <functional>
#include <limits>

namespace Acts {

using namespace Acts::UnitLiterals;

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
/// p the momentum magnitude and B the magnetic field
///
/// Depending on the extions chosen from the extension list, either a pure
/// Runge-Kutta numerical integration of the equation of motion (see above),
/// or transport through material dense volumes is done.
///
/// A nested State object guarantees for the thread-local caching of the
/// transport parameters. This state object is provided by the Propagator
/// and the stepper is used to interpret the non-trivial internal parameters,
/// e.g. the current global position is gathered through:
///
///   auto position = stepper.position(state);
///
template <typename bfield_t,
          typename extensionlist_t = StepperExtensionList<DefaultExtension>,
          typename auctioneer_t = detail::VoidAuctioneer>
class EigenStepper {
 public:
  /// Jacobian, Covariance and State defintions
  using Jacobian = BoundMatrix;
  using Covariance = BoundSymMatrix;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;
  using BField = bfield_t;

  using StateBase =
      StepperState<EigenStepper<bfield_t, extensionlist_t, auctioneer_t> >;

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator.
  ///
  /// The Stepper and the extensions are allowed to directly manipulate
  /// the state for perormance reasons, while all actors & aborters
  /// have to go through the public method interface
  struct State final : public StateBase {
    /// Access control to state: Stepper
    friend EigenStepper<bfield_t, extensionlist_t, auctioneer_t>;
    /// Access control to state: Default extension
    friend DefaultExtension;
    /// Access control to state: Dense extension
    friend DenseEnvironmentExtension;

    /// Constructor from the initial track parameters
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direciton w.r.t momentum
    /// @param [in] ssize is the maximum step size
    /// @param [in] stolerance is the stepping tolerance
    ///
    /// @note the covariance matrix is copied when needed
    template <typename parameters_t>
    explicit State(std::reference_wrapper<const GeometryContext> gctx,
                   std::reference_wrapper<const MagneticFieldContext> mctx,
                   const parameters_t& par, NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max(),
                   double stolerance = s_onSurfaceTolerance)
        : StateBase(gctx, mctx, par, ndir, ssize, stolerance),
          fieldCache(mctx) {}

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    typename BField::Cache fieldCache;

    /// List of algorithmic extensions
    extensionlist_t extension;

    /// Auctioneer for choosing the extension
    auctioneer_t auctioneer;

    /// @brief Storage of magnetic field and the sub steps during a RKN4 step
    struct {
      /// Magnetic field evaulations
      Vector3D B_first, B_middle, B_last;
      /// k_i of the RKN4 algorithm
      Vector3D k1, k2, k3, k4;
      /// k_i elements of the momenta
      std::array<double, 4> kQoP;
    } stepData;
  };

  using StepControl = ConstrainedStepControl<
      EigenStepper<bfield_t, extensionlist_t, auctioneer_t> >;
  StepControl stepControl;

  /// Constructor requires knowledge of the detector's magnetic field
  ///
  /// @param bField The magnetic field provided to the state
  EigenStepper(BField bField);

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reference surface of the bound parameters
  /// @param [in] navDir Navigation direction
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams, const BoundSymMatrix& cov,
      const Surface& surface, const NavigationDirection navDir = forward,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D getField(State& state, const Vector3D& pos) const {
    // get the field from the cell
    return m_bField.getField(pos, state.fieldCache);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double momentum(const State& state) const {
    return std::abs(1. / state.pars[eFreeQOverP]);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return state.q; }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Access to the current geometry context
  ///
  /// @param state [in] The stepping state (thread-local cache)
  const GeometryContext& geometryContext(const State& state) const {
    return state.geoContext;
  }

  /// Access to the navigation direction
  ///
  /// @param state [in] The stepping state (thread-local cache)
  NavigationDirection steppingDirection(const State& state) const {
    return state.navDir;
  }

  /// Access to the navigation direction
  ///
  /// @param state [in, out] The stepping state (thread-local cache)
  /// @param sdir [in] stepping direction
  void setSteppingDirection(State& state, NavigationDirection sdir) const {
    state.navDir = sdir;
  }

  /// Access to the stepping tolerance
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double steppingTolerance(const State& state) const { return state.tolerance; }

  /// Access to the accumulated path
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double accumulatedPath(const State& state) const {
    return state.pathAccumulated;
  }

  /// Reset to the accumulated path
  ///
  /// @param state [in, out] The stepping state (thread-local cache)
  void resetAccumulatedPath(State& state) const {
    state.pathAccumulated = 0;
    ;
  }

  /// Indicate if the covariance has to be transported
  ///
  /// @param state [in] The stepping state (thread-local cache)
  bool transportCovariance(const State& state) const {
    return state.covTransport;
  }

  /// Update surface status
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param surface [in] The surface provided
  /// @param bcheck [in] The boundary check for this status update
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck) const {
    return detail::updateSingleSurfaceStatus<EigenStepper>(*this, state,
                                                           surface, bcheck);
  }

  /// Overstep limit
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double overstepLimit(const State& /*state*/) const {
    // A dynamic overstep limit could sit here
    return -m_overstepLimit;
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  BoundState boundState(State& state, const Surface& surface,
                        bool transportCov = true) const;

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const;

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  void update(State& state, const FreeVector& parameters,
              const Covariance& covariance) const;

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  void update(State& state, const Vector3D& uposition,
              const Vector3D& udirection, double up, double time) const;

  /// @brief Convenience method for better readability
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] delta The change that may be applied to it
  ///
  void updateBoundVariance(State& state, BoundIndices bIndex,
                           double delta) const;

  /// @brief Convenience method for better readability
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] cov the bound covariance matrix to be updated
  ///
  void updateBoundCovariance(State& state, Covariance cov) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  ///
  /// @return the full transport jacobian
  void covarianceTransport(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @tparam surface_t the Surface type
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded to
  /// @note no check is done if the position is actually on the surface
  void covarianceTransport(State& state, const Surface& surface) const;

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const;

 private:
  /// Magnetic field inside of the detector
  BField m_bField;

  /// Overstep limit: could/should be dynamic
  double m_overstepLimit = 100_um;
};
}  // namespace Acts

#include "Acts/Propagator/EigenStepper.ipp"
