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
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Propagator/detail/StepperReturnState.hpp"

#include <cmath>
#include <functional>

namespace Acts {

/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
class StraightLineStepper {
 public:
  using Jacobian = BoundMatrix;
  using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;
  using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<CurvilinearTrackParameters, Jacobian, double>;
  using FreeState = std::tuple<FreeTrackParameters, std::variant<FreeMatrix, BoundToFreeMatrix>, double>;
  using BField = NullBField;

  /// State for track parameter propagation
  ///
  struct State {
    /// Delete the default constructor
    State() = delete;

    /// Constructor from the initial track parameters
    ///
    /// @tparam parameters_t the Type of the track parameters
    ///
    /// @param [in] gctx is the context object for the geometery
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir is the navigation direction
    /// @param [in] ssize is the (absolute) maximum step size
    /// @param [in] stolerance is the stepping tolerance
    template <typename parameters_t>
    explicit State(std::reference_wrapper<const GeometryContext> gctx,
                   std::reference_wrapper<const MagneticFieldContext> /*mctx*/,
                   const parameters_t& par, NavigationDirection ndir = forward,
                   double ssize = std::numeric_limits<double>::max(),
                   double stolerance = s_onSurfaceTolerance)
        : pos(par.position(gctx)),
          dir(par.unitDirection()),
          p(par.absoluteMomentum()),
          q(par.charge()),
          t(par.time()),
          navDir(ndir),
          stepSize(ndir * std::abs(ssize)),
          tolerance(stolerance),
          geoContext(gctx) {
      if (par.covariance()) {
        // Get the reference surface for navigation
        const auto& surface = par.referenceSurface();
        // set the covariance transport flag to true and copy
        covTransport = true;
        cov = BoundSymMatrix(*par.covariance());
        surface.initJacobianToGlobal(gctx, jacToGlobal, pos, dir,
                                     par.parameters());
      }
    }

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from runge kutta integration
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The full jacobian of the transport entire transport
    Jacobian jacobian = Jacobian::Identity();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Boolean to indiciate if you need covariance transport
    bool covTransport = false;
    Covariance cov = Covariance::Zero();

    /// Global particle position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1., 0., 0.);

    /// Momentum
    double p = 0.;

    /// Save the charge: neutral as default for SL stepper
    double q = 0.;

    /// Propagated time
    double t = 0.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize = std::numeric_limits<double>::max();

    // Previous step size for overstep estimation (ignored for SL stepper)
    double previousStepSize = 0.;

    /// The tolerance for the stepping
    double tolerance = s_onSurfaceTolerance;

    // Cache the geometry context of this propagation
    std::reference_wrapper<const GeometryContext> geoContext;
  };

  /// Always use the same propagation state type, independently of the initial
  /// track parameter type and of the target surface
  using state_type = State;

  /// Constructor
  StraightLineStepper() = default;

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] freeParams Parameters in free parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] navDir Navigation direction
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams, const BoundSymMatrix& cov,
      const Surface& surface, const NavigationDirection navDir = forward,
      const double stepSize = std::numeric_limits<double>::max()) const;

  /// Get the field for the stepping, this gives back a zero field
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D getField(State& /*state*/, const Vector3D& /*pos*/) const {
    // get the field from the cell
    return Vector3D(0., 0., 0.);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D position(const State& state) const { return state.pos; }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3D direction(const State& state) const { return state.dir; }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double momentum(const State& state) const { return state.p; }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return state.q; }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return state.t; }

  /// Overstep limit
  ///
  /// @param state The stepping state (thread-local cache)
  double overstepLimit(const State& /*state*/) const {
    return s_onSurfaceTolerance;
  }

  /// Update surface status
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param surface [in] The surface provided
  /// @param bcheck [in] The boundary check for this status update
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck) const {
    return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        *this, state, surface, bcheck);
  }

  /// Update step size
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
    detail::updateSingleStepSize<StraightLineStepper>(state, oIntersection,
                                                      release);
  }

  /// Set Step size - explicitely with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void setStepSize(State& state, double stepSize,
                   ConstrainedStep::Type stype = ConstrainedStep::actor) const {
    state.previousStepSize = state.stepSize;
    state.stepSize.update(stepSize, stype, true);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  void releaseStepSize(State& state) const {
    state.stepSize.release(ConstrainedStep::actor);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  template<typename start_parameters_t, typename end_parameters_t = start_parameters_t>
  auto 
  buildState(State& state, bool reinitialize) const
  {	  
	  using return_type = detail::return_state_type<start_parameters_t, end_parameters_t>;
	  if constexpr (end_parameters_t::is_local_representation)
	  {
		 return curvilinearState<return_type>(state, reinitialize);
	  }
	  else
	  {
		   return freeState<return_type>(state, reinitialize);
	  }
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief It does not check if the transported state is at the surface, this
  /// needs to be guaranteed by the propagator
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
  /// @brief This creates a curvilinear state.
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
  /// @param [in] time the updated time value
  void update(State& state, const Vector3D& uposition,
              const Vector3D& udirection, double up, double time) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @param [in,out] state State of the stepper
  void covarianceTransport(State& state) const;

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state - for the moment a dummy method
  ///
  /// @tparam surface_t the surface type - ignored here
  ///
  /// @param [in,out] state The stepper state
  /// @param [in] surface is the surface to which the covariance is
  ///        forwarded to
  /// @note no check is done if the position is actually on the surface
  ///
  void covarianceTransport(State& state, const Surface& surface) const;

  /// Perform a straight line propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///                The state contains the desired step size,
  ///                it can be negative during backwards track propagation,
  ///                and since we're using an adaptive algorithm, it can
  ///                be modified by the stepper class during propagation.
  ///
  /// @return the step size taken
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const {
    // use the adjusted step size
    const auto h = state.stepping.stepSize;
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = std::hypot(1., state.options.mass / state.stepping.p);
    // Update the track parameters according to the equations of motion
    state.stepping.pos += h * state.stepping.dir;
    state.stepping.t += h * dtds;
    // Propagate the jacobian
    if (state.stepping.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = ActsSymMatrixD<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * state.options.mass * state.options.mass *
                (state.stepping.q == 0. ? 1. : state.stepping.q) /
                (state.stepping.p * dtds);
      // Set the derivative factor the time
      state.stepping.derivative(3) = dtds;
      // Update jacobian and derivative
      state.stepping.jacTransport = D * state.stepping.jacTransport;
      state.stepping.derivative.template head<3>() = state.stepping.dir;
    }
    // state the path length
    state.stepping.pathAccumulated += h;

    // return h
    return h;
  }
  
private:

	/// @brief This function treats the modifications of the jacobian related to the projection onto a surface. Since a variation of the start parameters within a given uncertainty would lead to a variation of the end parameters, these need to be propagated onto the target surface. This an approximated approach to treat the (assumed) small change.
	///
	/// @param [in] state The current state
	/// @param [in, out] jacToFree The jacobian from the local start parameters to the propagated global end parameters
	/// @param [in] surface The surface onto which the projection should be performed
	/// @note The parameter @p surface is only required if projected to bound parameters. In the case of curvilinear parameters the geometry and the position is known and the calculation can be simplified
	///
	/// @return The projection jacobian from global end parameters to its local equivalent
	const FreeToBoundMatrix surfaceDerivative(State& state, const Surface* surface = nullptr) const
	{
		if(state.jacToGlobal.has_value())
		{
			// Set the surface projection contributions
			// If no surface is specified it is curvilinear
			if(surface == nullptr)
			{
				// Transport the covariance
				const ActsRowVectorD<3> normVec(state.dir);
				const BoundRowVector sfactors =
					normVec * (*state.jacToGlobal).template topLeftCorner<3, BoundParsDim>();
				*state.jacToGlobal -= state.derivative * sfactors;
				// Since the jacobian to local needs to calculated for the bound parameters here, it is convenient to do the same here
				return freeToCurvilinearJacobian(state);
			}
			// Else it is bound
			else
			{
				// Initialize the transport final frame jacobian
				FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
				// Initalize the jacobian to local, returns the transposed ref frame
			   auto rframeT = surface->initJacobianToLocal(state.geoContext, jacToLocal,
							   state.pos, state.dir);
				// Calculate the form factors for the derivatives
				const BoundRowVector sVec = surface->derivativeFactors(
					state.geoContext, state.pos, state.dir, rframeT, (*state.jacToGlobal));
				*state.jacToGlobal -= state.derivative * sVec;
				// Return the jacobian to local
				return jacToLocal;
			}
		}
		else
		{
			// Set the surface projection contributions
			// If no surface is specified it is curvilinear
			if(surface == nullptr)
			{
				// Transport the covariance
				const ActsRowVectorD<3> normVec(state.dir);
				const FreeRowVector sfactors =
					normVec * state.jacTransport.template topLeftCorner<3, FreeParsDim>();
				// Since the jacobian to local needs to calculated for the bound parameters here, it is convenient to do the same here
				return freeToCurvilinearJacobian(state) * (state.jacTransport - state.derivative * sfactors);
			}
			// Else it is bound
			else
			{
				// Initialize the transport final frame jacobian
				FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
				// Initalize the jacobian to local, returns the transposed ref frame
			   auto rframeT = surface->initJacobianToLocal(state.geoContext, jacToLocal,
							   state.pos, state.dir);
				// Calculate the form factors for the derivatives
				const FreeRowVector sVec = surface->derivativeFactors(
					state.geoContext, state.pos, state.dir, rframeT, state.jacTransport);
				// Return the jacobian to local
				return jacToLocal * (state.jacTransport - state.derivative * sVec);
			}
		}
	}

	/// @brief This function reinitialises the @p state member @p jacToGlobal.
	///
	/// @param [in, out] state The state object
	/// @param [in] surface The surface the represents the local parametrisation
	/// @note The surface is only required for bound parameters since it serves to derive the jacobian from it. In the case of curvilinear parameters this is not needed and can be evaluated without any surface.
     void reinitializeJacToGlobal(State& state, const Surface* surface = nullptr) const
     {
		using VectorHelpers::phi;
		using VectorHelpers::theta;
		
		 // Reset the jacobian
		state.jacToGlobal = BoundToFreeMatrix::Zero();
		
		// Fill the jacobian to global for next transport
		// If treating curvilinear parameters
		if(surface == nullptr)
		{
			auto& jac = *state.jacToGlobal;
			// TODO: This was calculated before - can it be reused?
			// Optimized trigonometry on the propagation direction
			const double x = state.dir(0);  // == cos(phi) * sin(theta)
			const double y = state.dir(1);  // == sin(phi) * sin(theta)
			const double z = state.dir(2);  // == cos(theta)
			// can be turned into cosine/sine
			const double cosTheta = z;
			const double sinTheta = sqrt(x * x + y * y);
			const double invSinTheta = 1. / sinTheta;
			const double cosPhi = x * invSinTheta;
			const double sinPhi = y * invSinTheta;

		  jac(0, eLOC_0) = -sinPhi;
		  jac(0, eLOC_1) = -cosPhi * cosTheta;
		  jac(1, eLOC_0) = cosPhi;
		  jac(1, eLOC_1) = -sinPhi * cosTheta;
		  jac(2, eLOC_1) = sinTheta;
		  jac(3, eT) = 1;
		  jac(4, ePHI) = -sinTheta * sinPhi;
		  jac(4, eTHETA) = cosTheta * cosPhi;
		  jac(5, ePHI) = sinTheta * cosPhi;
		  jac(5, eTHETA) = cosTheta * sinPhi;
		  jac(6, eTHETA) = -sinTheta;
		  jac(7, eQOP) = 1;
		}
		// If treating bound parameters
		else
		{
		  Vector2D loc{0., 0.};
		  surface->globalToLocal(state.geoContext, state.pos, state.dir, loc);
		  BoundVector pars;
		  pars << loc[eLOC_0], loc[eLOC_1], phi(state.dir), theta(state.dir),
			  state.q / state.p, state.t0 + state.dt;
		  surface->initJacobianToGlobal(state.geoContext, *state.jacToGlobal,
									   state.pos, state.dir, pars);
		}
	 }
	 
	/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
	///
	/// @param [in] state State that will be projected
	///
	/// @return Projection Jacobian
	FreeToBoundMatrix freeToCurvilinearJacobian(const State& state) const
	{
		// Optimized trigonometry on the propagation direction
		const double x = state.dir(0);  // == cos(phi) * sin(theta)
		const double y = state.dir(1);  // == sin(phi) * sin(theta)
		const double z = state.dir(2);  // == cos(theta)
		// can be turned into cosine/sine
		const double cosTheta = z;
		const double sinTheta = sqrt(x * x + y * y);
		const double invSinTheta = 1. / sinTheta;
		const double cosPhi = x * invSinTheta;
		const double sinPhi = y * invSinTheta;
		// prepare the jacobian to curvilinear
		FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
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
		  const double c = sqrt(y * y + z * z);
		  const double invC = 1. / c;
		  jacToCurv(0, 1) = -z * invC;
		  jacToCurv(0, 2) = y * invC;
		  jacToCurv(1, 0) = c;
		  jacToCurv(1, 1) = -x * y * invC;
		  jacToCurv(1, 2) = -x * z * invC;
		}
		// Time parameter
		jacToCurv(5, 3) = 1.;
		// Directional and momentum parameters for curvilinear
		jacToCurv(2, 4) = -sinPhi * invSinTheta;
		jacToCurv(2, 5) = cosPhi * invSinTheta;
		jacToCurv(3, 6) = -invSinTheta;
		jacToCurv(4, 7) = 1.;
		
		return jacToCurv;
	}
};

}  // namespace Acts
