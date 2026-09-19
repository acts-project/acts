// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <string>
#include <tuple>

namespace Acts {

class IVolumeMaterial;

/// @brief straight line stepper based on Surface intersection
///
/// The straight line stepper is a simple navigation stepper
/// to be used to navigate through the tracking geometry. It can be
/// used for simple material mapping, navigation validation
class StraightLineStepper final {
 public:
  /// Type alias for bound track parameters
  using BoundParameters = BoundTrackParameters;
  /// Type alias for transport jacobian matrix
  using Jacobian = BoundMatrix;
  /// Type alias for covariance matrix
  using Covariance = BoundMatrix;
  /// Type alias for magnetic field (null field for straight line propagation)
  using BField = NullBField;

  struct Config {};

  /// Configuration options for straight line propagation.
  struct Options : public StepperPlainOptions {
    /// Constructor from geometry and magnetic field contexts
    /// @param gctx The geometry context
    /// @param mctx The magnetic field context
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    /// Set plain stepper options
    /// @param options The plain options to set
    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  /// State for track parameter propagation
  ///
  struct State {
    /// Constructor from the initial bound track parameters
    ///
    /// @param [in] optionsIn The options for the stepper
    ///
    /// @note the covariance matrix is copied when needed
    explicit State(const Options& optionsIn) : options(optionsIn) {}

    /// Configuration options for the stepper
    Options options;

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from runge kutta integration
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// Boolean to indicate if you need covariance transport
    bool covTransport = false;
    /// Covariance matrix for track parameter uncertainties
    Covariance cov = Covariance::Zero();

    /// accumulated path length state
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    /// Previous step size for overstep estimation (ignored for straight line
    /// stepper)
    double previousStepSize = 0.;

    /// Statistics of the stepper
    StepperStatistics statistics;
  };

  /// Create a stepper state from propagation options
  /// @param options The propagation options
  /// @return A new stepper state initialized with the provided options
  State makeState(const Options& options) const;

  /// Initialize the stepper state from bound track parameters
  /// @param state The stepper state to initialize
  /// @param par The bound track parameters to initialize from
  void initialize(State& state, const BoundParameters& par) const;

  /// Initialize the stepper state from bound parameters and components
  /// @param state The stepper state to initialize
  /// @param boundParams The bound parameter vector
  /// @param cov Optional covariance matrix
  /// @param particleHypothesis The particle hypothesis (mass, charge, etc.)
  /// @param surface The reference surface
  void initialize(State& state, const BoundVector& boundParams,
                  const std::optional<BoundMatrix>& cov,
                  ParticleHypothesis particleHypothesis,
                  const Surface& surface) const;

  /// Get the field for the stepping, this gives back a zero field
  /// @return Always returns zero magnetic field vector for straight-line propagation
  Result<Vector3> getField(State& /*state*/, const Vector3& /*pos*/) const {
    // get the field from the cell
    return Result<Vector3>::success({0., 0., 0.});
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Current global position vector
  Vector3 position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Current normalized direction vector
  Vector3 direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// QoP direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Charge over momentum (q/p) value
  double qOverP(const State& state) const { return state.pars[eFreeQOverP]; }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Absolute momentum magnitude
  double absoluteMomentum(const State& state) const {
    return particleHypothesis(state).extractMomentum(qOverP(state));
  }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Current momentum vector
  Vector3 momentum(const State& state) const {
    return absoluteMomentum(state) * direction(state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Electric charge of the particle
  double charge(const State& state) const {
    return particleHypothesis(state).extractCharge(qOverP(state));
  }

  /// Particle hypothesis
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Reference to the particle hypothesis used
  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return state.particleHypothesis;
  }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The time coordinate from the free parameters vector
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Update surface status
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] boundaryTolerance The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] stype The step size type to be set
  /// @param [in] logger A logger instance
  /// @return Status of the intersection indicating whether surface was reached
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        *this, state, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
  }

  /// Update step size
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param target [in] The NavigationTarget
  /// @param direction [in] The propagation direction
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, const NavigationTarget& target,
                      Direction direction, ConstrainedStep::Type stype) const {
    static_cast<void>(direction);
    double stepSize = target.pathLength();
    updateStepSize(state, stepSize, stype);
  }

  /// Update step size - explicitly with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    state.previousStepSize = state.stepSize.value();
    state.stepSize.update(stepSize, stype);
  }

  /// Get the step size
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @param stype [in] The step size type to be returned
  /// @return Current step size for the specified constraint type
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Release the Step size
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    state.stepSize.release(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @return String representation of the current step size
  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  /// Get the step size constraints
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The step size constraints
  const ConstrainedStep& stepSize(const State& state) const {
    return state.stepSize;
  }

  /// Get the stepper statistics
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The statistics since the last initialization
  const StepperStatistics& statistics(const State& state) const {
    return state.statistics;
  }

  /// Get the path length
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The path length since the last initialization
  double pathLength(const State& state) const { return state.pathAccumulated; }

  /// Check if the state carries a covariance
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return True if the covariance is transported
  bool hasCovariance(const State& state) const { return state.covTransport; }

  /// Get the covariance at the anchor
  ///
  /// The anchor is the frame of the last initialization, transport or update.
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The covariance at the anchor
  const Covariance& covariance(const State& state) const { return state.cov; }

  /// Set the covariance at the anchor
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] covariance The new covariance at the anchor
  void setCovariance(State& state, const Covariance& covariance) const {
    state.cov = covariance;
  }

  /// Get the bound parameters at the current position
  ///
  /// The parameters carry the covariance at the anchor if the state has one.
  ///
  /// @note It does not check if the state is on @p surface or anchored on it
  ///
  /// @param [in] state The stepping state (thread-local cache)
  /// @param [in] surface The surface of the parameters
  /// @return The bound parameters, or a failure if the position cannot be
  ///         expressed on @p surface
  Result<BoundParameters> boundParameters(const State& state,
                                          const Surface& surface) const;

  /// @brief If necessary fill additional members needed for
  /// transportToCurvilinear
  ///
  /// Compute path length derivatives in case they have not been computed
  /// yet, which is the case if no step has been executed yet.
  ///
  /// @param [in, out] state The stepping state (thread-local cache)
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const {
    // test whether the accumulated path has still its initial value.
    if (state.pathAccumulated != 0) {
      return true;
    }

    // dr/ds :
    state.derivative.template head<3>() = direction(state);
    // dt / ds
    state.derivative(eFreeTime) = fastHypot(
        1., state.particleHypothesis.mass() / absoluteMomentum(state));
    // d (dr/ds) / ds : == 0
    state.derivative.template segment<3>(4) = Acts::Vector3::Zero().transpose();
    // d qop / ds  == 0
    state.derivative(eFreeQOverP) = 0.;

    return true;
  }

  /// Get the curvilinear parameters at the current position
  ///
  /// The parameters carry the covariance at the anchor if the state has one.
  ///
  /// @param [in] state The stepping state (thread-local cache)
  /// @return The curvilinear parameters
  BoundParameters curvilinearParameters(const State& state) const;

  /// Method to update a stepper state to the some parameters
  ///
  /// This anchors the state on @p surface.
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] freeParams Free parameters that will be written into @p state
  /// @param [in] boundParams Corresponding bound parameters used to update jacToGlobal in @p state
  /// @param [in] covariance Covariance that will be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update the stepper state
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] qop the updated qop value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qop, double time) const;

  /// Transport the covariance to the curvilinear frame at the current position
  ///
  /// This anchors the state on the curvilinear frame. Without a covariance
  /// the state does not change.
  ///
  /// @param [in,out] state State of the stepper
  /// @return The jacobian from the previous anchor to the curvilinear frame
  Jacobian transportToCurvilinear(State& state) const;

  /// Transport the covariance to a surface at the current position
  ///
  /// This anchors the state on @p surface. Without a covariance the state
  /// does not change.
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface The surface to transport the covariance to
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  /// @return The jacobian from the previous anchor to @p surface, or a failure
  ///         if the parameters cannot be expressed on @p surface
  Result<Jacobian> transportToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// Perform a straight line propagation step
  ///
  /// @param [in,out] state State of the stepper
  /// @param propDir is the direction of propagation
  /// @param material is the optional volume material we are stepping through.
  //         This is simply ignored if `nullptr`.
  /// @return the result of the step
  ///
  /// @note The state contains the desired step size. It can be negative during
  ///       backwards track propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const {
    static_cast<void>(material);

    // use the adjusted step size
    const auto h = state.stepSize.value() * propDir;
    const auto m = state.particleHypothesis.mass();
    const auto p = absoluteMomentum(state);
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = fastHypot(1, m / p);
    // Update the track parameters according to the equations of motion
    Vector3 dir = direction(state);
    state.pars.template segment<3>(eFreePos0) += h * dir;
    state.pars[eFreeTime] += h * dtds;

    // Propagate the jacobian
    if (state.covTransport) {
      // The step transport matrix in global coordinates
      FreeMatrix D = FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = SquareMatrix<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // d(t)/d(q/p) = h m^2 (q/p) / (q^2 dt/ds), with the q^2 folded into p
      // via p = |q| / |q/p|.
      D(3, 7) = h * m * m / (p * p * state.pars[eFreeQOverP] * dtds);
      // Set the derivative factor the time
      state.derivative(3) = dtds;
      // Update jacobian and derivative
      state.jacTransport = D * state.jacTransport;
      state.derivative.template head<3>() = dir;
    }

    // state the path length
    state.pathAccumulated += h;
    ++state.nSteps;

    ++state.statistics.nAttemptedSteps;
    ++state.statistics.nSuccessfulSteps;
    if (propDir != Direction::fromScalarZeroAsPositive(h)) {
      ++state.statistics.nReverseSteps;
    }
    state.statistics.pathLength += h;
    state.statistics.absolutePathLength += std::abs(h);

    return h;
  }
};

template <>
struct SupportsBoundParameters<StraightLineStepper> : public std::true_type {};

}  // namespace Acts
