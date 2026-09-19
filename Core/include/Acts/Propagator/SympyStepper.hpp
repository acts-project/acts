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
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/MaterialEffectsAccumulator.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"

namespace Acts {

class IVolumeMaterial;

/// Stepper implementation using sympy-generated expressions.
class SympyStepper final {
 public:
  /// Type alias for bound track parameters
  using BoundParameters = BoundTrackParameters;
  /// Type alias for jacobian matrix
  using Jacobian = BoundMatrix;
  /// Type alias for covariance matrix
  using Covariance = BoundMatrix;

  /// Configuration for the sympy stepper.
  struct Config {
    /// Magnetic field provider
    std::shared_ptr<const MagneticFieldProvider> bField;
  };

  /// Runtime options for sympy propagation.
  struct Options : public StepperPlainOptions {
    /// Whether to perform dense output
    bool doDense = true;
    /// Maximum radiation length fraction per step
    double maxXOverX0Step = 1;

    /// Constructor
    /// @param gctx Geometry context
    /// @param mctx Magnetic field context
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    /// Set plain stepper options
    /// @param options Plain stepper options
    void setPlainOptions(const StepperPlainOptions& options) {
      static_cast<StepperPlainOptions&>(*this) = options;
    }
  };

  /// @brief State for track parameter propagation
  ///
  /// It contains the stepping information and is provided thread local
  /// by the propagator
  struct State {
    /// Constructor from the initial bound track parameters
    ///
    /// @param [in] optionsIn is the configuration of the stepper
    /// @param [in] fieldCacheIn is the cache object for the magnetic field
    ///
    /// @note the covariance matrix is copied when needed
    State(const Options& optionsIn, MagneticFieldProvider::Cache fieldCacheIn)
        : options(optionsIn), fieldCache(std::move(fieldCacheIn)) {}

    // Declaration order matters: members used by `step()` are kept in one
    // contiguous run of cache lines, the rest come last.

    /// Configuration options for the stepper
    Options options;

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Bound-to-free jacobian from the last reference surface, transported
    /// along with the track.
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Covariance matrix (and indicator)
    /// associated with the initial error on track parameters
    bool covTransport = false;

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// dt/ds, handed to the vacuum kernel rather than formed in it. Constant
    /// while q/p is, so it is refreshed wherever q/p moves: initialize(),
    /// update() and a dense step.
    double dtds = 1;

    /// Adaptive step size of the runge-kutta integration
    ConstrainedStep stepSize;

    /// Last performed step (for overstep limit calculation)
    double previousStepSize = 0.;

    /// Magnetic field at the current position, reused as the next step's
    /// first sample. Reset whenever the parameters are set from outside.
    std::optional<Vector3> field;

    /// Accumulated path length state
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// Statistics of the stepper
    StepperStatistics statistics;

    /// Accumulator for material effects along the trajectory
    detail::MaterialEffectsAccumulator materialEffectsAccumulator;

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    MagneticFieldProvider::Cache fieldCache;

    // Only used when transporting to a surface:

    /// Covariance matrix for error propagation
    Covariance cov = Covariance::Zero();
  };

  /// Constructor requires knowledge of the detector's magnetic field
  /// @param bField The magnetic field provider
  explicit SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField);

  /// @brief Constructor with configuration
  /// @param config The configuration of the stepper
  explicit SympyStepper(const Config& config);

  /// Create a state object
  /// @param options Stepper options
  /// @return State object
  State makeState(const Options& options) const;

  /// Initialize the state from bound track parameters
  /// @param state The state to initialize
  /// @param par The bound track parameters
  void initialize(State& state, const BoundParameters& par) const;

  /// Initialize the state from bound parameters
  /// @param state The state to initialize
  /// @param boundParams Bound track parameters vector
  /// @param cov Covariance matrix
  /// @param particleHypothesis Particle hypothesis
  /// @param surface Reference surface
  void initialize(State& state, const BoundVector& boundParams,
                  const std::optional<BoundMatrix>& cov,
                  ParticleHypothesis particleHypothesis,
                  const Surface& surface) const;

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  /// @return Magnetic field vector
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    // get the field from the cell
    return m_bField->getField(pos, state.fieldCache);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Position vector
  Vector3 position(const State& state) const {
    return state.pars.template segment<3>(eFreePos0);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Direction vector
  Vector3 direction(const State& state) const {
    return state.pars.template segment<3>(eFreeDir0);
  }

  /// QoP direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Charge over momentum
  double qOverP(const State& state) const { return state.pars[eFreeQOverP]; }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Absolute momentum
  double absoluteMomentum(const State& state) const {
    return particleHypothesis(state).extractMomentum(qOverP(state));
  }

  /// Momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Momentum vector
  Vector3 momentum(const State& state) const {
    return absoluteMomentum(state) * direction(state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Particle charge
  double charge(const State& state) const {
    return particleHypothesis(state).extractCharge(qOverP(state));
  }

  /// Particle hypothesis
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Particle hypothesis
  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return state.particleHypothesis;
  }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return Time
  double time(const State& state) const { return state.pars[eFreeTime]; }

  /// Update surface status
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] surface The surface provided
  /// @param [in] index The surface intersection index
  /// @param [in] navDir The navigation direction
  /// @param [in] boundaryTolerance The boundary check for this status update
  /// @param [in] surfaceTolerance Surface tolerance used for intersection
  /// @param [in] stype The step size type to be set
  /// @param [in] logger A @c Logger instance
  /// @return Intersection status
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return detail::updateSingleSurfaceStatus<SympyStepper>(
        *this, state, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
  }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
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
  /// @return Step size
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param [in] stype The step size type to be released
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    state.stepSize.release(stype);
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @return String representation of step size
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
  /// @param [in, out] state The state of the stepper
  /// @return true if nothing is missing after this call, false otherwise.
  bool prepareCurvilinearState(State& state) const;

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
  /// @param [in] covariance The covariance that will be written into @p state
  /// @param [in] surface The surface used to update the jacToGlobal
  void update(State& state, const FreeVector& freeParams,
              const BoundVector& boundParams, const Covariance& covariance,
              const Surface& surface) const;

  /// Method to update the stepper state
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] qOverP the updated qOverP value
  /// @param [in] time the updated time value
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double qOverP, double time) const;

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
  /// does not change. The state must be on @p surface, and it stays unchanged
  /// if it is not.
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface The surface to transport the covariance to
  /// @param [in] freeToBoundCorrection Correction for non-linearity effect during transform from free to bound
  /// @return The jacobian from the previous anchor to @p surface, or a failure
  ///         if the state is not on @p surface
  Result<Jacobian> transportToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const;

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state State of the stepper
  /// @param propDir is the direction of propagation
  /// @param material is the optional volume material we are stepping through.
  //         This is simply ignored if `nullptr`.
  /// @return the result of the step
  ///
  /// @note The state contains the desired step size. It can be negative during
  ///       backwards track propagation, and since we're using an adaptive
  ///       algorithm, it can be modified by the stepper class during
  ///       propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const;

 protected:
  /// Magnetic field inside of the detector
  std::shared_ptr<const MagneticFieldProvider> m_bField;
};

template <>
struct SupportsBoundParameters<SympyStepper> : public std::true_type {};

}  // namespace Acts
