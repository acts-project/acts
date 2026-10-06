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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/ButcherTableau.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorTraits.hpp"
#include "Acts/Propagator/StepperOptions.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"

#include <optional>
#include <stdexcept>
#include <string>

namespace Acts {

class IVolumeMaterial;

/// @brief Reference stepper with a configurable Runge-Kutta tableau
///
/// The stepper integrates the equation of motion in vacuum for the free
/// parameters y = (r, t, T, q/p),
///
///   dy/ds = f(y) = (T, 1/beta, (q/p) T x B(r), 0),
///
/// with any explicit @ref ButcherTableau. It is not optimised. It is meant as a
/// reference to compare other steppers against.
///
/// The transport jacobian is the exact derivative of the discrete step. The
/// stepper differentiates each stage with the chain rule,
///
///   dY_i = I + h sum_j a_ij dK_j,  dK_i = J_f(Y_i) dY_i,
///   D = I + h sum_i b_i dK_i,
///
/// so it works for every tableau. Here dY_i, dK_i and D are the derivatives
/// of the stage, of the stage slope and of the step by the free parameters at
/// the start of the step. J_f = df/dy has the non-zero blocks
///
///   d(r)/dT = I,  d(t)/d(q/p) = d(1/beta)/d(q/p),
///   d(T)/dr = (q/p) [T x] dB/dr,  d(T)/dT = -(q/p) [B x],
///   d(T)/d(q/p) = T x B.
///
/// The stepper computes dB/dr with central finite differences of the field,
/// so the jacobian is exact up to the O(epsilon^2) error of the gradient, see
/// @ref Options::fieldGradientEpsilon. Without the gradient term
/// (@ref Options::includeFieldGradient) the jacobian is not exact.
///
/// After the step, the stepper normalises the direction to unit length. This
/// is part of the discrete step, so the stepper also applies its derivative
/// (I - T T^T / |T|^2) / |T| to the direction rows of D.
///
/// With embedded weights the stepper adapts the step size to
/// @ref StepperPlainOptions::stepTolerance. The error estimate is the maximum
/// norm of the difference between the solution and the embedded solution,
/// before the direction is normalised. It uses the position (mm), the time
/// (mm) and the dimensionless direction components with the same absolute
/// tolerance. In vacuum q/p does not change, so it does not contribute. The
/// step size control follows Hairer, Nørsett, Wanner, Solving Ordinary
/// Differential Equations I, 2nd ed., Section II.4.
///
/// As for @ref HelixStepper, a step to the straight-line distance of a surface
/// can pass it, and the next step goes back. A target aborter must accept an
/// intersection that far behind the track.
///
/// @note The order of a tableau only holds in a smooth field. An interpolated
///       field map is continuous but not differentiable at the cell faces.
///       A step that crosses a face has a local error of order h^2 for every
///       tableau. The number of crossings does not decrease with h, so in a
///       field map every tableau is only accurate to about second order.
/// @note The stepper propagates in vacuum only. It ignores volume material.
class GenericRungeKuttaStepper final {
 public:
  /// Type alias for bound track parameters
  using BoundParameters = BoundTrackParameters;
  /// Type alias for jacobian matrix
  using Jacobian = BoundMatrix;
  /// Type alias for covariance matrix
  using Covariance = BoundMatrix;

  /// Magnetic field and its spatial gradient at one position
  struct FieldAndGradient {
    /// Magnetic field vector
    Vector3 field = Vector3::Zero();
    /// Spatial gradient of the field, with gradient(i, j) = dB_i / dx_j
    SquareMatrix3 gradient = SquareMatrix3::Zero();
  };

  /// Configuration for the Runge-Kutta stepper.
  struct Config {
    /// Magnetic field provider
    std::shared_ptr<const MagneticFieldProvider> bField;

    /// Runge-Kutta tableau. It is shared, so that all steppers can use the
    /// static built-in tableaus without a copy.
    std::shared_ptr<const ButcherTableau> tableau =
        ButcherTableau::dormandPrince54();
  };

  /// Runtime options for Runge-Kutta propagation.
  struct Options : public StepperPlainOptions {
    /// Constructor
    /// @param gctx Geometry context
    /// @param mctx Magnetic field context
    Options(const GeometryContext& gctx, const MagneticFieldContext& mctx)
        : StepperPlainOptions(gctx, mctx) {}

    /// Adapt the step size to the error estimate. The stepper uses the given
    /// step size if this is false or if the tableau has no embedded weights.
    bool adaptiveStepSize = true;

    /// Include the field gradient in the transport jacobian
    bool includeFieldGradient = true;

    /// Distance of the field lookups for the finite-difference gradient.
    ///
    /// The truncation error of the central difference is of order epsilon^2.
    /// The round-off error is of order eps_machine |B| / epsilon, which is
    /// about 1e-14 |B| per mm for 10 um. This is far below the gradient of a
    /// real field. 10 um is also small compared with a field-map cell, so the
    /// lookups usually stay in one cell.
    double fieldGradientEpsilon = 10 * UnitConstants::um;

    /// Include the maximum difference of the jacobian and the embedded
    /// jacobian in the error estimate. The embedded jacobian is
    /// I + h sum_i bEmbedded_i dK_i. Its entries have mixed units, and the
    /// stepper compares all of them with the same tolerance.
    bool jacobianInErrorEstimate = false;

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
    /// Constructor from the options and the field cache
    ///
    /// @param [in] optionsIn is the configuration of the stepper
    /// @param [in] fieldCacheIn is the cache object for the magnetic field
    State(const Options& optionsIn, MagneticFieldProvider::Cache fieldCacheIn)
        : options(optionsIn), fieldCache(std::move(fieldCacheIn)) {}

    /// Configuration options for the stepper
    Options options;

    /// Internal free vector parameters
    FreeVector pars = FreeVector::Zero();

    /// The propagation derivative
    FreeVector derivative = FreeVector::Zero();

    /// Jacobian from local to the global frame
    BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();

    /// Pure transport jacobian part from the Runge-Kutta steps
    FreeMatrix jacTransport = FreeMatrix::Identity();

    /// Covariance matrix for track parameter uncertainties, set if the
    /// covariance is transported
    std::optional<Covariance> cov;

    /// Particle hypothesis
    ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

    /// Adaptive step size of the Runge-Kutta steps
    ConstrainedStep stepSize;

    /// Last performed step (for overstep limit calculation)
    double previousStepSize = 0.;

    /// Magnetic field at the current position, reused as the first stage of
    /// the next step. Reset whenever the position is set from outside.
    std::optional<FieldAndGradient> field;

    /// Whether @ref field contains the gradient
    bool fieldHasGradient = false;

    /// Accumulated path length state
    double pathAccumulated = 0.;

    /// Total number of performed steps
    std::size_t nSteps = 0;

    /// Statistics of the stepper
    StepperStatistics statistics;

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    MagneticFieldProvider::Cache fieldCache;
  };

  /// Constructor requires knowledge of the detector's magnetic field
  /// @param bField The magnetic field provider
  explicit GenericRungeKuttaStepper(
      std::shared_ptr<const MagneticFieldProvider> bField);

  /// @brief Constructor with configuration
  /// @param config The configuration of the stepper
  explicit GenericRungeKuttaStepper(const Config& config);

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
    return detail::updateSingleSurfaceStatus<GenericRungeKuttaStepper>(
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
  bool hasCovariance(const State& state) const { return state.cov.has_value(); }

  /// Get the covariance at the anchor
  ///
  /// The anchor is the frame of the last initialization, transport or update.
  ///
  /// @param state [in] The stepping state (thread-local cache)
  /// @return The covariance at the anchor, or no value if the state does not
  ///         carry a covariance
  const std::optional<Covariance>& covariance(const State& state) const {
    return state.cov;
  }

  /// Set the covariance at the anchor
  ///
  /// The state must already carry a covariance, because the stepper only
  /// transports the Jacobian while it has one.
  ///
  /// @param [in,out] state The stepping state (thread-local cache)
  /// @param [in] covariance The new covariance at the anchor
  /// @throws std::logic_error if the state does not carry a covariance
  void setCovariance(State& state, const Covariance& covariance) const {
    if (!state.cov.has_value()) {
      throw std::logic_error(
          "Cannot set the covariance of a state without a covariance");
    }
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
  /// @param [in, out] state State of the stepper
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
  ///                        if the state carries a covariance
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
  /// @param material is ignored, the stepper propagates in vacuum only
  /// @return the result of the step
  ///
  /// @note The state contains the desired step size. It can be negative during
  ///       backwards track propagation, and since the step size is adaptive,
  ///       it can be modified by the stepper class during propagation.
  Result<double> step(State& state, Direction propDir,
                      const IVolumeMaterial* material) const;

  /// Get the field and optionally its gradient at a position.
  ///
  /// @param [in,out] state The stepper state with the field cache
  /// @param [in] pos The position of the lookup
  /// @param [in] withGradient Whether the gradient is needed
  /// @return The field and the gradient, which is zero without @p withGradient
  Result<FieldAndGradient> getFieldAndGradient(State& state, const Vector3& pos,
                                               bool withGradient) const;

 private:
  /// Magnetic field inside of the detector
  std::shared_ptr<const MagneticFieldProvider> m_bField;

  /// Runge-Kutta tableau
  std::shared_ptr<const ButcherTableau> m_tableau;
};

template <>
struct SupportsBoundParameters<GenericRungeKuttaStepper>
    : public std::true_type {};

}  // namespace Acts
