// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorInitializeArguments.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace Acts {

class IVolumeMaterial;
class Surface;
class TrackingVolume;

namespace Experimental {

/// @brief Navigator that finds surface crossings on the chord of each step
///
/// The standard @c Navigator predicts the surfaces ahead of the track on a
/// straight line and steers the stepper towards them. When the track bends
/// strongly before it reaches a surface, the prediction is wrong: candidates
/// are missed, or a surface that the track never reaches is kept as target.
/// Low-momentum tracks that loop in the magnetic field are the extreme case.
///
/// This navigator does not predict. It lets the stepper take a free step, then
/// intersects the straight segment between the previous and the current
/// position (the chord) with the candidates of the current volume. The chord
/// lies within the sagitta of the step from the integrated trajectory. The
/// crossings found on it become the targets in turn: the stepper steps back
/// to the first one on the real trajectory, the bounds are checked there, and
/// the following ones are reached from there. After the last one, the next
/// free step starts a new chord. A portal crossing switches the volume and
/// starts a new chord at the portal.
///
/// The candidates of a chord come from the navigation policy of the volume,
/// queried with the chord start and direction (and optionally from the chord
/// end looking back), plus the portals of the volume, the extended surfaces
/// of the volume and the additional surfaces of the options. It requires a
/// Gen3 tracking geometry.
///
/// Optionally (on by default), when no crossing is pending, the nearest
/// candidate ahead on the straight line is returned as a step limit. This hint
/// is only a step length: it is re-evaluated after every step, and a reached
/// hint counts only if the chord up to it has no earlier crossing.
///
/// The navigator needs no field information. The bounds tolerance on the
/// chord is derived from the sagitta, using the curvature that the propagator
/// passes to every call.
///
/// @note A chord cannot see a trajectory that crosses a surface and comes back
///       within one step. The depth of such an excursion is below the sagitta
///       of the step, which can be bounded with the maximum step size of the
///       stepper.
/// @note The propagator checks its aborters (path limit, target surface)
///       before it asks the navigator for the next target, so crossings in
///       the step that ends the propagation are not reported. The forward
///       hint makes such a step rare; without it, bound the step size.
class ChordNavigator final {
 public:
  /// Where the candidates of a chord come from
  enum class CandidateSearch {
    /// Navigation policy of the volume, plus its portals
    Policy,
    /// Every surface and portal of the volume (reference, slow)
    AllVolumeSurfaces,
  };

  /// Configuration of the navigator
  struct Config {
    /// The tracking geometry, must be Gen3
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

    /// Where the candidates of a chord come from
    CandidateSearch candidateSearch = CandidateSearch::Policy;

    /// Query the navigation policy a second time from the chord end, looking
    /// back. Covers chords whose straight-line extension meets the policy's
    /// lookup surface closer to the end than to the start.
    bool queryFromChordEnd = true;

    /// Skip the surface search of chords that miss the volume's surfaces.
    ///
    /// At construction, the r-z extent of the surfaces of every volume (all
    /// non-portal surfaces, sensitive and passive) is computed from their
    /// polyhedron representation in @c envelopeGeometryContext. A chord whose
    /// r-z box, padded by the sagitta tolerance, the surface tolerance and
    /// @c envelopeTolerance, misses that extent is tested against the
    /// portals, the extended and the additional surfaces only: the navigation
    /// policy is not queried. This assumes that the policy only offers
    /// portals and surfaces of the volume.
    ///
    /// The forward hint uses the candidates of the last chord, so after a
    /// skipped chord it is taken from the portals only.
    bool envelopePreCheck = false;

    /// Padding of the envelope, for alignment differences to
    /// @c envelopeGeometryContext and approximations of curved bounds
    double envelopeTolerance = 1. * UnitConstants::mm;

    /// Geometry context in which the envelopes are computed
    GeometryContext envelopeGeometryContext =
        GeometryContext::dangerouslyDefaultConstruct();

    /// When no crossing is pending, return the nearest candidate ahead on the
    /// straight line as a step limit. The hint is re-evaluated after every
    /// step and never decides which surfaces are crossed: a hint reached on
    /// the real trajectory counts only if the chord up to it has no earlier
    /// crossing. It makes the steps of straight tracks end on the surfaces,
    /// bounds the steps of steppers without step control, and keeps the step
    /// that ends a propagation at a path limit free of unchecked crossings.
    bool forwardHint = true;

    /// Restart the chord at the current position when the state was moved
    /// after a surface was reached (by an actor, e.g. a Kalman update), and
    /// drop the crossings still pending from before the move. Without it, the
    /// next chord spans the jump and can return crossings the trajectory never
    /// makes. Only meant to be disabled to study that behaviour.
    bool restartChordAfterActorMove = true;

    /// Factor on the sagitta of the chord for the bounds tolerance of a chord
    /// crossing. The tolerance is divided by the incidence cosine, because the
    /// offset between chord and trajectory slides along the surface.
    double sagittaToleranceFactor = 2.;

    /// Smallest incidence cosine used to scale the bounds tolerance
    double minIncidence = 0.05;

    /// Warn once per run when a chord deviates more than this from the
    /// trajectory. The chord then no longer tells where the track went: a
    /// stepper without step control (e.g. the HelixStepper) needs a maximum
    /// step size on curved tracks.
    double maxSagittaWarning = 1. * UnitConstants::mm;

    /// Bounds tolerance of the extended surfaces of the options, used in place
    /// of the tolerance they request. The crossings are found on the real
    /// trajectory, so a measurement surface is reached without relaxing its
    /// bounds. An infinite tolerance would report crossings of the unbounded
    /// plane anywhere in the volume.
    BoundaryTolerance extendedSurfaceTolerance = BoundaryTolerance::None();

    /// Number of candidates to reserve memory for
    std::size_t candidatePreReserve = 50;
  };

  /// Options of the navigator
  struct Options : public NavigatorPlainOptions {
    /// @param gctx The geometry context of the navigation
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// Set the plain navigation options
    /// @param options The plain navigator options to copy
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// Counters of the chord search, for diagnostics and cost studies
  struct ChordStatistics {
    /// Number of chords tested
    std::size_t nChords = 0;
    /// Number of candidates tested on the chords
    std::size_t nCandidates = 0;
    /// Number of crossings found on the chords
    std::size_t nCrossings = 0;
    /// Number of reached surfaces rejected by the bounds check
    std::size_t nRejectedOutOfBounds = 0;
    /// Number of chords whose surface search was skipped by the envelope
    /// pre-check, see Config::envelopePreCheck
    std::size_t nEnvelopeSkips = 0;
  };

  /// A crossing found on a chord
  struct Crossing {
    /// The target, with an infinite tolerance for the approach
    NavigationTarget target = NavigationTarget::None();
    /// Tolerance of the bounds check once the surface is reached
    BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();
    /// Position of the crossing on the chord
    Vector3 chordPosition = Vector3::Zero();
    /// Path length along the chord
    double chordPathLength = 0.;
    /// Bounds tolerance added on the chord for the sagitta
    double slide = 0.;
  };

  /// A surface reached since the last chord search
  struct ReachedSurface {
    /// The surface
    const Surface* surface = nullptr;
    /// Where it was reached
    Vector3 position = Vector3::Zero();
  };

  /// An extended surface with its volume resolved
  struct ResolvedExtendedSurface {
    /// The surface
    const Surface* surface = nullptr;
    /// Volume that holds the surface
    const TrackingVolume* volume = nullptr;
  };

  /// An additional surface of the options with its bookkeeping
  struct AdditionalSurfaceState {
    /// The entry of the options
    const AdditionalSurface* entry = nullptr;
    /// Whether the propagation reached the surface
    bool reached = false;
  };

  /// Navigation state, one per propagation
  struct State {
    /// @param options_ The navigation options
    explicit State(const Options& options_) : options(options_) {}

    /// Navigation options
    Options options;

    /// Surface the navigation started on
    const Surface* startSurface = nullptr;
    /// Surface the navigation targets
    const Surface* targetSurface = nullptr;
    /// Surface the propagation is on, if any
    const Surface* currentSurface = nullptr;
    /// Current tracking volume
    const TrackingVolume* currentVolume = nullptr;

    /// Latest geometric curvature passed by the propagator
    Vector3 curvature = Vector3::Zero();

    /// Policy states of the volumes
    NavigationPolicyStateManager policyStateManager;

    /// Candidates of the last chord, grouped by surface
    NavigationStream candidates;
    /// Volume the candidates were collected in
    const TrackingVolume* candidatesVolume = nullptr;

    /// The step limit handed out last, see Config::forwardHint
    NavigationTarget hint = NavigationTarget::None();
    /// Whether the current target is the hint
    bool hintActive = false;

    /// Start of the chord of the current free step
    std::optional<Vector3> chordStart;
    /// Position at the last reached surface, before the actors ran. No step
    /// happens until the next target is asked for, so a different position
    /// then means an actor (e.g. a Kalman update) moved the state.
    std::optional<Vector3> reachedPosition;

    /// Crossings found on the last chord, sorted along the chord
    std::vector<Crossing> crossings;
    /// Index of the next crossing to hand out
    std::size_t nextCrossing = 0;
    /// Surfaces reached since the last chord search, so the next chord does
    /// not count them again
    std::vector<ReachedSurface> recentlyReached;

    /// Extended surfaces resolved at initialization
    std::vector<ResolvedExtendedSurface> extendedSurfaces;
    /// Additional surfaces of the options
    std::vector<AdditionalSurfaceState> additionalSurfaces;

    /// Flag to break the navigation
    bool navigationBreak = false;
    /// Whether the large-sagitta warning was issued in this run
    bool sagittaWarned = false;

    /// Navigation statistics
    NavigatorStatistics statistics;
    /// Chord search statistics
    ChordStatistics chordStatistics;
  };

  /// @param cfg The navigator configuration
  /// @param logger The logger
  /// @throws std::invalid_argument without a Gen3 tracking geometry
  explicit ChordNavigator(Config cfg,
                          std::shared_ptr<const Logger> logger =
                              getDefaultLogger("ChordNavigator",
                                               Logging::Level::INFO));

  /// @param options The navigation options
  /// @return A new navigation state
  State makeState(const Options& options) const;

  /// @param state The navigation state
  /// @return The surface the propagation is on, if any
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// @param state The navigation state
  /// @return The current tracking volume
  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  /// @param state The navigation state
  /// @return The material of the current volume, if any
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const;

  /// @param state The navigation state
  /// @return The surface the navigation started on
  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  /// @param state The navigation state
  /// @return The surface the navigation targets
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  /// @param state The navigation state
  /// @return Whether the propagation left the world
  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  /// @param state The navigation state
  /// @return Whether the navigation is broken
  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// Initialize the navigation state for a run
  /// @param state The navigation state
  /// @param args The initialization arguments of this run
  /// @return Success, or the reason the start volume could not be resolved
  [[nodiscard]] Result<void> initialize(
      State& state, const NavigatorInitializeArguments& args) const;

  /// Get the next target
  ///
  /// After a free step, the chord from the previous position is checked and
  /// the first crossing on it is returned. Without a crossing, a new chord
  /// starts at @p position and no target is returned, so the stepper takes a
  /// free step.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param curvature The current geometric curvature (inverse length)
  /// @return The next target, or none for a free step
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction,
                              const Vector3& curvature = Vector3::Zero()) const;

  /// Check whether the current target is still valid
  ///
  /// A crossing found on a chord stays the target until it is reached or the
  /// propagator reports it unreachable. Dropping it early would lose it,
  /// because the next chord starts past it. The forward hint (see
  /// Config::forwardHint) only bounds one step and is never kept.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param curvature The current geometric curvature (inverse length)
  /// @return False if the current target is the forward hint, else true
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction,
                        const Vector3& curvature = Vector3::Zero()) const;

  /// Handle a reached target
  ///
  /// Checks the bounds of a reached surface at the real position and switches
  /// the volume at a portal. The next chord starts here.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param surface The reached surface
  /// @param curvature The current geometric curvature (inverse length)
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction, const Surface& surface,
                            const Vector3& curvature = Vector3::Zero()) const;

 private:
  /// The r-z extent of the non-portal surfaces of a volume
  struct VolumeEnvelope {
    double rMin = std::numeric_limits<double>::max();
    double rMax = std::numeric_limits<double>::lowest();
    double zMin = std::numeric_limits<double>::max();
    double zMax = std::numeric_limits<double>::lowest();

    /// @return Whether the volume holds no surface
    bool empty() const { return rMin > rMax; }
  };

  /// Whether the chord from @p start to @p end, padded by @p padding, misses
  /// the envelope of @p volume. False without envelope pre-check.
  bool missesEnvelope(const TrackingVolume& volume, const Vector3& start,
                      const Vector3& end, double padding) const;

  /// Collect the candidates of the current volume for a chord
  ///
  /// @param padding Distance by which the trajectory can deviate from the
  ///        chord, for the envelope pre-check
  void collectCandidates(State& state, const Vector3& start, const Vector3& end,
                         const Vector3& direction, bool queryFromEnd,
                         double padding) const;

  /// Nearest candidate ahead on the straight line, see Config::forwardHint
  NavigationTarget forwardHint(State& state, const Vector3& position,
                               const Vector3& direction) const;

  /// Handle a reached crossing or hint: switch the volume at a portal, check
  /// the bounds of a surface
  void reachSurface(State& state, const Vector3& position,
                    const Vector3& direction, const NavigationTarget& target,
                    const BoundaryTolerance& tolerance) const;

  /// Find the crossings of the chord from @p start to @p end
  void findCrossings(State& state, const Vector3& start,
                     const Vector3& end) const;

  /// Enter @p volume, creating its policy state
  void enterVolume(State& state, const TrackingVolume* volume,
                   const Vector3& position, const Vector3& direction) const;

  /// Leave the current volume, removing its policy state
  void leaveVolume(State& state) const;

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::shared_ptr<const Logger> m_logger;

  /// Envelopes of the volumes, filled with Config::envelopePreCheck
  std::unordered_map<const TrackingVolume*, VolumeEnvelope> m_envelopes;
};

}  // namespace Experimental
}  // namespace Acts
