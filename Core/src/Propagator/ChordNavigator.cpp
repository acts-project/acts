// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/ChordNavigator.hpp"

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Propagator/NavigatorError.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>
#include <string>

namespace Acts::Experimental {

namespace {

std::string volInfo(const ChordNavigator::State& state) {
  return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                         : "No Volume") +
         " | ";
}

/// Widen a bounds tolerance by an absolute distance
BoundaryTolerance widen(const BoundaryTolerance& tolerance, double extra) {
  if (tolerance.isInfinite() || extra <= 0.) {
    return tolerance;
  }
  if (tolerance.hasAbsoluteEuclidean()) {
    return BoundaryTolerance::AbsoluteEuclidean(
        tolerance.asAbsoluteEuclidean().tolerance + extra);
  }
  // None, and the chi2 tolerances which have no absolute scale to add to
  return BoundaryTolerance::AbsoluteEuclidean(extra);
}

}  // namespace

ChordNavigator::ChordNavigator(Config cfg,
                               std::shared_ptr<const Logger> navLogger)
    : m_cfg{std::move(cfg)}, m_logger{std::move(navLogger)} {
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("ChordNavigator: no tracking geometry");
  }
  if (m_cfg.trackingGeometry->geometryVersion() !=
      TrackingGeometry::GeometryVersion::Gen3) {
    throw std::invalid_argument(
        "ChordNavigator: requires a Gen3 tracking geometry");
  }
  if (m_cfg.sagittaToleranceFactor < 0. || m_cfg.minIncidence <= 0.) {
    throw std::invalid_argument(
        "ChordNavigator: sagittaToleranceFactor must be >= 0 and "
        "minIncidence > 0");
  }
  if (m_cfg.envelopeTolerance < 0.) {
    throw std::invalid_argument(
        "ChordNavigator: envelopeTolerance must be >= 0");
  }

  if (m_cfg.envelopePreCheck) {
    const GeometryContext& gctx = m_cfg.envelopeGeometryContext;
    m_cfg.trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
      VolumeEnvelope envelope;
      for (const Surface& surface : volume->surfaces()) {
        const Extent extent = surface.polyhedronRepresentation(gctx).extent();
        envelope.rMin =
            std::min(envelope.rMin, extent.min(AxisDirection::AxisR));
        envelope.rMax =
            std::max(envelope.rMax, extent.max(AxisDirection::AxisR));
        envelope.zMin =
            std::min(envelope.zMin, extent.min(AxisDirection::AxisZ));
        envelope.zMax =
            std::max(envelope.zMax, extent.max(AxisDirection::AxisZ));
      }
      ACTS_VERBOSE("Envelope of " << volume->volumeName() << ": r ["
                                  << envelope.rMin << ", " << envelope.rMax
                                  << "], z [" << envelope.zMin << ", "
                                  << envelope.zMax << "]");
      m_envelopes.emplace(volume, envelope);
    });
  }
}

bool ChordNavigator::missesEnvelope(const TrackingVolume& volume,
                                    const Vector3& start, const Vector3& end,
                                    double padding) const {
  if (!m_cfg.envelopePreCheck) {
    return false;
  }
  const auto found = m_envelopes.find(&volume);
  if (found == m_envelopes.end()) {
    // Not in the geometry the envelopes were built from: no skipping
    return false;
  }
  const VolumeEnvelope& envelope = found->second;
  if (envelope.empty()) {
    return true;
  }

  // r-z box of the chord: z between the end points, r between the closest
  // approach of the segment to the z axis and the farther end point
  const double zMin = std::min(start.z(), end.z());
  const double zMax = std::max(start.z(), end.z());
  if (zMax < envelope.zMin - padding || zMin > envelope.zMax + padding) {
    return true;
  }
  const Vector2 p0 = start.head<2>();
  const Vector2 d = end.head<2>() - p0;
  const double rMax = std::max(p0.norm(), end.head<2>().norm());
  const double d2 = d.squaredNorm();
  const double t = d2 > 0. ? std::clamp(-p0.dot(d) / d2, 0., 1.) : 0.;
  const double rMin = (p0 + t * d).norm();
  return rMax < envelope.rMin - padding || rMin > envelope.rMax + padding;
}

ChordNavigator::State ChordNavigator::makeState(const Options& options) const {
  State state(options);
  state.candidates.reserve(m_cfg.candidatePreReserve);
  return state;
}

const IVolumeMaterial* ChordNavigator::currentVolumeMaterial(
    const State& state) const {
  if (state.currentVolume == nullptr) {
    return nullptr;
  }
  return state.currentVolume->volumeMaterial();
}

Result<void> ChordNavigator::initialize(
    State& state, const NavigatorInitializeArguments& args) const {
  const GeometryContext& gctx = state.options.geoContext;
  const Vector3& position = args.position;
  const Vector3& direction = args.direction;

  state.startSurface = args.startSurface;
  state.targetSurface = args.targetSurface;
  state.currentSurface = nullptr;
  state.currentVolume = nullptr;
  state.curvature = args.curvature;
  state.policyStateManager.reset();
  state.candidates.reset();
  state.chordStart.reset();
  state.reachedPosition.reset();
  state.crossings.clear();
  state.nextCrossing = 0;
  state.recentlyReached.clear();
  state.hintActive = false;
  state.candidatesVolume = nullptr;
  state.navigationBreak = false;
  state.sagittaWarned = false;

  // Extended surfaces are scoped to the volume that holds them
  state.extendedSurfaces.clear();
  for (const ExtendedSurface& extended : state.options.extendedSurfaces) {
    if (extended.surface == nullptr) {
      throw std::invalid_argument(
          "ChordNavigator: extended surface is nullptr");
    }
    const GeometryIdentifier geoId = extended.surface->geometryId();
    if (m_cfg.trackingGeometry->findSurface(geoId) != extended.surface) {
      throw std::invalid_argument(
          std::format("ChordNavigator: the extended surface {} is not part of "
                      "the tracking geometry",
                      geoId));
    }
    const TrackingVolume* volume = m_cfg.trackingGeometry->findVolume(
        geoId.withSensitive(0).withBoundary(0));
    if (volume == nullptr) {
      throw std::invalid_argument(std::format(
          "ChordNavigator: no volume found for the extended surface {}",
          geoId));
    }
    state.extendedSurfaces.push_back({extended.surface, volume});
  }

  state.additionalSurfaces.clear();
  for (const AdditionalSurface& additional : state.options.additionalSurfaces) {
    if (additional.surface == nullptr) {
      throw std::invalid_argument(
          "ChordNavigator: additional surface is nullptr");
    }
    state.additionalSurfaces.push_back({&additional, false});
  }

  if (state.startSurface != nullptr &&
      !state.startSurface->isOnSurface(gctx, position, direction,
                                       BoundaryTolerance::Infinite(),
                                       state.options.surfaceTolerance)) {
    ACTS_DEBUG("Not on the start surface " << state.startSurface->geometryId()
                                           << " at " << toString(position));
    return NavigatorError::NotOnExpectedSurface;
  }

  const TrackingVolume* startVolume = args.startVolume;
  if (startVolume == nullptr) {
    const auto resolved = m_cfg.trackingGeometry->resolveLowestTrackingVolume(
        gctx, position, direction, state.startSurface,
        state.options.surfaceTolerance);
    if (!resolved.ok()) {
      ACTS_DEBUG(
          "Could not resolve the start volume: " << resolved.error().message());
      state.navigationBreak = true;
      return resolved.error();
    }
    startVolume = *resolved;
  }
  if (startVolume == nullptr) {
    ACTS_DEBUG("No start volume at " << toString(position));
    state.navigationBreak = true;
    return NavigatorError::NoStartVolume;
  }

  enterVolume(state, startVolume, position, direction);
  state.currentSurface = state.startSurface;

  // The first free step starts a chord here
  state.chordStart = position;

  ACTS_VERBOSE(volInfo(state) << "Initialized at " << toString(position));

  return Result<void>::success();
}

NavigationTarget ChordNavigator::nextTarget(State& state,
                                            const Vector3& position,
                                            const Vector3& direction,
                                            const Vector3& curvature) const {
  state.curvature = curvature;
  state.currentSurface = nullptr;

  if (state.navigationBreak || state.currentVolume == nullptr) {
    return NavigationTarget::None();
  }

  state.hintActive = false;

  // No step is taken between a reached surface and this call, so a moved
  // state was moved by an actor (e.g. a Kalman update). Searching a chord
  // across that jump would find crossings the trajectory never makes, and the
  // crossings still pending belong to the trajectory before the jump: they
  // are dropped, the chord restarts here, and the reached surface is recorded
  // where the state is now.
  if (state.reachedPosition.has_value()) {
    const Vector3 reached = *state.reachedPosition;
    state.reachedPosition.reset();
    const double moved = (position - reached).norm();
    if (m_cfg.restartChordAfterActorMove &&
        moved > state.options.surfaceTolerance) {
      ACTS_VERBOSE(volInfo(state)
                   << "State moved by " << moved
                   << " mm since the surface was reached, the chord restarts "
                      "here and "
                   << state.crossings.size() - state.nextCrossing
                   << " pending crossing(s) are dropped");
      state.crossings.clear();
      state.nextCrossing = 0;
      state.chordStart = position;
      for (ReachedSurface& entry : state.recentlyReached) {
        if ((entry.position - reached).norm() <=
            state.options.surfaceTolerance) {
          entry.position = position;
        }
      }
    }
  }

  // A step ended here: check its chord
  if (state.chordStart.has_value() &&
      (position - *state.chordStart).norm() > state.options.surfaceTolerance) {
    const Vector3 start = *state.chordStart;
    state.chordStart.reset();
    findCrossings(state, start, position);
  }

  // Hand out the crossings in order. The propagator asks again after a reached
  // crossing, or after one it dropped as unreachable.
  const GeometryContext& gctx = state.options.geoContext;
  while (state.nextCrossing < state.crossings.size()) {
    const Crossing& crossing = state.crossings[state.nextCrossing++];
    NavigationTarget target = crossing.target;

    // The propagator intersects the target from the current position along
    // the current direction; pick the solution at the chord crossing
    const auto intersections = target.surface().intersect(
        gctx, position, direction, BoundaryTolerance::Infinite(),
        state.options.surfaceTolerance);
    std::optional<std::size_t> best;
    double bestDistance = std::numeric_limits<double>::max();
    for (const auto [index, intersection] : enumerate(intersections)) {
      if (!intersection.isValid()) {
        continue;
      }
      const double distance =
          (intersection.position() - crossing.chordPosition).norm();
      if (distance < bestDistance) {
        bestDistance = distance;
        best = index;
      }
    }
    if (!best.has_value()) {
      ACTS_VERBOSE(volInfo(state)
                   << "Crossing of " << target.surface().geometryId()
                   << " has no straight-line solution, skip");
      continue;
    }
    target.intersection() = intersections.at(*best);
    target.intersectionIndex() = static_cast<IntersectionIndex>(*best);

    ACTS_VERBOSE(volInfo(state)
                 << "Target set to crossing " << target.surface().geometryId()
                 << " at path " << target.pathLength());
    return target;
  }

  // Nothing pending: the next step starts a new chord here
  state.crossings.clear();
  state.nextCrossing = 0;
  if (!state.chordStart.has_value()) {
    state.chordStart = position;
  }

  if (!m_cfg.forwardHint) {
    return NavigationTarget::None();
  }
  return forwardHint(state, position, direction);
}

NavigationTarget ChordNavigator::forwardHint(State& state,
                                             const Vector3& position,
                                             const Vector3& direction) const {
  const GeometryContext& gctx = state.options.geoContext;
  const double surfaceTolerance = state.options.surfaceTolerance;

  // The candidates of the last chord lie around here; collect new ones after
  // a volume switch or at the start
  if (state.candidatesVolume != state.currentVolume) {
    collectCandidates(state, position, position, direction, false, 0.);
  }

  NavigationTarget hint = NavigationTarget::None();
  const Surface* previous = nullptr;
  for (const NavigationTarget& candidate : state.candidates.candidates()) {
    const Surface& surface = candidate.surface();
    if (&surface == previous) {
      continue;
    }
    previous = &surface;
    const auto intersections =
        surface.intersect(gctx, position, direction,
                          candidate.boundaryTolerance(), surfaceTolerance);
    for (const auto [index, intersection] : enumerate(intersections)) {
      if (!intersection.isValid() ||
          intersection.pathLength() <= surfaceTolerance ||
          (!hint.isNone() && intersection.pathLength() >= hint.pathLength())) {
        continue;
      }
      const auto intersectionIndex = static_cast<IntersectionIndex>(index);
      hint = candidate.isPortalTarget()
                 ? NavigationTarget(intersection, intersectionIndex,
                                    candidate.portal(),
                                    candidate.boundaryTolerance())
                 : NavigationTarget(intersection, intersectionIndex, surface,
                                    candidate.boundaryTolerance());
    }
  }

  if (!hint.isNone()) {
    state.hint = hint;
    state.hintActive = true;
    ACTS_VERBOSE(volInfo(state)
                 << "Step limited by the hint " << hint.surface().geometryId()
                 << " at path " << hint.pathLength());
  }
  return hint;
}

bool ChordNavigator::checkTargetValid(const State& state,
                                      const Vector3& /*position*/,
                                      const Vector3& /*direction*/,
                                      const Vector3& /*curvature*/) const {
  // A hint only bounds one step, a crossing stays until reached or dropped
  return !state.hintActive;
}

void ChordNavigator::handleSurfaceReached(State& state, const Vector3& position,
                                          const Vector3& direction,
                                          const Surface& surface,
                                          const Vector3& curvature) const {
  state.curvature = curvature;
  if (state.navigationBreak || state.currentVolume == nullptr) {
    return;
  }

  state.reachedPosition = position;

  if (state.hintActive) {
    state.hintActive = false;
    // The hint was predicted on a straight line: check the chord up to here
    // for crossings before it
    if (state.chordStart.has_value() && (position - *state.chordStart).norm() >
                                            state.options.surfaceTolerance) {
      const Vector3 start = *state.chordStart;
      state.chordStart.reset();
      findCrossings(state, start, position);
      std::erase_if(state.crossings, [&](const Crossing& crossing) {
        return &crossing.target.surface() == &surface &&
               (crossing.chordPosition - position).norm() <=
                   crossing.slide + 2. * state.options.surfaceTolerance;
      });
      if (!state.crossings.empty()) {
        // Step back to them first; this surface is found again on the way
        ACTS_VERBOSE(volInfo(state)
                     << "Hint " << surface.geometryId() << " reached after "
                     << state.crossings.size() << " other crossing(s)");
        return;
      }
    }
    const NavigationTarget hint = state.hint;
    reachSurface(state, position, direction, hint, hint.boundaryTolerance());
    if (!state.chordStart.has_value()) {
      state.chordStart = position;
    }
    return;
  }

  // Otherwise the reached target is the crossing handed out last
  if (state.nextCrossing > 0 && state.nextCrossing <= state.crossings.size() &&
      &state.crossings[state.nextCrossing - 1].target.surface() == &surface) {
    const Crossing crossing = state.crossings[state.nextCrossing - 1];
    reachSurface(state, position, direction, crossing.target,
                 crossing.boundaryTolerance);
  } else {
    ACTS_VERBOSE(volInfo(state) << "Reached " << surface.geometryId()
                                << ", which is not a pending target");
  }
  // The remaining crossings of the chord are handed out next. Near-coincident
  // crossings can be ordered differently on the chord than on the trajectory,
  // so the next one may lie slightly behind; the stepper then steps back to it.
}

void ChordNavigator::reachSurface(State& state, const Vector3& position,
                                  const Vector3& direction,
                                  const NavigationTarget& target,
                                  const BoundaryTolerance& tolerance) const {
  const GeometryContext& gctx = state.options.geoContext;
  const Surface& surface = target.surface();

  if (target.isPortalTarget()) {
    state.currentSurface = &surface;
    const auto resolved =
        target.portal().resolveVolume(gctx, position, direction);
    if (!resolved.ok()) {
      ACTS_ERROR(volInfo(state) << "Failed to resolve the volume through "
                                << "portal " << surface.geometryId() << ": "
                                << resolved.error().message());
      state.navigationBreak = true;
      return;
    }
    leaveVolume(state);
    enterVolume(state, *resolved, position, direction);
    ++state.statistics.nVolumeSwitches;
    if (state.currentVolume == nullptr) {
      ACTS_VERBOSE("Left the world through portal " << surface.geometryId());
      state.navigationBreak = true;
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Entered volume through portal " << surface.geometryId());
    }
    // The other crossings belong to the volume just left; the next chord
    // starts at the portal
    state.crossings.clear();
    state.nextCrossing = 0;
    state.recentlyReached.clear();
    state.chordStart = position;
    return;
  }

  // The step lands within the surface tolerance, which also blurs the position
  // on the surface: a crossing on a shared module edge belongs to both modules
  state.recentlyReached.push_back({&surface, position});
  if (surface.isOnSurface(gctx, position, direction,
                          widen(tolerance, state.options.surfaceTolerance),
                          state.options.surfaceTolerance)) {
    state.currentSurface = &surface;
    ACTS_VERBOSE(volInfo(state) << "Reached surface " << surface.geometryId());
    for (AdditionalSurfaceState& additional : state.additionalSurfaces) {
      if (additional.entry->surface == &surface) {
        additional.reached = true;
      }
    }
  } else {
    ++state.chordStatistics.nRejectedOutOfBounds;
    ACTS_VERBOSE(volInfo(state) << "Reached surface " << surface.geometryId()
                                << " outside its bounds, no hit");
  }
}

void ChordNavigator::collectCandidates(State& state, const Vector3& start,
                                       const Vector3& end,
                                       const Vector3& direction,
                                       bool queryFromEnd,
                                       double padding) const {
  const GeometryContext& gctx = state.options.geoContext;
  const TrackingVolume& volume = *state.currentVolume;

  state.candidates.reset();

  // Envelope pre-check: a chord that misses all surfaces of the volume needs
  // only the portals, the extended and the additional surfaces
  const bool searchSurfaces = !missesEnvelope(
      volume, start, end,
      padding + state.options.surfaceTolerance + m_cfg.envelopeTolerance);
  if (!searchSurfaces) {
    ++state.chordStatistics.nEnvelopeSkips;
    ACTS_VERBOSE(volInfo(state)
                 << "Chord misses the surface envelope, portals only");
  }

  // Extended and additional surfaces first: the de-duplication below keeps the
  // first entry of a surface
  for (const ResolvedExtendedSurface& extended : state.extendedSurfaces) {
    if (extended.volume == &volume) {
      state.candidates.addSurfaceCandidate(*extended.surface,
                                           m_cfg.extendedSurfaceTolerance);
    }
  }
  for (const AdditionalSurfaceState& additional : state.additionalSurfaces) {
    if (additional.reached && additional.entry->dropAfterReached) {
      continue;
    }
    if (additional.entry->volume != nullptr &&
        additional.entry->volume != &volume) {
      continue;
    }
    state.candidates.addSurfaceCandidate(*additional.entry->surface,
                                         additional.entry->boundaryTolerance);
  }

  // Portals always, so the volume can be left whatever the policy offers
  for (const Portal& portal : volume.portals()) {
    state.candidates.addPortalCandidate(portal);
  }

  if (!searchSurfaces) {
    // only the candidates above
  } else if (m_cfg.candidateSearch == CandidateSearch::AllVolumeSurfaces) {
    for (const Surface& surface : volume.surfaces()) {
      state.candidates.addSurfaceCandidate(surface, BoundaryTolerance::None());
    }
  } else if (volume.navigationPolicy() != nullptr) {
    AppendOnlyNavigationStream appendOnly{state.candidates};
    NavigationPolicyState policyState = state.policyStateManager.currentState();
    NavigationArguments args{.position = start,
                             .direction = direction,
                             .curvature = state.curvature,
                             .tolerance = BoundaryTolerance::None()};
    volume.initializeNavigationCandidates(gctx, args, policyState, appendOnly,
                                          logger());
    if (queryFromEnd) {
      args.position = end;
      args.direction = -direction;
      volume.initializeNavigationCandidates(gctx, args, policyState, appendOnly,
                                            logger());
    }
  }

  // A surface can be offered more than once (both policy queries, extended
  // surfaces); sorting groups its entries and keeps the first one first
  std::ranges::stable_sort(
      state.candidates.candidates(), std::less{},
      [](const NavigationTarget& candidate) { return &candidate.surface(); });
  state.candidatesVolume = &volume;
}

void ChordNavigator::findCrossings(State& state, const Vector3& start,
                                   const Vector3& end) const {
  const GeometryContext& gctx = state.options.geoContext;
  const double surfaceTolerance = state.options.surfaceTolerance;

  state.crossings.clear();
  state.nextCrossing = 0;

  const Vector3 chord = end - start;
  const double length = chord.norm();
  const Vector3 direction = chord / length;
  // Largest distance between the chord and the trajectory of the step
  const double sagitta = 0.125 * length * length * state.curvature.norm();
  if (sagitta > m_cfg.maxSagittaWarning && !state.sagittaWarned) {
    state.sagittaWarned = true;
    ACTS_WARNING(volInfo(state)
                 << "Chord of " << length << " mm deviates up to " << sagitta
                 << " mm from the trajectory; bound the step size of the "
                    "stepper");
  }

  collectCandidates(state, start, end, direction, m_cfg.queryFromChordEnd,
                    m_cfg.sagittaToleranceFactor * sagitta);

  ++state.chordStatistics.nChords;

  // Entries of the same surface are adjacent, keep the first one
  const Surface* previous = nullptr;
  for (const NavigationTarget& candidate : state.candidates.candidates()) {
    const Surface& surface = candidate.surface();
    if (&surface == previous) {
      continue;
    }
    previous = &surface;
    ++state.chordStatistics.nCandidates;

    const auto intersections =
        surface.intersect(gctx, start, direction, BoundaryTolerance::Infinite(),
                          surfaceTolerance);
    for (const auto [index, intersection] : enumerate(intersections)) {
      if (!intersection.isValid()) {
        continue;
      }
      const double pathLength = intersection.pathLength();
      // The surface the chord starts on is excluded by the near limit
      if (pathLength <= surfaceTolerance ||
          pathLength > length + surfaceTolerance) {
        continue;
      }

      // The trajectory crosses the surface within the sagitta of the chord,
      // slid along the surface by the incidence
      const Vector3 normal =
          surface.normal(gctx, intersection.position(), direction);
      const double incidence =
          std::max(std::abs(normal.dot(direction)), m_cfg.minIncidence);
      const double slide = m_cfg.sagittaToleranceFactor * sagitta / incidence;
      // A crossing handled since the last chord, at the same place, is not
      // counted again
      const bool handled = std::ranges::any_of(
          state.recentlyReached, [&](const ReachedSurface& reached) {
            return reached.surface == &surface &&
                   (reached.position - intersection.position()).norm() <=
                       slide + s_onSurfaceTolerance + surfaceTolerance;
          });
      if (handled) {
        continue;
      }
      const BoundaryTolerance chordTolerance =
          widen(candidate.boundaryTolerance(), slide);
      if (!chordTolerance.isInfinite()) {
        const auto local =
            surface.globalToLocal(gctx, intersection.position(), direction,
                                  std::numeric_limits<double>::max());
        if (!local.ok() || !surface.insideBounds(*local, chordTolerance)) {
          continue;
        }
      }

      Crossing crossing;
      if (candidate.isPortalTarget()) {
        crossing.target = NavigationTarget(
            intersection, static_cast<IntersectionIndex>(index),
            candidate.portal(), BoundaryTolerance::Infinite());
      } else {
        crossing.target = NavigationTarget(
            intersection, static_cast<IntersectionIndex>(index), surface,
            BoundaryTolerance::Infinite());
      }
      crossing.boundaryTolerance = candidate.boundaryTolerance();
      crossing.chordPosition = intersection.position();
      crossing.chordPathLength = pathLength;
      crossing.slide = slide;
      state.crossings.push_back(crossing);
    }
  }

  state.recentlyReached.clear();
  // Crossings at exactly the same path length are left in arbitrary order:
  // such ties need coincident surfaces, which a valid geometry does not have
  std::ranges::sort(state.crossings, std::less{}, &Crossing::chordPathLength);
  state.chordStatistics.nCrossings += state.crossings.size();

  ACTS_VERBOSE(volInfo(state)
               << "Chord of " << length << " mm, sagitta " << sagitta << " mm: "
               << state.chordStatistics.nCandidates << " candidates so far, "
               << state.crossings.size() << " crossings");
}

void ChordNavigator::enterVolume(State& state, const TrackingVolume* volume,
                                 const Vector3& position,
                                 const Vector3& direction) const {
  state.currentVolume = volume;
  if (volume == nullptr) {
    return;
  }
  const INavigationPolicy* policy = volume->navigationPolicy();
  if (policy != nullptr && !policy->isStateless()) {
    policy->createState(state.options.geoContext,
                        {.position = position,
                         .direction = direction,
                         .curvature = state.curvature},
                        state.policyStateManager, logger());
  }
}

void ChordNavigator::leaveVolume(State& state) const {
  if (state.currentVolume == nullptr) {
    return;
  }
  const INavigationPolicy* policy = state.currentVolume->navigationPolicy();
  if (policy != nullptr && !policy->isStateless()) {
    policy->popState(state.policyStateManager, logger());
  }
}

}  // namespace Acts::Experimental
