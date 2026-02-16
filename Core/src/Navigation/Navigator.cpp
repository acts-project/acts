// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/Navigator.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Propagator/NavigatorError.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <cassert>
#include <sstream>

namespace Acts {

Navigator::Navigator(Config cfg, std::shared_ptr<const Logger> _logger)
    : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Navigator: No tracking geometry provided.");
  }
  m_geometryVersion = m_cfg.trackingGeometry->geometryVersion();
}

Navigator::State Navigator::makeState(const Options& options) const {
  State state(options);
  return state;
}

const Surface* Navigator::currentSurface(const State& state) const {
  return state.currentSurface;
}

const TrackingVolume* Navigator::currentVolume(const State& state) const {
  return state.currentVolume;
}

const IVolumeMaterial* Navigator::currentVolumeMaterial(
    const State& state) const {
  if (state.currentVolume == nullptr) {
    return nullptr;
  }
  return state.currentVolume->volumeMaterial();
}

const Surface* Navigator::startSurface(const State& state) const {
  return state.startSurface;
}

const Surface* Navigator::targetSurface(const State& state) const {
  return state.targetSurface;
}

bool Navigator::endOfWorldReached(const State& state) const {
  return state.currentVolume == nullptr;
}

bool Navigator::navigationBreak(const State& state) const {
  return state.navigationBreak;
}

Result<void> Navigator::initialize(State& state, const Vector3& position,
                                   const Vector3& direction,
                                   Direction propagationDirection) const {
  static_cast<void>(propagationDirection);

  ACTS_VERBOSE(volInfo(state) << "Initialization.");

  auto printGeometryVersion = [](auto ver) {
    using enum TrackingGeometry::GeometryVersion;
    switch (ver) {
      case Gen1:
        return "Gen1";
      case Gen3:
        return "Gen3";
      default:
        throw std::runtime_error("Unknown geometry version.");
    }
  };
  ACTS_VERBOSE(volInfo(state) << "Geometry version is: "
                              << printGeometryVersion(m_geometryVersion));

  state.reset();

  if (m_geometryVersion == GeometryVersion::Gen3) {
    // Empirical pre-allocation of candidates for the next navigation
    // iteration.
    // @TODO: Make this user configurable through the configuration
    state.stream.candidates().reserve(50);

    state.freeCandidates.clear();
    state.freeCandidates.reserve(state.options.freeSurfaces.size());
    for (const Surface* candidate : state.options.freeSurfaces) {
      state.freeCandidates.emplace_back(candidate, false);
    }
  }

  state.startSurface = state.options.startSurface;
  state.targetSurface = state.options.targetSurface;

  // @TODO: Implement fast initialization with Gen3. This requires the volume
  // lookup to work properly

  // Fast Navigation initialization for start condition:
  // - short-cut through object association, saves navigation in the
  // - geometry and volume tree search for the lowest volume
  if (state.startSurface != nullptr &&
      state.startSurface->associatedLayer() != nullptr) {
    ACTS_VERBOSE(
        volInfo(state)
        << "Fast start initialization through association from Surface.");

    state.startLayer = state.startSurface->associatedLayer();
    state.startVolume = state.startLayer->trackingVolume();
  } else if (state.startVolume != nullptr) {
    ACTS_VERBOSE(
        volInfo(state)
        << "Fast start initialization through association from Volume.");

    state.startLayer =
        state.startVolume->associatedLayer(state.options.geoContext, position);
  } else {
    ACTS_VERBOSE(volInfo(state) << "Slow start initialization through search.");
    ACTS_VERBOSE(volInfo(state)
                 << "Starting from position " << toString(position)
                 << " and direction " << toString(direction));

    // current volume and layer search through global search
    state.startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
        state.options.geoContext, position);

    if (state.startVolume != nullptr) {
      state.startLayer = state.startVolume->associatedLayer(
          state.options.geoContext, position);
    } else {
      ACTS_DEBUG(volInfo(state)
                 << "No start volume resolved. Nothing left to do.");
      state.navigationBreak = true;
      return Result<void>::failure(NavigatorError::NoStartVolume);
    }
  }

  state.currentVolume = state.startVolume;
  state.currentLayer = state.startLayer;
  state.currentSurface = state.startSurface;

  if (state.currentVolume != nullptr) {
    ACTS_VERBOSE(volInfo(state) << "Start volume resolved "
                                << state.currentVolume->geometryId());

    if (!state.currentVolume->inside(state.options.geoContext, position,
                                     state.options.surfaceTolerance)) {
      ACTS_DEBUG(volInfo(state)
                 << "We did not end up inside the expected volume. position = "
                 << position.transpose());

      return Result<void>::failure(NavigatorError::NotInsideExpectedVolume);
    }

    if (const auto* policy = state.currentVolume->navigationPolicy();
        policy != nullptr) {
      ACTS_VERBOSE(volInfo(state)
                   << "Creating initial navigation policy state for volume.");
      policy->createState(state.options.geoContext,
                          {.position = position, .direction = direction},
                          state.policyStateManager, logger());
    }
  }
  if (state.currentLayer != nullptr) {
    ACTS_VERBOSE(volInfo(state) << "Start layer resolved "
                                << state.currentLayer->geometryId());
  }
  if (state.currentSurface != nullptr) {
    ACTS_VERBOSE(volInfo(state) << "Start surface resolved "
                                << state.currentSurface->geometryId());

    if (!state.currentSurface->isOnSurface(
            state.options.geoContext, position, direction,
            BoundaryTolerance::Infinite(), state.options.surfaceTolerance)) {
      ACTS_DEBUG(volInfo(state)
                 << "We did not end up on the expected surface. surface = "
                 << state.currentSurface->geometryId()
                 << " position = " << position.transpose()
                 << " direction = " << direction.transpose());

      return Result<void>::failure(NavigatorError::NotOnExpectedSurface);
    }
  }

  return Result<void>::success();
}

NavigationTarget Navigator::nextTarget(State& state, const Vector3& position,
                                       const Vector3& direction) const {
  // Reset the current surface
  state.currentSurface = nullptr;

  if (inactive(state)) {
    return NavigationTarget::None();
  }

  ACTS_VERBOSE(volInfo(state) << "Entering Navigator::nextTarget.");

  NavigationTarget nextTarget = tryGetNextTarget(state, position, direction);
  if (!nextTarget.isNone()) {
    return nextTarget;
  }

  state.reset();
  ++state.statistics.nRenavigations;

  // We might have punched through a boundary and entered another volume
  // so we have to reinitialize
  state.currentVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
      state.options.geoContext, position);

  if (state.currentVolume == nullptr) {
    ACTS_VERBOSE(volInfo(state) << "No volume found, stop navigation.");
    state.navigationBreak = true;
    return NavigationTarget::None();
  }

  if (m_geometryVersion == GeometryVersion::Gen3) {
    // If we re-navigated, we need to recreate the navigation policy state for
    // the new volume if it exists
    const auto* policy = state.currentVolume->navigationPolicy();
    if (policy == nullptr) {
      ACTS_ERROR(volInfo(state) << "No navigation policy for new volume, this "
                                   "should not happen. Stop navigation.");
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "Creating navigation policy state for new "
                                   "volume after renavigation.");
    policy->createState(state.options.geoContext,
                        {.position = position, .direction = direction},
                        state.policyStateManager, logger());
  }

  state.currentLayer =
      state.currentVolume->associatedLayer(state.options.geoContext, position);

  ACTS_VERBOSE(volInfo(state) << "Resolved volume and layer.");

  // Rerun the targeting
  nextTarget = tryGetNextTarget(state, position, direction);
  if (!nextTarget.isNone()) {
    return nextTarget;
  }

  ACTS_VERBOSE(
      volInfo(state)
      << "No targets found again, we got really lost! Stop navigation.");
  state.navigationBreak = true;
  return NavigationTarget::None();
}

bool Navigator::checkTargetValid(State& state, const Vector3& position,
                                 const Vector3& direction) const {
  ACTS_VERBOSE(volInfo(state) << "Entering Navigator::checkTargetValid.");

  if (state.navigationStage == Stage::initial) {
    return false;
  }

  if (state.currentVolume != nullptr &&
      state.currentVolume->navigationPolicy() != nullptr) {
    ACTS_VERBOSE(volInfo(state) << "Checking policy validity for volume");

    auto policyState = state.policyStateManager.currentState();
    bool isValid = state.currentVolume->navigationPolicy()->isValid(
        state.options.geoContext,
        {.position = position, .direction = direction}, policyState, logger());
    return isValid;
  }

  return true;
}

void Navigator::handleSurfaceReached(State& state, const Vector3& position,
                                     const Vector3& direction,
                                     const Surface& surface) const {
  if (inactive(state)) {
    return;
  }

  ACTS_VERBOSE(volInfo(state) << "Entering Navigator::handleSurfaceReached.");

  state.currentSurface = &surface;

  ACTS_VERBOSE(volInfo(state)
               << "Current surface: " << state.currentSurface->geometryId());

  // handling portals in gen3 configuration
  if (m_geometryVersion == GeometryVersion::Gen3) {
    if (state.navCandidate().isPortalTarget() &&
        &state.navCandidate().surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling portal status.");

      // Switch to the next volume using the portal
      const Portal* portal = &state.navCandidate().portal();
      auto res =
          portal->resolveVolume(state.options.geoContext, position, direction);
      if (!res.ok()) {
        ACTS_ERROR(volInfo(state) << "Failed to resolve volume through portal: "
                                  << res.error().message());
        return;
      }

      if (state.currentVolume != nullptr &&
          state.currentVolume->navigationPolicy() != nullptr) {
        ACTS_VERBOSE(volInfo(state)
                     << "Popping navigation policy state for volume on exit");
        state.currentVolume->navigationPolicy()->popState(
            state.policyStateManager, logger());
      }

      state.currentVolume = res.value();

      // partial reset
      state.resetAfterVolumeSwitch();

      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Volume updated.");

        const auto* policy = state.currentVolume->navigationPolicy();

        if (policy == nullptr) {
          ACTS_ERROR(
              volInfo(state)
              << "No navigation policy for new volume, this should not happen");
          return;
        }

        ACTS_VERBOSE(volInfo(state)
                     << "Creating navigation policy state for new "
                        "volume after portal transition.");
        policy->createState(state.options.geoContext,
                            {.position = position, .direction = direction},
                            state.policyStateManager, logger());

        // this is set only for the check target validity since gen3 does not
        // care
        state.navigationStage = Stage::surfaceTarget;
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "No more volume to progress to, stopping navigation.");
        state.navigationBreak = true;
      }
    }
    // Mark reached free candidates
    else if (&state.navCandidate().surface() == &surface &&
             surface.geometryId() == GeometryIdentifier{}) {
      auto freeItr = std::ranges::find_if(
          state.freeCandidates,
          [&surface](const std::pair<const Surface*, bool>& cand) {
            return &surface == cand.first;
          });
      if (freeItr != state.freeCandidates.end()) {
        freeItr->second = true;
      }
    }
    return;
  }

  if (state.navigationStage == Stage::surfaceTarget &&
      &state.navSurface().surface() == &surface) {
    ACTS_VERBOSE(volInfo(state) << "Handling surface status.");

    return;
  }

  if (state.navigationStage == Stage::layerTarget &&
      &state.navLayer().surface() == &surface) {
    ACTS_VERBOSE(volInfo(state) << "Handling layer status.");

    // Switch to the next layer
    state.currentLayer = &state.navLayer().layer();
    state.navigationStage = Stage::surfaceTarget;

    // partial reset
    state.resetAfterLayerSwitch();

    return;
  }

  if (state.navigationStage == Stage::boundaryTarget &&
      &state.navBoundary().surface() == &surface) {
    ACTS_VERBOSE(volInfo(state) << "Handling boundary status.");

    // Switch to the next volume using the boundary
    const BoundarySurface* boundary = &state.navBoundary().boundarySurface();
    assert(boundary != nullptr && "Retrieved boundary surface is nullptr");
    state.currentVolume =
        boundary->attachedVolume(state.options.geoContext, position, direction);

    // partial reset
    state.resetAfterVolumeSwitch();

    if (state.currentVolume != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Volume updated.");
      state.navigationStage = Stage::layerTarget;
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "No more volume to progress to, stopping navigation.");
      state.navigationBreak = true;
    }

    return;
  }

  ACTS_ERROR(volInfo(state) << "Surface reached but unknown state.");
}

NavigationTarget Navigator::getNextTargetGen1(State& state,
                                              const Vector3& position,
                                              const Vector3& direction) const {
  // Try targeting the surfaces - then layers - then boundaries
  if (state.navigationStage == Stage::surfaceTarget) {
    if (!state.navSurfaceIndex.has_value()) {
      // First time, resolve the surfaces
      resolveSurfaces(state, position, direction);
      state.navSurfaceIndex = 0;
    } else {
      ++state.navSurfaceIndex.value();
    }
    if (state.navSurfaceIndex.value() < state.navSurfaces.size()) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
      return state.navSurface();
    } else {
      // This was the last surface, switch to layers
      ACTS_VERBOSE(volInfo(state) << "Target layers.");
      state.navigationStage = Stage::layerTarget;
    }
  }

  if (state.navigationStage == Stage::layerTarget) {
    if (!state.navLayerIndex.has_value()) {
      // First time, resolve the layers
      resolveLayers(state, position, direction);
      state.navLayerIndex = 0;
    } else {
      ++state.navLayerIndex.value();
    }
    if (state.navLayerIndex.value() < state.navLayers.size()) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
      return state.navLayer();
    } else {
      // This was the last layer, switch to boundaries
      ACTS_VERBOSE(volInfo(state) << "Target boundaries.");
      state.navigationStage = Stage::boundaryTarget;
    }
  }

  if (state.navigationStage == Stage::boundaryTarget) {
    if (!state.navBoundaryIndex.has_value()) {
      // First time, resolve the boundaries
      resolveBoundaries(state, position, direction);
      state.navBoundaryIndex = 0;
    } else {
      ++state.navBoundaryIndex.value();
    }
    if (state.navBoundaryIndex.value() < state.navBoundaries.size()) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
      return state.navBoundary();
    } else {
      // This was the last boundary, we have to leave the volume somehow,
      // renavigate
      ACTS_VERBOSE(volInfo(state) << "Boundary targets exhausted. Renavigate.");
      return NavigationTarget::None();
    }
  }

  ACTS_VERBOSE(volInfo(state) << "Unknown state. No target found. Renavigate.");
  return NavigationTarget::None();
}

NavigationTarget Navigator::getNextTargetGen3(State& state,
                                              const Vector3& position,
                                              const Vector3& direction) const {
  if (state.currentVolume == nullptr) {
    ACTS_VERBOSE(volInfo(state) << "No volume to get next target.");
    return NavigationTarget::None();
  }

  auto policyState = state.policyStateManager.currentState();
  bool isValid = state.currentVolume->navigationPolicy()->isValid(
      state.options.geoContext, {.position = position, .direction = direction},
      policyState, logger());

  ACTS_VERBOSE(volInfo(state) << "Current policy says navigation sequence is "
                              << (isValid ? "VALID" : "INVALID"));

  ACTS_VERBOSE(volInfo(state)
               << "Current candidate index is "
               << (state.navCandidateIndex.has_value()
                       ? std::to_string(state.navCandidateIndex.value())
                       : "n/a"));

  if (!isValid || !state.navCandidateIndex.has_value()) {
    // first time, resolve the candidates
    resolveCandidates(state, position, direction);
    state.navCandidateIndex = 0;
  } else {
    ++state.navCandidateIndex.value();
  }
  if (state.navCandidateIndex.value() < state.navCandidates.size()) {
    ACTS_VERBOSE(volInfo(state)
                 << "Target set to next candidate " << state.navCandidate());
    return state.navCandidate();
  } else {
    ACTS_VERBOSE(volInfo(state) << "Candidate targets exhausted. Renavigate.");
    return NavigationTarget::None();
  }
}

NavigationTarget Navigator::tryGetNextTarget(State& state,
                                             const Vector3& position,
                                             const Vector3& direction) const {
  // Try different approach to get navigation target for gen1 and gen3
  // configuration

  // This is common, in gen1 we start by surfaces and in gen3 we always look
  // for surfaces
  if (state.navigationStage == Stage::initial) {
    ACTS_VERBOSE(volInfo(state) << "Target surfaces.");
    state.navigationStage = Stage::surfaceTarget;
  }

  if (m_geometryVersion == GeometryVersion::Gen1) {
    return getNextTargetGen1(state, position, direction);

  } else {  // gen3 handling of the next target

    return getNextTargetGen3(state, position, direction);
  }
}

void Navigator::resolveCandidates(State& state, const Vector3& position,
                                  const Vector3& direction) const {
  if (state.currentVolume == nullptr) {
    ACTS_VERBOSE(volInfo(state) << "No volume to resolve candidates.");
    return;
  }
  ACTS_VERBOSE(volInfo(state) << "Searching for compatible candidates.");

  state.stream.reset();
  AppendOnlyNavigationStream appendOnly{state.stream};
  NavigationArguments args;
  args.position = position;
  args.direction = direction;

  const INavigationPolicy* policy = state.currentVolume->navigationPolicy();
  if (policy == nullptr) {
    ACTS_ERROR(volInfo(state) << "No navigation policy found for volume. "
                                 "Cannot resolve navigation candidates.");
    throw std::runtime_error(
        "Navigator: No navigation policy found for current volume.");
  }

  auto policyState = state.policyStateManager.currentState();
  state.currentVolume->initializeNavigationCandidates(
      state.options.geoContext, args, policyState, appendOnly, logger());

  ACTS_VERBOSE(volInfo(state) << "Found " << state.stream.candidates().size()
                              << " navigation candidates.");
  if (!state.options.externalSurfaces.empty()) {
    for (const GeometryIdentifier& geoId : state.options.externalSurfaces) {
      // Don't add any surface which is not in the same volume (volume bits)
      // or sub volume (extra bits)
      if (geoId.volume() != state.currentVolume->geometryId().volume() ||
          geoId.extra() != state.currentVolume->geometryId().extra()) {
        continue;
      }
      const Surface* surface = m_cfg.trackingGeometry->findSurface(geoId);
      assert(surface != nullptr);
      ACTS_VERBOSE(volInfo(state) << "Try to navigate to " << surface->type()
                                  << " surface " << geoId);
      appendOnly.addSurfaceCandidate(*surface, BoundaryTolerance::Infinite());
    };
  }
  bool pruneFreeCand{false};
  if (!state.freeCandidates.empty()) {
    for (const auto& [surface, wasReached] : state.freeCandidates) {
      /// Don't process already reached surfaces again
      if (wasReached) {
        continue;
      }
      if (!state.options.freeSurfaceSelector.connected() ||
          state.options.freeSurfaceSelector(state.options.geoContext,
                                            *state.currentVolume, position,
                                            direction, *surface)) {
        ACTS_VERBOSE(volInfo(state)
                     << "Append free " << surface->type() << " surface  \n"
                     << surface->toStream(state.options.geoContext));
        appendOnly.addSurfaceCandidate(*surface, BoundaryTolerance::Infinite());
        pruneFreeCand = !state.options.freeSurfaceSelector.connected();
      }
    };
  }
  state.stream.initialize(state.options.geoContext, {position, direction},
                          BoundaryTolerance::None(),
                          state.options.surfaceTolerance);

  ACTS_VERBOSE(volInfo(state) << "Now " << state.stream.candidates().size()
                              << " navigation candidates after initialization");

  state.navCandidates.clear();

  double farLimit = state.options.farLimit;
  // If the user has not provided the selection delegate, then
  // just apply a simple candidate pruning. Constrain the maximum
  // reach of the navigation to the last portal in the state
  if (pruneFreeCand) {
    farLimit = state.options.nearLimit;
    for (const auto& candidate : state.stream.candidates()) {
      if (candidate.isPortalTarget()) {
        farLimit = std::max(farLimit, candidate.intersection().pathLength() +
                                          state.options.surfaceTolerance);
      }
    }
  }

  for (auto& candidate : state.stream.candidates()) {
    if (!detail::checkPathLength(candidate.intersection().pathLength(),
                                 state.options.nearLimit, farLimit, logger())) {
      continue;
    }

    state.navCandidates.emplace_back(candidate);
  }

  // Sort the candidates with the path length
  std::ranges::sort(state.navCandidates, [](const auto& a, const auto& b) {
    return a.intersection().pathLength() < b.intersection().pathLength();
  });

  // Print the navigation candidates

  if (logger().doPrint(Logging::VERBOSE)) {
    std::ostringstream os;
    os << "Navigation candidates: " << state.navCandidates.size() << "\n";
    for (auto& candidate : state.navCandidates) {
      os << " -- " << candidate << " at "
         << candidate.pathLength() / UnitConstants::mm << "mm\n";
    }

    logger().log(Logging::VERBOSE, os.str());
  }
}

void Navigator::resolveSurfaces(State& state, const Vector3& position,
                                const Vector3& direction) const {
  // Gen-1 surface resolution
  ACTS_VERBOSE(volInfo(state) << "Searching for compatible surfaces.");

  const Layer* currentLayer = state.currentLayer;

  if (currentLayer == nullptr) {
    ACTS_VERBOSE(volInfo(state) << "No layer to resolve surfaces.");
    return;
  }

  const Surface* layerSurface = &currentLayer->surfaceRepresentation();

  NavigationOptions<Surface> navOpts;
  navOpts.resolveSensitive = m_cfg.resolveSensitive;
  navOpts.resolveMaterial = m_cfg.resolveMaterial;
  navOpts.resolvePassive = m_cfg.resolvePassive;
  navOpts.startObject = state.currentSurface;
  navOpts.endObject = state.targetSurface;
  navOpts.nearLimit = state.options.nearLimit;
  navOpts.farLimit = state.options.farLimit;

  if (!state.options.externalSurfaces.empty()) {
    const auto layerId = layerSurface->geometryId().layer();
    for (const GeometryIdentifier& id : state.options.externalSurfaces) {
      if (id.layer() == layerId) {
        navOpts.externalSurfaces.push_back(id);
      }
    }
  }

  // Request the compatible surfaces
  state.navSurfaces = currentLayer->compatibleSurfaces(
      state.options.geoContext, position, direction, navOpts);
  // Sort the surfaces by path length.
  // Special care is taken for the external surfaces which should always
  // come first, so they are preferred to be targeted and hit first.
  std::ranges::sort(state.navSurfaces, [&state](const NavigationTarget& a,
                                                const NavigationTarget& b) {
    // Prefer to sort by path length. We assume surfaces are at the same
    // distance if the difference is smaller than the tolerance.
    if (std::abs(a.pathLength() - b.pathLength()) >
        state.options.surfaceTolerance) {
      return NavigationTarget::pathLengthOrder(a, b);
    }
    // If the path length is practically the same, sort by geometry.
    // First we check if one of the surfaces is external.
    bool aIsExternal = a.boundaryTolerance().isInfinite();
    bool bIsExternal = b.boundaryTolerance().isInfinite();
    if (aIsExternal == bIsExternal) {
      // If both are external or both are not external, sort by geometry
      // identifier
      return a.surface().geometryId() < b.surface().geometryId();
    }
    // If only one is external, it should come first
    return aIsExternal;
  });
  // For now we implicitly remove overlapping surfaces.
  // For track finding it might be useful to discover overlapping surfaces
  // and check for compatible measurements. This is under investigation
  // and might be implemented in the future.
  auto toBeRemoved =
      std::ranges::unique(state.navSurfaces, [&](const auto& a, const auto& b) {
        return std::abs(a.pathLength() - b.pathLength()) <
               state.options.surfaceTolerance;
      });
  if (toBeRemoved.begin() != toBeRemoved.end()) {
    ACTS_VERBOSE(volInfo(state)
                 << "Removing "
                 << std::distance(toBeRemoved.begin(), toBeRemoved.end())
                 << " overlapping surfaces.");
  }
  state.navSurfaces.erase(toBeRemoved.begin(), toBeRemoved.end());

  // Print surface information
  if (logger().doPrint(Logging::VERBOSE)) {
    std::ostringstream os;
    os << state.navSurfaces.size();
    os << " surface candidates found at path(s): ";
    for (auto& sfc : state.navSurfaces) {
      os << sfc.pathLength() << "  ";
    }
    logger().log(Logging::VERBOSE, os.str());
  }

  if (state.navSurfaces.empty()) {
    ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
  }
}

void Navigator::resolveLayers(State& state, const Vector3& position,
                              const Vector3& direction) const {
  ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

  NavigationOptions<Layer> navOpts;
  navOpts.resolveSensitive = m_cfg.resolveSensitive;
  navOpts.resolveMaterial = m_cfg.resolveMaterial;
  navOpts.resolvePassive = m_cfg.resolvePassive;
  navOpts.startObject = state.currentLayer;
  navOpts.nearLimit = state.options.nearLimit;
  navOpts.farLimit = state.options.farLimit;

  // Request the compatible layers
  state.navLayers = state.currentVolume->compatibleLayers(
      state.options.geoContext, position, direction, navOpts);
  std::ranges::sort(state.navLayers, NavigationTarget::pathLengthOrder);

  // Print layer information
  if (logger().doPrint(Logging::VERBOSE)) {
    std::ostringstream os;
    os << state.navLayers.size();
    os << " layer candidates found at path(s): ";
    for (auto& lc : state.navLayers) {
      os << lc.pathLength() << "  ";
    }
    logger().log(Logging::VERBOSE, os.str());
  }

  if (state.navLayers.empty()) {
    ACTS_VERBOSE(volInfo(state) << "No layer candidates found.");
  }
}

void Navigator::resolveBoundaries(State& state, const Vector3& position,
                                  const Vector3& direction) const {
  ACTS_VERBOSE(volInfo(state) << "Searching for compatible boundaries.");

  NavigationOptions<Surface> navOpts;
  navOpts.startObject = state.currentSurface;
  navOpts.nearLimit = state.options.nearLimit;
  navOpts.farLimit = state.options.farLimit;

  ACTS_VERBOSE(volInfo(state)
               << "Try to find boundaries, we are at: " << toString(position)
               << ", dir: " << toString(direction));

  // Request the compatible boundaries
  state.navBoundaries = state.currentVolume->compatibleBoundaries(
      state.options.geoContext, position, direction, navOpts, logger());
  std::ranges::sort(state.navBoundaries, NavigationTarget::pathLengthOrder);

  // Print boundary information
  if (logger().doPrint(Logging::VERBOSE)) {
    std::ostringstream os;
    os << state.navBoundaries.size();
    os << " boundary candidates found at path(s): ";
    for (const auto& bc : state.navBoundaries) {
      os << bc.pathLength() << "  ";
    }
    logger().log(Logging::VERBOSE, os.str());
  }

  if (state.navBoundaries.empty()) {
    ACTS_VERBOSE(volInfo(state) << "No boundary candidates found.");
  }
}

bool Navigator::inactive(const State& state) const {
  // Turn the navigator into void when you are instructed to do nothing
  if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
      !m_cfg.resolvePassive) {
    return true;
  }

  if (state.navigationBreak) {
    return true;
  }

  return false;
}

std::string Navigator::volInfo(const State& state) const {
  if (state.currentVolume == nullptr) {
    return "No Volume | ";
  }
  std::stringstream sstr{};
  sstr << state.currentVolume->volumeName() << " ("
       << state.currentVolume->geometryId() << ") | ";
  return sstr.str();
}

}  // namespace Acts
