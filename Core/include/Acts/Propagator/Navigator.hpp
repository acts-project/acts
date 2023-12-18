// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/NavigationHelpers.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <cstddef>
#include <limits>
#include <sstream>
#include <string>

namespace Acts {

/// @brief struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryCheck boundaryCheck = BoundaryCheck(true);

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// object to check against: at start
  const object_t* startObject = nullptr;
  /// object to check against: at end
  const object_t* endObject = nullptr;

  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The minimum distance for a surface to be considered
  double nearLimit = 0;
  /// The maximum distance for a surface to be considered
  double farLimit = std::numeric_limits<double>::max();
};

/// @brief Steers the propagation through the geometry by adjusting the step
///        size and providing the next surface to be targeted.
///
/// The Navigator is part of the propagation and responsible for steering
/// the step size in order to encounter all the relevant surfaces which are
/// intersected by the trajectory.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation flow is hit, it is
/// set to the navigation state, such that other actors can deal with it.
///
/// The current target surface is referenced by an index which points into
/// the navigation candidates. The navigation candidates are ordered by the
/// path length to the surface. If a surface is hit, the
/// `state.navigation.currentSurface` pointer is set. This actors to observe
/// that we are on a surface.
///
class Navigator {
 public:
  /// @brief Configuration for this Navigator
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    /// @note that this includes representative surfaces which might
    ///       overlap with other surfaces
    bool resolvePassive = false;

    /// Which boundary checks to perform for surface approach
    BoundaryCheck boundaryCheckSurfaceApproach = BoundaryCheck(true);
    /// Whether to perform boundary checks for layer resolving (improves
    /// navigation for bended tracks)
    BoundaryCheck boundaryCheckLayerResolving = BoundaryCheck(true);

    /// Whether to reinitialize the navigation candidates on layer hit.
    /// This improves navigation for bended tracks but is more expensive.
    bool reinitializeOnLayerHit = false;
  };

  using IntersectionCandidates = std::vector<detail::IntersectionCandidate>;

  using ExternalSurfaces = std::multimap<uint64_t, GeometryIdentifier>;

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    /// The vector of intersection candidates to work through
    IntersectionCandidates candidates = {};
    /// The current candidate index of the navigation state
    std::size_t candidateIndex = 0;

    /// Provides easy access to the active navigation candidate
    const detail::IntersectionCandidate& candidate() const {
      return candidates.at(candidateIndex);
    }

    /// Externally provided surfaces which are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    // Starting geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the starting surface (e.g. MaterialInteraction).
    const TrackingVolume* startVolume = nullptr;
    const Layer* startLayer = nullptr;
    const Surface* startSurface = nullptr;

    // Target geometry information of the navigation which should only be set
    // while initialization. NOTE: This information is mostly used by actors to
    // check if we are on the target surface (e.g. MaterialInteraction).
    const Surface* targetSurface = nullptr;

    // Current geometry information of the navigation which is set during
    // initialization and potentially updated after each step.
    const Surface* currentSurface = nullptr;
    const Layer* currentLayer = nullptr;
    const TrackingVolume* currentVolume = nullptr;

    /// Indicator if the target is reached
    bool targetReached = false;
    /// If a break has been detected
    bool navigationBreak = false;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {}

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void currentSurface(State& state, const Surface* surface) const {
    state.currentSurface = surface;
  }

  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  void insertExternalSurface(State& state, GeometryIdentifier geoid) const {
    state.externalSurfaces.insert(
        std::pair<uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
  }

  /// @brief Initialize call - start of navigation
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE("Initialization.");

    // Depending on the start information set there are different modes for
    // finding the start layer and volume.
    if (state.navigation.startSurface != nullptr &&
        state.navigation.startSurface->associatedLayer() != nullptr) {
      // Everything was given by the user
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      state.navigation.startLayer =
          state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume =
          state.navigation.startLayer->trackingVolume();
    } else if (state.navigation.startVolume != nullptr) {
      // Only the layer has to be determined
      ACTS_VERBOSE(
          "Fast start initialization through association from Volume.");
      state.navigation.startLayer =
          state.navigation.startVolume->associatedLayer(
              state.geoContext, stepper.position(state.stepping));
    } else {
      // Otherwise we only rely on the position and direction
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
      state.navigation.startVolume =
          m_cfg.trackingGeometry->lowestTrackingVolume(
              state.geoContext, stepper.position(state.stepping));
    }

    // Special handling if the start surface is a boundary surface since the
    // start volume might not be correct.
    if (state.navigation.startSurface != nullptr) {
      for (const auto& boundary :
           state.navigation.startVolume->boundarySurfaces()) {
        if (state.navigation.startSurface ==
            &boundary->surfaceRepresentation()) {
          state.navigation.startVolume = boundary->attachedVolume(
              state.geoContext, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping));
          ACTS_VERBOSE("We are starting from a boundary surface "
                       << state.navigation.startSurface->geometryId());
          break;
        }
      }
    }

    // Initialize current volume, layer and surface
    {
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      } else {
        ACTS_VERBOSE("Start volume not resolved.");
      }

      state.navigation.currentLayer = state.navigation.startLayer;
      if (state.navigation.currentLayer != nullptr) {
        ACTS_VERBOSE(volInfo(state)
                     << "Start layer resolved to "
                     << state.navigation.currentLayer->geometryId() << " .");
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start layer set.");
      }

      state.navigation.currentSurface = state.navigation.startSurface;
      if (state.navigation.currentSurface != nullptr) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current surface set to start surface "
                     << state.navigation.currentSurface->geometryId());
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    // Initialize navigation candidates for the start volume
    reinitializeCandidates(state, stepper);
  }

  /// @brief Navigator pre step call
  ///
  /// This determines the next surface to be targeted and sets the step length
  /// accordingly.
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::preStep.");

    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    bool reinitialized = false;
    // Check next navigation candidate
    while (state.navigation.candidateIndex !=
           state.navigation.candidates.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << (state.navigation.candidates.size() -
                       state.navigation.candidateIndex)
                   << " out of " << state.navigation.candidates.size()
                   << " surfaces remain to try.");

      const auto& candidate = state.navigation.candidate();
      const auto& intersection = candidate.intersection;
      const Surface& surface = *intersection.object();
      BoundaryCheck boundaryCheck = candidate.boundaryCheck;

      ACTS_VERBOSE(volInfo(state) << "Next surface candidate will be "
                                  << surface.geometryId());

      // Determine the surface status and set the step length accordingly
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, surface, intersection.index(),
          state.options.direction, boundaryCheck,
          state.options.surfaceTolerance, logger());

      // We should never be on surface before we actually targeted a surface
      if (surfaceStatus == IntersectionStatus::onSurface) {
        ACTS_ERROR(volInfo(state)
                   << "We are on surface " << surface.geometryId()
                   << " before trying to reach it. This should not happen. "
                      "Good luck.");
        ++state.navigation.candidateIndex;
        continue;
      }

      if (surfaceStatus == IntersectionStatus::reachable) {
        ACTS_VERBOSE(volInfo(state)
                     << "Surface reachable, step size updated to "
                     << stepper.outputStepSize(state.stepping));
        break;
      }

      // Handle an unreachable candidates. This can happen for sensitive
      // surfaces as we may not encounter them. If this happens for layers or
      // boundaries we have to reinitialize the navigation.
      if (candidate.template checkType<Surface>()) {
        // Skip if this is an ordinary surface
        ACTS_VERBOSE(volInfo(state) << "Surface unreachable, skip.");
        ++state.navigation.candidateIndex;

        assert(state.navigation.candidateIndex !=
                   state.navigation.candidates.size() &&
               "No more candidates.");
      } else {
        // Renavigate as this is a sign of an invalidated navigation stream
        ACTS_VERBOSE(volInfo(state)
                     << "Layer/Boundary unreachable, renavigate.");
        // If we tried to reinitialize already we break the look and error
        if (reinitialized) {
          ACTS_ERROR(volInfo(state) << "Renavigation failed. Good luck.");
          // Set navigation break and release the navigation step size
          state.navigation.navigationBreak = true;
          stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);
          break;
        }
        reinitializeCandidates(state, stepper);
        reinitialized = true;
      }
    }

    // There should always be an active candidate we are trying to navigate to.
    // Otherwise we are in trouble.
    // Final boundary or target surface not found.
    if (state.navigation.candidateIndex == state.navigation.candidates.size()) {
      ACTS_ERROR(volInfo(state) << "Exhausted navigation candidates.");
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping, ConstrainedStep::actor);
    }

    // Navigator preStep always resets the current surface
    state.navigation.currentSurface = nullptr;
  }

  /// @brief Navigator post step call
  ///
  /// This determines if we hit the next navigation candidate and deals with it
  /// accordingly. It sets the current surface, enters layers and changes
  /// volumes.
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return Boolean to indicate if we continue with the actors and
  ///         aborters or if we should target again.
  template <typename propagator_state_t, typename stepper_t>
  bool postStep(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::postStep.");

    // Check if the navigator is inactive
    if (inactive(state)) {
      return true;
    }

    assert(state.navigation.currentSurface == nullptr &&
           "Current surface must be reset.");

    const auto& candidate = state.navigation.candidate();
    const auto& intersection = candidate.intersection;
    const Surface& surface = *intersection.object();

    // Determine the surface status without boundary check.
    // TODO This will also update the step length which is unnecessary
    Intersection3D::Status surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, surface, intersection.index(), state.options.direction,
        BoundaryCheck(true), state.options.surfaceTolerance, logger());

    // Check if we are on surface otherwise wait for the next step.
    if (surfaceStatus != IntersectionStatus::onSurface) {
      ACTS_VERBOSE(volInfo(state) << "Staying focussed on surface.");
      return true;
    }

    // We are on surface and switch already to the next candidate which is
    // necessary for the next `preStep` call.
    ++state.navigation.candidateIndex;

    /*
    // Determine the surface status with boundary check.
    // TODO second intersection should not be necessary
    surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, surface, intersection.index(), state.options.direction,
        BoundaryCheck(true), state.options.surfaceTolerance, logger());

    // Check if we are still on surface, otherwise we are out of bounds.
    if (surfaceStatus != IntersectionStatus::onSurface) {
      ACTS_VERBOSE(volInfo(state)
                   << "Surface successfully hit, but outside bounds.");
      return true;
    }
    */

    ACTS_VERBOSE(volInfo(state) << "Surface successfully hit, storing it.");
    // Update navigation state, so actors and aborters can access it
    state.navigation.currentSurface = &surface;

    if (state.navigation.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state)
                   << "Current surface set to " << surface.geometryId());
    }

    // Depending on what kind of surface we intersected we need to update the
    // navigation state further.
    if (candidate.template checkType<Surface>()) {
      // In case of a surface we are done already

      ACTS_VERBOSE(volInfo(state) << "This is a surface");

      assert(state.navigation.candidateIndex !=
                 state.navigation.candidates.size() &&
             "No more candidates.");
    } else if (candidate.template checkType<Layer>()) {
      // In case of a layer we need to update our navigation candidates with the
      // surfaces we might hit in the layer

      ACTS_VERBOSE(volInfo(state)
                   << "This is a layer. Initialize layer candidates");

      state.navigation.currentLayer = candidate.template object<Layer>();

      if (m_cfg.reinitializeOnLayerHit) {
        reinitializeCandidates(state, stepper);
      } else {
        // Note that this is hacky as we do not reintersect the existing
        // candidates.

        initializeLayerSurfaceCandidates(state, stepper);

        std::sort(state.navigation.candidates.begin() +
                      state.navigation.candidateIndex,
                  state.navigation.candidates.end(),
                  detail::IntersectionCandidate::forwardOrder);
      }
    } else if (candidate.template checkType<BoundarySurface>()) {
      // In case of a boundary we need to switch volume and reinitialize

      ACTS_VERBOSE(volInfo(state)
                   << "This is a boundary. Reinitialize navigation");

      const auto* boundary = candidate.template object<BoundarySurface>();

      state.navigation.currentLayer = nullptr;
      state.navigation.currentVolume = boundary->attachedVolume(
          state.geoContext, stepper.position(state.stepping),
          state.options.direction * stepper.direction(state.stepping));

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      reinitializeCandidates(state, stepper);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }

    return true;
  }

 private:
  /// Helper method to initialize navigation candidates for the current volume
  /// boundaries.
  template <typename propagator_state_t, typename stepper_t>
  void initializeVolumeBoundaryCandidates(propagator_state_t& state,
                                          const stepper_t& stepper) const {
    const TrackingVolume* volume = state.navigation.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume boundaries");

    if (volume == nullptr) {
      state.navigation.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    // The navigation options
    NavigationOptions<Surface> navOpts;
    navOpts.boundaryCheck = BoundaryCheck(true);
    // Exclude the current surface in case it's a boundary
    navOpts.startObject = state.navigation.currentSurface;
    navOpts.nearLimit = state.options.surfaceTolerance;
    navOpts.farLimit = std::numeric_limits<double>::max();

    ACTS_VERBOSE(volInfo(state)
                 << "Try to find boundaries, we are at: "
                 << stepper.position(state.stepping).transpose()
                 << ", dir: " << stepper.direction(state.stepping).transpose());

    auto boundaries = volume->compatibleBoundaries(
        state.geoContext, stepper.position(state.stepping),
        state.options.direction * stepper.direction(state.stepping), navOpts,
        logger());

    // Screen output where they are
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream oss;
      oss << boundaries.size();
      oss << " boundary candidates found at path(s): ";
      for (const auto& [boundary, _] : boundaries) {
        oss << boundary.pathLength() << "  ";
      }
      ACTS_VERBOSE(oss.str());
    }

    for (const auto& [boundary, object] : boundaries) {
      state.navigation.candidates.emplace_back(boundary, object,
                                               BoundaryCheck(true));
    }
  }

  /// Helper method to initialize navigation candidates for the current volume
  /// layers.
  template <typename propagator_state_t, typename stepper_t>
  void initializeVolumeLayerCandidates(propagator_state_t& state,
                                       const stepper_t& stepper) const {
    const TrackingVolume* volume = state.navigation.currentVolume;
    ACTS_VERBOSE(volInfo(state) << "Initialize volume layers");

    if (volume == nullptr) {
      state.navigation.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    // Create the navigation options
    // - and get the compatible layers, start layer will be excluded
    NavigationOptions<Layer> navOpts;
    navOpts.boundaryCheck = m_cfg.boundaryCheckLayerResolving;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.navigation.currentLayer;
    navOpts.nearLimit = state.options.surfaceTolerance;
    navOpts.farLimit = std::numeric_limits<double>::max();

    auto layers = volume->compatibleLayers(
        state.geoContext, stepper.position(state.stepping),
        state.options.direction * stepper.direction(state.stepping), navOpts);

    // Screen output where they are
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream oss;
      oss << layers.size();
      oss << " layer candidates found at path(s): ";
      for (const auto& [layer, _] : layers) {
        oss << layer.pathLength() << "  ";
      }
      ACTS_VERBOSE(oss.str());
    }

    for (const auto& [layer, object] : layers) {
      state.navigation.candidates.emplace_back(layer, object,
                                               BoundaryCheck(true));
    }
  }

  /// Helper method to initialize navigation candidates for the current layer.
  template <typename propagator_state_t, typename stepper_t>
  void initializeLayerSurfaceCandidates(propagator_state_t& state,
                                        const stepper_t& stepper) const {
    const Layer* layer = state.navigation.currentLayer;
    assert(layer != nullptr && "Current layer must be set.");

    ACTS_VERBOSE(volInfo(state) << "Initialize layer surface candidates for "
                                << layer->geometryId());

    // Use navigation parameters and NavigationOptions
    NavigationOptions<Surface> navOpts;
    navOpts.boundaryCheck = BoundaryCheck(true);
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.navigation.currentSurface;
    navOpts.endObject = state.navigation.targetSurface;
    navOpts.nearLimit = state.options.surfaceTolerance;
    navOpts.farLimit = std::numeric_limits<double>::max();

    if (!state.navigation.externalSurfaces.empty()) {
      auto layerID = layer->geometryId().layer();
      auto externalSurfaceRange =
          state.navigation.externalSurfaces.equal_range(layerID);
      navOpts.externalSurfaces.reserve(
          state.navigation.externalSurfaces.count(layerID));
      for (auto itSurface = externalSurfaceRange.first;
           itSurface != externalSurfaceRange.second; itSurface++) {
        navOpts.externalSurfaces.push_back(itSurface->second);
      }
    }

    auto surfaces = layer->compatibleSurfaces(
        state.geoContext, stepper.position(state.stepping),
        state.options.direction * stepper.direction(state.stepping), navOpts);

    // Screen output where they are
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream oss;
      oss << surfaces.size();
      oss << " surface candidates found at path(s): ";
      for (const auto& surface : surfaces) {
        oss << surface.pathLength() << "  ";
      }
      ACTS_VERBOSE(oss.str());
    }

    for (const auto& surface : surfaces) {
      state.navigation.candidates.emplace_back(
          surface, surface.object(), m_cfg.boundaryCheckSurfaceApproach);
    }
  }

  /// Helper method to reset and reinitialize the navigation candidates.
  template <typename propagator_state_t, typename stepper_t>
  void reinitializeCandidates(propagator_state_t& state,
                              const stepper_t& stepper) const {
    state.navigation.candidates.clear();
    state.navigation.candidateIndex = 0;

    if (state.navigation.currentLayer != nullptr) {
      initializeLayerSurfaceCandidates(state, stepper);
    }

    initializeVolumeLayerCandidates(state, stepper);
    initializeVolumeBoundaryCandidates(state, stepper);

    std::sort(state.navigation.candidates.begin(),
              state.navigation.candidates.end(),
              detail::IntersectionCandidate::forwardOrder);
  }

  /// Checks if a navigation break had been triggered or navigator is
  /// misconfigured.
  ///
  /// @tparam propagator_state_t The state type of the propagator
  ///
  /// @param [in,out] state is the propagation state object
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t>
  bool inactive(propagator_state_t& state) const {
    // Void behavior in case no tracking geometry is present
    if (!m_cfg.trackingGeometry) {
      return true;
    }
    // turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }
    // handle navigation break
    if (state.navigation.navigationBreak) {
      return true;
    }
    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume != nullptr
                ? state.navigation.currentVolume->volumeName()
                : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
