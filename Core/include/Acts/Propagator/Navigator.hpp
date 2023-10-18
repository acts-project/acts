// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryCheck boundaryCheck = true;

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

  /// Target surface to exclude
  const Surface* targetSurface = nullptr;
  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The maximum path limit for this navigation step
  double pathLimit = std::numeric_limits<double>::max();

  /// The overstep tolerance for this navigation step
  double overstepLimit = 0;

  /// Force intersection with boundaries
  bool forceIntersectBoundaries = false;
};

/// Navigator class
///
/// This is an Actor to be added to the ActorList in order to navigate
/// through the static tracking geometry setup.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation  flow is hit, it is
/// set to the propagation satate, such that other actors can deal with it.
/// This navigation actor thus always needs to run first!
/// It does two things: it figures out the order of volumes, layers and
/// surfaces.
///
/// The current target surface is the surface pointed to by of the index
/// for the surfaces, layers or volume boundaries.
/// If a surface is found, the state.navigation.currentSurface
/// pointer is set. This  enables subsequent actors to react. Secondly, this
/// actor uses the ordered indices to figure out which surface, layer or
/// volume boundary is _supposed_ to be hit next. It then sets the maximum
/// step size to the path length found out by straight line intersection. If
/// the state is not on surface, it also  re-computes the step size, to make
/// sure we end up at the desired surface.
///
class Navigator {
 public:
  using Surfaces = std::vector<const Surface*>;

  using NavigationSurfaces =
      boost::container::small_vector<SurfaceIntersection, 10>;

  using NavigationLayers =
      boost::container::small_vector<LayerIntersection, 10>;

  using NavigationBoundaries =
      boost::container::small_vector<BoundaryIntersection, 4>;

  using ExternalSurfaces = std::multimap<uint64_t, GeometryIdentifier>;

  /// The navigation stage
  enum struct Stage : int {
    undefined = 0,
    surfaceTarget = 1,
    layerTarget = 2,
    boundaryTarget = 3
  };

  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// Configuration for this Navigator
    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;

    /// Whether to perform boundary checks for layer resolving (improves
    /// navigation for bended tracks)
    BoundaryCheck boundaryCheckLayerResolving = true;
  };

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface index of the navigation state
    std::size_t navSurfaceIndex = navSurfaces.size();

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer index of the navigation state
    std::size_t navLayerIndex = navLayers.size();

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary index of the navigation state
    std::size_t navBoundaryIndex = navBoundaries.size();

    auto navSurface() const { return navSurfaces.at(navSurfaceIndex); }
    auto navLayer() const { return navLayers.at(navLayerIndex); }
    auto navBoundary() const { return navBoundaries.at(navBoundaryIndex); }

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    /// Navigation state: the world volume
    const TrackingVolume* worldVolume = nullptr;

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;
    /// Navigation state: the target layer
    const Layer* targetLayer = nullptr;
    /// Navigation state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator for start layer treatment
    bool startLayerResolved = false;
    /// Indicator if the target is reached
    bool targetReached = false;
    /// Indicator that the last VolumeHierarchy surface was reached
    /// skip the next layer targeting to the next boundary/volume
    bool lastHierarchySurfaceReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;
    // The navigation stage (@todo: integrate break, target)
    Stage navigationStage = Stage::undefined;
    /// Force intersection with boundaries
    bool forceIntersectBoundaries = false;
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

  /// Reset state
  ///
  /// @param state is the state
  /// @param geoContext is the geometry context
  /// @param pos is the global position
  /// @param dir is the direction of navigation
  /// @param ssurface is the new starting surface
  /// @param tsurface is the target surface
  void resetState(State& state, const GeometryContext& geoContext,
                  const Vector3& pos, const Vector3& dir,
                  const Surface* ssurface, const Surface* tsurface) const {
    // Reset everything first
    state = State();

    // Set the start, current and target objects
    state.startSurface = ssurface;
    if (ssurface->associatedLayer() != nullptr) {
      state.startLayer = ssurface->associatedLayer();
    }
    if (state.startLayer->trackingVolume() != nullptr) {
      state.startVolume = state.startLayer->trackingVolume();
    }
    state.currentSurface = state.startSurface;
    state.currentVolume = state.startVolume;
    state.targetSurface = tsurface;

    // Get the compatible layers (including the current layer)
    NavigationOptions<Layer> navOpts;
    navOpts.resolveSensitive = true;
    navOpts.resolveMaterial = true;
    navOpts.resolvePassive = true;
    state.navLayers =
        state.currentVolume->compatibleLayers(geoContext, pos, dir, navOpts);

    // Set the index to the first
    state.navLayerIndex = 0;
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

  bool endOfWorldReached(State& state) const {
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

  /// @brief Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    // Set the world volume if it is not set
    if (not state.navigation.worldVolume) {
      state.navigation.worldVolume =
          m_cfg.trackingGeometry->highestTrackingVolume();
    }

    // We set the current surface to the start surface
    // for eventual post-update action, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    state.navigation.currentSurface = state.navigation.startSurface;
    if (state.navigation.currentSurface) {
      ACTS_VERBOSE(volInfo(state)
                   << "Current surface set to start surface "
                   << state.navigation.currentSurface->geometryId());
    }

    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.navigation.startSurface &&
        state.navigation.startSurface->associatedLayer()) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Surface.");
      // assign the current layer and volume by association
      state.navigation.startLayer =
          state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume =
          state.navigation.startLayer->trackingVolume();
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
    } else if (state.navigation.startVolume) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Volume.");
      state.navigation.startLayer =
          state.navigation.startVolume->associatedLayer(
              state.geoContext, stepper.position(state.stepping));
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Slow start initialization through search.");
      // current volume and layer search through global search
      ACTS_VERBOSE(volInfo(state)
                   << "Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
      state.navigation.startVolume =
          m_cfg.trackingGeometry->lowestTrackingVolume(
              state.geoContext, stepper.position(state.stepping));
      state.navigation.startLayer =
          state.navigation.startVolume
              ? state.navigation.startVolume->associatedLayer(
                    state.geoContext, stepper.position(state.stepping))
              : nullptr;
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.startVolume) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      }
    }
  }

  /// @brief Navigator pre step call
  ///
  /// Call options
  /// (a) there are still surfaces to be resolved: handle those
  /// (b) there no surfaces but still layers to be resolved, handle those
  /// (c) there are no surfaces nor layers to be resolved, handle boundary
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Check if the navigator is inactive
    if (inactive(state, stepper)) {
      return;
    }

    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::preStep.");

    // Initialize the target and target volume
    if (state.navigation.targetSurface and not state.navigation.targetVolume) {
      // Find out about the target as much as you can
      initializeTarget(state, stepper);
    }
    // Try targeting the surfaces - then layers - then boundaries
    if (state.navigation.navigationStage <= Stage::surfaceTarget and
        targetSurfaces(state, stepper)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
    } else if (state.navigation.navigationStage <= Stage::layerTarget and
               targetLayers(state, stepper)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
    } else if (targetBoundaries(state, stepper)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "No further navigation action, proceed to target.");
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping);
    }

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
  }

  /// @brief Navigator post step call, will be called in two modes
  ///
  /// (a) It initializes the Navigation stream if start volume is
  ///     not yet defined:
  ///  - initialize the volume
  ///  - establish the start layer and start volume
  ///  - set the current surface to the start surface
  ///
  /// (b) It establishes the currentSurface status during
  ///     the propagation flow, currentSurface can be
  ///  - surfaces still to be handled within a layer
  ///  - layers still to be handled within a volume
  ///  - boundaries still to be handled to exit a volume
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Check if the navigator is inactive
    if (inactive(state, stepper)) {
      return;
    }

    // Set the navigation stage
    state.navigation.navigationStage = Stage::undefined;

    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::postStep.");

    // Navigator post step always starts without current surface
    state.navigation.currentSurface = nullptr;

    // (b) Status call within propagation loop
    // Try finding status of surfaces
    if (surfaceStatus(state, stepper, state.navigation.navSurfaces,
                      state.navigation.navSurfaceIndex)) {
      ACTS_VERBOSE(volInfo(state) << "Post step: in surface handling.");
      if (state.navigation.currentSurface) {
        ACTS_VERBOSE(volInfo(state)
                     << "On surface: switch forward or release.");
        if (++state.navigation.navSurfaceIndex ==
            state.navigation.navSurfaces.size()) {
          // this was the last surface, check if we have layers
          if (!state.navigation.navLayers.empty()) {
            ++state.navigation.navLayerIndex;
          } else if (state.navigation.startLayer != nullptr and
                     state.navigation.currentSurface->associatedLayer() ==
                         state.navigation.startLayer) {
            // this was the start layer, switch to layer target next
            state.navigation.navigationStage = Stage::layerTarget;
            return;
          } else {
            // no layers, go to boundary
            state.navigation.navigationStage = Stage::boundaryTarget;
            return;
          }
        }
      }
      // Set the navigation stage to surface target
      state.navigation.navigationStage = Stage::surfaceTarget;
      ACTS_VERBOSE(volInfo(state) << "Staying focussed on surface.");
      // Try finding status of layer
    } else if (surfaceStatus(state, stepper, state.navigation.navLayers,
                             state.navigation.navLayerIndex)) {
      ACTS_VERBOSE(volInfo(state) << "Post step: in layer handling.");
      if (state.navigation.currentSurface != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "On layer: update layer information.");
        if (resolveSurfaces(state, stepper)) {
          // Set the navigation stage back to surface handling
          state.navigation.navigationStage = Stage::surfaceTarget;
          return;
        }
      } else {
        // Set the navigation stage to layer target
        state.navigation.navigationStage = Stage::layerTarget;
        ACTS_VERBOSE(volInfo(state) << "Staying focussed on layer.");
      }
      // Try finding status of boundaries
    } else if (surfaceStatus(state, stepper, state.navigation.navBoundaries,
                             state.navigation.navBoundaryIndex)) {
      ACTS_VERBOSE(volInfo(state) << "Post step: in boundary handling.");

      // Are we on the boundary - then overwrite the stage
      if (state.navigation.currentSurface != nullptr) {
        // Set the navigation stage back to surface handling
        ACTS_VERBOSE(volInfo(state)
                     << "On boundary: update volume information.");
        // We are on a boundary, reset all information
        state.navigation.navSurfaces.clear();
        state.navigation.navSurfaceIndex = state.navigation.navSurfaces.size();
        state.navigation.navLayers.clear();
        state.navigation.navLayerIndex = state.navigation.navLayers.size();
        state.navigation.lastHierarchySurfaceReached = false;
        // Update volume information
        // get the attached volume information
        auto boundary = state.navigation.navBoundary().object();
        state.navigation.currentVolume = boundary->attachedVolume(
            state.geoContext, stepper.position(state.stepping),
            stepper.direction(state.stepping), state.options.direction);
        // No volume anymore : end of known world
        if (!state.navigation.currentVolume) {
          ACTS_VERBOSE(
              volInfo(state)
              << "No more volume to progress to, stopping navigation.");
          // Navigation break & release navigation stepping
          state.navigation.navigationBreak = true;
          stepper.releaseStepSize(state.stepping);
          return;
        } else {
          ACTS_VERBOSE(volInfo(state) << "Volume updated.");
          // Forget the boundary information
          state.navigation.navBoundaries.clear();
          state.navigation.navBoundaryIndex =
              state.navigation.navBoundaries.size();
        }
      } else {
        // Set the navigation stage back to boundary target
        state.navigation.navigationStage = Stage::boundaryTarget;
        ACTS_VERBOSE(volInfo(state) << "Staying focussed on boundary.");
      }
    } else if (state.navigation.currentVolume ==
               state.navigation.targetVolume) {
      if (state.navigation.targetSurface == nullptr) {
        ACTS_WARNING(volInfo(state)
                     << "No further navigation action, proceed to "
                        "target. This is very likely an error");
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "No further navigation action, proceed to target.");
      }
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping);
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Status could not be determined - good luck.");
    }
  }

 private:
  /// @brief Status call for test surfaces (surfaces, layers, boundaries)
  ///
  /// If there are surfaces to be handled, check if the current
  /// state is on the surface
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  /// @tparam navigation_surfaces_t Type of the propagagor
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  /// @param [in] navSurfaces is the navigation status objects
  /// @param [in] navIndex test surface fore the status test
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t,
            typename navigation_surfaces_t>
  bool surfaceStatus(propagator_state_t& state, const stepper_t& stepper,
                     const navigation_surfaces_t& navSurfaces,
                     std::size_t navIndex) const {
    // No surfaces, status check will be done on layer
    if (navSurfaces.empty() or navIndex == navSurfaces.size()) {
      return false;
    }
    // Take the current surface
    auto surface = navSurfaces.at(navIndex).representation();
    // Check if we are at a surface
    // If we are on the surface pointed at by the index, we can make
    // it the current one to pass it to the other actors
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, *surface, state.options.direction, true,
        state.options.targetTolerance, logger());
    if (surfaceStatus == Intersection3D::Status::onSurface) {
      ACTS_VERBOSE(volInfo(state)
                   << "Status Surface successfully hit, storing it.");
      // Set in navigation state, so actors and aborters can access it
      state.navigation.currentSurface = surface;
      if (state.navigation.currentSurface) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current surface set to surface "
                     << state.navigation.currentSurface->geometryId());
      }
    }
    // Return a positive status: either on it, or on the way
    return true;
  }

  /// Loop over surface candidates here:
  ///  - if an intersect is  valid but not yet reached
  ///    then return with updated step size
  ///  - if an intersect is not valid, switch to next
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool targetSurfaces(propagator_state_t& state,
                      const stepper_t& stepper) const {
    if (state.navigation.navigationBreak) {
      return false;
    }
    // Make sure resolve Surfaces is called on the start layer
    if (state.navigation.startLayer and
        not state.navigation.startLayerResolved) {
      ACTS_VERBOSE(volInfo(state) << "Start layer to be resolved.");
      // We provide the layer to the resolve surface method in this case
      state.navigation.startLayerResolved = true;
      bool startResolved =
          resolveSurfaces(state, stepper, state.navigation.startLayer);
      if (not startResolved and
          state.navigation.startLayer == state.navigation.targetLayer) {
        ACTS_VERBOSE(volInfo(state)
                     << "Start is target layer, nothing left to do.");
        // set the navigation break
        state.navigation.navigationBreak = true;
        stepper.releaseStepSize(state.stepping);
      }
      return startResolved;
    }

    // The call that we are on a layer and have not yet resolved the surfaces
    // No surfaces, do not return to stepper
    if (state.navigation.navSurfaces.empty() or
        state.navigation.navSurfaceIndex ==
            state.navigation.navSurfaces.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No surfaces present, target at layer first.");
      return false;
    }
    auto layerID = state.navigation.navSurface().object()->geometryId().layer();
    std::pair<ExternalSurfaces::iterator, ExternalSurfaces::iterator>
        externalSurfaceRange =
            state.navigation.externalSurfaces.equal_range(layerID);
    // Loop over the remaining navigation surfaces
    while (state.navigation.navSurfaceIndex !=
           state.navigation.navSurfaces.size()) {
      // Screen output how much is left to try
      ACTS_VERBOSE(volInfo(state)
                   << (state.navigation.navSurfaces.size() -
                       state.navigation.navSurfaceIndex)
                   << " out of " << state.navigation.navSurfaces.size()
                   << " surfaces remain to try.");
      // Take the surface
      auto surface = state.navigation.navSurface().object();
      // Screen output which surface you are on
      ACTS_VERBOSE(volInfo(state) << "Next surface candidate will be "
                                  << surface->geometryId());
      // Estimate the surface status
      bool boundaryCheck = true;
      for (auto it = externalSurfaceRange.first;
           it != externalSurfaceRange.second; it++) {
        if (surface->geometryId() == it->second) {
          boundaryCheck = false;
          break;
        }
      }
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, *surface, state.options.direction, boundaryCheck,
          state.options.targetTolerance, logger());
      if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE(volInfo(state)
                     << "Surface reachable, step size updated to "
                     << stepper.outputStepSize(state.stepping));
        return true;
      }
      ++state.navigation.navSurfaceIndex;
      continue;
    }

    // Reached the end of the surface iteration
    if (state.navigation.navSurfaceIndex ==
        state.navigation.navSurfaces.size()) {
      // first clear the surface cache
      state.navigation.navSurfaces.clear();
      state.navigation.navSurfaceIndex = state.navigation.navSurfaces.size();

      if (state.navigation.navLayerIndex != state.navigation.navLayers.size()) {
        ACTS_VERBOSE(volInfo(state)
                     << "Last surface on layer reached, switching layer.");
        // now switch to the next layer
        ++state.navigation.navLayerIndex;
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "Last surface on layer reached, and no layer.");
        // first clear the surface cache
        state.navigation.lastHierarchySurfaceReached = true;
        state.navigation.navigationBreak =
            (state.navigation.currentVolume == state.navigation.targetVolume);
      }
    }
    // Do not return to the propagator
    return false;
  }

  /// @brief Target layer candidates.
  ///
  /// We are now trying to advance to the next layer (with surfaces)
  /// Check if we are on the representing surface of the layer pointed
  /// at by navLayerIndex. If so, we unpack the compatible surfaces
  /// (determined by straight line intersect), and set up the index
  /// so that the next postStep() call will enter the surface
  /// check mode above. If no surfaces are found, we skip the layer.
  /// If we unpack a surface, the step size is set to the path length
  /// to the first surface, as determined by straight line intersect.
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool targetLayers(propagator_state_t& state, const stepper_t& stepper) const {
    using namespace UnitLiterals;

    if (state.navigation.navigationBreak ||
        state.navigation.lastHierarchySurfaceReached) {
      return false;
    }

    // if there are no layers, go back to the navigator (not stepper yet)
    if (state.navigation.navLayers.empty()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No layers present, resolve volume first.");

      // check if current volume has BVH, or layers
      if (state.navigation.currentVolume->hasBoundingVolumeHierarchy()) {
        // has hierarchy, use that, skip layer resolution
        NavigationOptions<Surface> navOpts;
        navOpts.resolveSensitive = m_cfg.resolveSensitive;
        navOpts.resolveMaterial = m_cfg.resolveMaterial;
        navOpts.resolvePassive = m_cfg.resolvePassive;
        navOpts.endObject = state.navigation.targetSurface;
        navOpts.overstepLimit = stepper.overstepLimit(state.stepping);
        double opening_angle = 0;

        // Preliminary version of the frustum opening angle estimation.
        // Currently not used (only rays), but will be.

        /*
        Vector3 pos = stepper.position(state.stepping);
        double mom = stepper.momentum(state.stepping) / UnitConstants::GeV;
        double q = stepper.charge(state.stepping);
        Vector3 dir = stepper.direction(state.stepping);
        Vector3 B = stepper.getField(state.stepping, pos);
        if (B.squaredNorm() > 1e-9) {
          // ~ non-zero field
          double ir = (dir.cross(B).norm()) * q / mom;
          double s;
          if (state.options.direction == Direction::Forward) {
            s = state.stepping.stepSize.max();
          } else {
            s = state.stepping.stepSize.min();
          }
          opening_angle = std::atan((1 - std::cos(s * ir)) / std::sin(s * ir));
        }

        ACTS_VERBOSE(volInfo(state) << "Estimating opening angle for frustum
        nav:"); ACTS_VERBOSE(volInfo(state) << "pos: " << pos.transpose());
        ACTS_VERBOSE(volInfo(state) << "dir: " << dir.transpose());
        ACTS_VERBOSE(volInfo(state) << "B: " << B.transpose() << " |B|: " <<
        B.norm()); ACTS_VERBOSE(volInfo(state) << "step mom: " <<
        stepper.momentum(state.stepping)); ACTS_VERBOSE(volInfo(state) << "=>
        opening angle: " << opening_angle);
        */

        auto protoNavSurfaces =
            state.navigation.currentVolume->compatibleSurfacesFromHierarchy(
                state.geoContext, stepper.position(state.stepping),
                state.options.direction * stepper.direction(state.stepping),
                opening_angle, navOpts);
        if (!protoNavSurfaces.empty()) {
          // did we find any surfaces?

          // Check: are we on the first surface?
          // TODO magic number `1_um`
          if ((state.navigation.currentSurface == nullptr &&
               state.navigation.navSurfaces.empty()) ||
              protoNavSurfaces.front().pathLength() > 1_um) {
            // we are not, go on
            // state.navigation.navSurfaces = std::move(protoNavSurfaces);
            state.navigation.navSurfaces.clear();
            state.navigation.navSurfaces.insert(
                state.navigation.navSurfaces.begin(), protoNavSurfaces.begin(),
                protoNavSurfaces.end());

            state.navigation.navSurfaceIndex = 0;
            state.navigation.navLayers = {};
            state.navigation.navLayerIndex = state.navigation.navLayers.size();
            // The stepper updates the step size ( single / multi component)
            stepper.updateStepSize(state.stepping,
                                   state.navigation.navSurface(),
                                   state.options.direction, true);
            ACTS_VERBOSE(volInfo(state)
                         << "Navigation stepSize updated to "
                         << stepper.outputStepSize(state.stepping));
            return true;
          }
        }
      }

      if (resolveLayers(state, stepper)) {
        // The layer resolving worked
        return true;
      }
    }
    // loop over the available navigation layer candidates
    while (state.navigation.navLayerIndex !=
           state.navigation.navLayers.size()) {
      // The layer surface
      auto layerSurface = state.navigation.navLayer().representation();
      // We are on the layer
      if (state.navigation.currentSurface == layerSurface) {
        ACTS_VERBOSE(volInfo(state) << "We are on a layer, resolve Surfaces.");
        // If you found surfaces return to the propagator
        if (resolveSurfaces(state, stepper)) {
          return true;
        } else {
          // Try the next one
          ++state.navigation.navLayerIndex;
          continue;
        }
      }
      // Try to step towards it
      auto layerStatus = stepper.updateSurfaceStatus(
          state.stepping, *layerSurface, state.options.direction, true,
          state.options.targetTolerance, logger());
      if (layerStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE(volInfo(state) << "Layer reachable, step size updated to "
                                    << stepper.outputStepSize(state.stepping));
        return true;
      }
      ACTS_VERBOSE(volInfo(state)
                   << "Layer intersection not valid, skipping it.");
      ++state.navigation.navLayerIndex;
    }

    // Re-initialize target at last layer, only in case it is the target volume
    // This avoids a wrong target volume estimation
    if (state.navigation.currentVolume == state.navigation.targetVolume) {
      initializeTarget(state, stepper);
    }
    // Screen output
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << "Last layer";
      if (state.navigation.currentVolume == state.navigation.targetVolume) {
        os << " (final volume) done, proceed to target.";
      } else {
        os << " done, target volume boundary.";
      }
      logger().log(Logging::VERBOSE, os.str());
    }
    // Set the navigation break if necessary
    state.navigation.navigationBreak =
        (state.navigation.currentVolume == state.navigation.targetVolume);
    return false;
  }

  /// @brief Navigation through volumes
  ///
  /// This is the boundary check routine. If the code above set up the
  /// boundary surface index, we advance through them here. If we are on
  /// the boundary surface, we set the current surface to the boundary
  /// surface, and get the volume pointed at by the boundary surface.  Next
  /// we unpack the layers from that volume. If the volume contains layers
  /// we set the step size to the straight line path length to the first
  /// layer.  If we don't find a next volume, the navigationBreak
  /// indicator is set.  This ends the navigation. Finally, the boundary
  /// index is cleared, so that the subsequent call goes back to
  /// the layer iteration logic.
  ///
  /// If we are not on the current boundary surface, we try the next one.
  /// The index is advanced and the step size is set. If no straight
  /// line intersect is found, the boundary surface is skipped.
  /// If we are out of boundary surfaces, the navigation is terminated.
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool targetBoundaries(propagator_state_t& state,
                        const stepper_t& stepper) const {
    if (state.navigation.navigationBreak) {
      return false;
    }

    if (!state.navigation.currentVolume) {
      ACTS_VERBOSE(volInfo(state)
                   << "No sufficient information to resolve boundary, "
                      "stopping navigation.");
      stepper.releaseStepSize(state.stepping);
      return false;
    } else if (state.navigation.currentVolume ==
               state.navigation.targetVolume) {
      ACTS_VERBOSE(volInfo(state)
                   << "In target volume: no need to resolve boundary, "
                      "stopping navigation.");
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping);
      return true;
    }

    // Re-initialize target at volume boundary
    initializeTarget(state, stepper);

    // Helper function to find boundaries
    auto findBoundaries = [&]() -> bool {
      // The navigation options
      NavigationOptions<Surface> navOpts;
      // Exclude the current surface in case it's a boundary
      navOpts.startObject = state.navigation.currentSurface;
      navOpts.pathLimit =
          stepper.getStepSize(state.stepping, ConstrainedStep::aborter);
      navOpts.overstepLimit = stepper.overstepLimit(state.stepping);
      navOpts.forceIntersectBoundaries =
          state.navigation.forceIntersectBoundaries;

      ACTS_VERBOSE(volInfo(state)
                   << "Try to find boundaries, we are at: "
                   << stepper.position(state.stepping).transpose() << ", dir: "
                   << stepper.direction(state.stepping).transpose());

      // Evaluate the boundary surfaces
      state.navigation.navBoundaries =
          state.navigation.currentVolume->compatibleBoundaries(
              state.geoContext, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping),
              navOpts, logger());
      // The number of boundary candidates
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << state.navigation.navBoundaries.size();
        os << " boundary candidates found at path(s): ";
        for (auto& bc : state.navigation.navBoundaries) {
          os << bc.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }
      // Set the begin index
      state.navigation.navBoundaryIndex = 0;
      if (not state.navigation.navBoundaries.empty()) {
        // Set to the first and return to the stepper
        stepper.updateStepSize(state.stepping, state.navigation.navBoundary(),
                               state.options.direction, true);
        ACTS_VERBOSE(volInfo(state) << "Navigation stepSize updated to "
                                    << stepper.outputStepSize(state.stepping));
        return true;
      }
      return false;
    };

    // No boundaries are assigned yet, find them
    if (state.navigation.navBoundaries.empty() and findBoundaries()) {
      return true;
    }

    // Loop over the boundary surface
    while (state.navigation.navBoundaryIndex !=
           state.navigation.navBoundaries.size()) {
      // That is the current boundary surface
      auto boundarySurface = state.navigation.navBoundary().representation();
      // Step towards the boundary surfrace
      auto boundaryStatus = stepper.updateSurfaceStatus(
          state.stepping, *boundarySurface, state.options.direction, true,
          state.options.targetTolerance, logger());
      if (boundaryStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE(volInfo(state)
                     << "Boundary reachable, step size updated to "
                     << stepper.outputStepSize(state.stepping));
        return true;
      } else {
        ACTS_VERBOSE("Boundary "
                     << (state.navigation.navBoundaries.size() -
                         state.navigation.navBoundaryIndex)
                     << " out of " << state.navigation.navBoundaries.size()
                     << " not reachable anymore, switching to next.");
        ACTS_VERBOSE("Targeted boundary surface was: \n"
                     << std::tie(*boundarySurface, state.geoContext));
      }
      // Increase the index to the next one
      ++state.navigation.navBoundaryIndex;
    }
    // We have to leave the volume somehow, so try again
    state.navigation.navBoundaries.clear();
    ACTS_VERBOSE(volInfo(state) << "Boundary navigation lost, re-targetting.");
    state.navigation.forceIntersectBoundaries = true;
    if (findBoundaries()) {
      // Resetting intersection check for boundary surfaces
      state.navigation.forceIntersectBoundaries = false;
      return true;
    }

    // Tried our best, but couldn't do anything
    return false;
  }

  /// @brief Navigation (re-)initialisation for the target
  ///
  /// @note This is only called a few times every propagation/extrapolation
  ///
  /// As a straight line estimate can lead you to the wrong destination
  /// Volume, this will be called at:
  /// - initialization
  /// - attempted volume switch
  /// Target finding by association will not be done again
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  void initializeTarget(propagator_state_t& state,
                        const stepper_t& stepper) const {
    if (state.navigation.targetVolume and
        state.stepping.pathAccumulated == 0.) {
      ACTS_VERBOSE(volInfo(state)
                   << "Re-initialzing cancelled as it is the first step.");
      return;
    }
    // Fast Navigation initialization for target:
    if (state.navigation.targetSurface &&
        state.navigation.targetSurface->associatedLayer() &&
        !state.navigation.targetVolume) {
      ACTS_VERBOSE(volInfo(state)
                   << "Fast target initialization through association.");
      ACTS_VERBOSE(volInfo(state)
                   << "Target surface set to "
                   << state.navigation.targetSurface->geometryId());
      state.navigation.targetLayer =
          state.navigation.targetSurface->associatedLayer();
      state.navigation.targetVolume =
          state.navigation.targetLayer->trackingVolume();
    } else if (state.navigation.targetSurface) {
      // screen output that you do a re-initialization
      if (state.navigation.targetVolume) {
        ACTS_VERBOSE(volInfo(state)
                     << "Re-initialization of target volume triggered.");
      }
      // Slow navigation initialization for target:
      // target volume and layer search through global search
      auto targetIntersection =
          state.navigation.targetSurface
              ->intersect(
                  state.geoContext, stepper.position(state.stepping),
                  state.options.direction * stepper.direction(state.stepping),
                  false, state.options.targetTolerance)
              .closest();
      if (targetIntersection) {
        ACTS_VERBOSE(volInfo(state)
                     << "Target estimate position ("
                     << targetIntersection.position().x() << ", "
                     << targetIntersection.position().y() << ", "
                     << targetIntersection.position().z() << ")");
        /// get the target volume from the intersection
        auto tPosition = targetIntersection.position();
        state.navigation.targetVolume =
            m_cfg.trackingGeometry->lowestTrackingVolume(state.geoContext,
                                                         tPosition);
        state.navigation.targetLayer =
            state.navigation.targetVolume
                ? state.navigation.targetVolume->associatedLayer(
                      state.geoContext, tPosition)
                : nullptr;
        if (state.navigation.targetVolume) {
          ACTS_VERBOSE(volInfo(state)
                       << "Target volume estimated : "
                       << state.navigation.targetVolume->volumeName());
        }
      }
    }
  }

  /// @brief Resolve the surfaces of this layer
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  /// @param [in] cLayer is the already assigned current layer, e.g. start layer
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool resolveSurfaces(propagator_state_t& state, const stepper_t& stepper,
                       const Layer* cLayer = nullptr) const {
    // get the layer and layer surface
    auto layerSurface = cLayer ? state.navigation.startSurface
                               : state.navigation.navLayer().representation();
    auto navLayer = cLayer ? cLayer : state.navigation.navLayer().object();
    // are we on the start layer
    bool onStart = (navLayer == state.navigation.startLayer);
    auto startSurface = onStart ? state.navigation.startSurface : layerSurface;
    // Use navigation parameters and NavigationOptions
    NavigationOptions<Surface> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = startSurface;
    navOpts.endObject = state.navigation.targetSurface;

    std::vector<GeometryIdentifier> externalSurfaces;
    if (!state.navigation.externalSurfaces.empty()) {
      auto layerID = layerSurface->geometryId().layer();
      auto externalSurfaceRange =
          state.navigation.externalSurfaces.equal_range(layerID);
      navOpts.externalSurfaces.reserve(
          state.navigation.externalSurfaces.count(layerID));
      for (auto itSurface = externalSurfaceRange.first;
           itSurface != externalSurfaceRange.second; itSurface++) {
        navOpts.externalSurfaces.push_back(itSurface->second);
      }
    }
    // Check the limit
    navOpts.pathLimit =
        stepper.getStepSize(state.stepping, ConstrainedStep::aborter);
    // No overstepping on start layer, otherwise ask the stepper
    navOpts.overstepLimit = (cLayer != nullptr)
                                ? state.options.targetTolerance
                                : stepper.overstepLimit(state.stepping);

    // get the surfaces
    state.navigation.navSurfaces = navLayer->compatibleSurfaces(
        state.geoContext, stepper.position(state.stepping),
        state.options.direction * stepper.direction(state.stepping), navOpts);
    // the number of layer candidates
    if (!state.navigation.navSurfaces.empty()) {
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << state.navigation.navSurfaces.size();
        os << " surface candidates found at path(s): ";
        for (auto& sfc : state.navigation.navSurfaces) {
          os << sfc.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }

      // set the index
      state.navigation.navSurfaceIndex = 0;
      // The stepper updates the step size ( single / multi component)
      stepper.updateStepSize(state.stepping, state.navigation.navSurface(),
                             state.options.direction, true);
      ACTS_VERBOSE(volInfo(state) << "Navigation stepSize updated to "
                                  << stepper.outputStepSize(state.stepping));
      return true;
    }
    state.navigation.navSurfaceIndex = state.navigation.navSurfaces.size();
    ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
    return false;
  }

  /// @brief Navigation through layers
  ///
  /// Resolve layers.
  ///
  /// This initializes the layer candidates when starting
  /// or when entering a new volume
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool resolveLayers(propagator_state_t& state,
                     const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

    // Check if we are in the start volume
    auto startLayer =
        (state.navigation.currentVolume == state.navigation.startVolume)
            ? state.navigation.startLayer
            : nullptr;
    // Create the navigation options
    // - and get the compatible layers, start layer will be excluded
    NavigationOptions<Layer> navOpts;
    navOpts.boundaryCheck = m_cfg.boundaryCheckLayerResolving;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = startLayer;
    // Set also the target surface
    navOpts.targetSurface = state.navigation.targetSurface;
    navOpts.pathLimit =
        stepper.getStepSize(state.stepping, ConstrainedStep::aborter);
    navOpts.overstepLimit = stepper.overstepLimit(state.stepping);
    // Request the compatible layers
    state.navigation.navLayers =
        state.navigation.currentVolume->compatibleLayers(
            state.geoContext, stepper.position(state.stepping),
            state.options.direction * stepper.direction(state.stepping),
            navOpts);

    // Layer candidates have been found
    if (!state.navigation.navLayers.empty()) {
      // Screen output where they are
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << state.navigation.navLayers.size();
        os << " layer candidates found at path(s): ";
        for (auto& lc : state.navigation.navLayers) {
          os << lc.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }
      // Set the index to the first
      state.navigation.navLayerIndex = 0;
      // Setting the step size towards first
      if (state.navigation.startLayer &&
          state.navigation.navLayer().object() != state.navigation.startLayer) {
        ACTS_VERBOSE(volInfo(state) << "Target at layer.");
        // The stepper updates the step size ( single / multi component)
        stepper.updateStepSize(state.stepping, state.navigation.navLayer(),
                               state.options.direction, true);
        ACTS_VERBOSE(volInfo(state) << "Navigation stepSize updated to "
                                    << stepper.outputStepSize(state.stepping));
        return true;
      }
    }

    // Set the index to the end of the list
    state.navigation.navLayerIndex = state.navigation.navLayers.size();

    // Screen output - no layer candidates found
    ACTS_VERBOSE(volInfo(state) << "No compatible layer candidates found.");
    // Release the step size
    stepper.releaseStepSize(state.stepping);
    return false;
  }

  /// Inactive
  ///
  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool inactive(propagator_state_t& state, const stepper_t& stepper) const {
    // Void behavior in case no tracking geometry is present
    if (!m_cfg.trackingGeometry) {
      return true;
    }
    // turn the navigator into void when you are intructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    // Navigation break handling
    // This checks if a navigation break had been triggered:
    // - If so & the target exists or was hit - it simply returns
    // - If a target exists and was not yet hit, it checks for it
    // -> return is always to the stepper
    if (state.navigation.navigationBreak) {
      // target exists and reached, or no target exists
      if (state.navigation.targetReached || !state.navigation.targetSurface) {
        return true;
      }
      auto targetStatus = stepper.updateSurfaceStatus(
          state.stepping, *state.navigation.targetSurface,
          state.options.direction, true, state.options.targetTolerance,
          logger());
      // the only advance could have been to the target
      if (targetStatus == Intersection3D::Status::onSurface) {
        // set the target surface
        state.navigation.currentSurface = state.navigation.targetSurface;
        ACTS_VERBOSE(volInfo(state)
                     << volInfo(state)
                     << "Current surface set to target surface "
                     << state.navigation.currentSurface->geometryId());
        return true;
      }
    }
    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume
                ? state.navigation.currentVolume->volumeName()
                : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
