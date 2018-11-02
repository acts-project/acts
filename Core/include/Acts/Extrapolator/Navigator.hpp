
// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iterator>
#include <sstream>
#include <string>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"
#include <boost/algorithm/string.hpp>

namespace Acts {

using Cstep = detail::ConstrainedStep;
/// @breif struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions
{

  /// The navigation direction
  NavigationDirection navDir = forward;

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

  double pathLimit = std::numeric_limits<double>::max();

  /// Constructor
  ///
  /// @param nDir Navigation direction prescription
  /// @param bcheck Boundary check for the navigation action
  /// @param sobject Start object to check against
  /// @param eobject End object to check against
  /// @param maxStepLength Maximal step length to check against
  NavigationOptions(NavigationDirection ndir,
                    BoundaryCheck       bcheck,
                    bool                resolves = true,
                    bool                resolvem = true,
                    bool                resolvep = false,
                    const object_t*     sobject  = nullptr,
                    const object_t*     eobject  = nullptr)
    : navDir(ndir)
    , boundaryCheck(std::move(bcheck))
    , resolveSensitive(resolves)
    , resolveMaterial(resolvem)
    , resolvePassive(resolvep)
    , startObject(sobject)
    , endObject(eobject)
    , pathLimit(ndir * std::numeric_limits<double>::max())
  {
  }
};

/// Navigator class
///
/// This is an Actor to be added to the ActorList in order to navigate
/// through the static tracking geometry setup.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation  flow is hit, it is
/// set to the propagation satate, such that other actors can deal wit it.
/// This navigation actor thus always needs to run first!
/// It does two things: it figures out the order of volumes, layers and
/// surfaces. For each propagation step, the operator() runs, which checks if
/// the current surface (or layer/volume boundary) is reached via isOnSurface.
/// The current target surface is the surface pointed to by of the iterators
/// for the surfaces, layers or volume boundaries.
/// If a surface is found, the state.navigation.currentSurface
/// pointer is set. This  enables subsequent actors to react. Secondly, this
/// actor uses the ordered  iterators  to figure out which surface, layer or
/// volume boundary is _supposed_ to be hit next. It then sets the maximum
/// step size to the path length found out by straight line intersection. If
/// the isOnSurface call fails, it also  re-computes the step size, to make
/// sure we end up at the desired surface.
///
class Navigator
{

public:
  using Surfaces    = std::vector<const Surface*>;
  using SurfaceIter = std::vector<const Surface*>::iterator;

  using NavigationSurfaces    = std::vector<SurfaceIntersection>;
  using NavigationSurfaceIter = NavigationSurfaces::iterator;

  using NavigationLayers    = std::vector<LayerIntersection>;
  using NavigationLayerIter = NavigationLayers::iterator;

  using NavigationBoundaries   = std::vector<BoundaryIntersection>;
  using NavigationBoundaryIter = NavigationBoundaries::iterator;

  using ExternalSurfaces = std::multimap<const Layer*, const Surface*>;

  /// The navigation stage
  enum struct Stage : int {
    undefined      = 0,
    surfaceTarget  = 1,
    layerTarget    = 2,
    boundaryTarget = 3
  };

  /// Constructor with shared tracking geometry
  ///
  /// @param tGeometry The tracking geometry for the navigator
  Navigator(std::shared_ptr<const TrackingGeometry> tGeometry = nullptr)
    : trackingGeometry(std::move(tGeometry))
  {
  }

  /// Tracking Geometry for this Navigator
  std::shared_ptr<const TrackingGeometry> trackingGeometry;

  /// The tolerance used to defined "reached"
  double tolerance = s_onSurfaceTolerance;

  /// Configuration for this Navigator
  /// stop at every sensitive surface (whether it has material or not)
  bool resolveSensitive = true;
  /// stop at every material surface (whether it is passive or not)
  bool resolveMaterial = true;
  /// stop at every surface regardless what it is
  bool resolvePassive = false;

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State
  {
    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface iterator of the navigation state
    NavigationSurfaceIter navSurfaceIter = navSurfaces.end();

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer iterator of the navigation state
    NavigationLayerIter navLayerIter = navLayers.end();

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary iterator of the navigation state
    NavigationBoundaryIter navBoundaryIter = navBoundaries.end();

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    /// Navigation sate: the world volume
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
    /// Navigation state : a break has been detected
    bool navigationBreak = false;
    // The navigation stage (@todo: integrate break, target)
    Stage navigationStage = Stage::undefined;
  };

  /// Unique typedef to publish to the Propagator
  using state_type = State;

  /// @brief Navigator status call, will be called in two modes
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
  ///
  /// @param[in,out] state is the mutable propagator state object
  template <typename propagator_state_t>
  void
  status(propagator_state_t& state) const
  {
    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    // Set the navigation stage
    state.navigation.navigationStage = Stage::undefined;

    // Call the navigation helper prior to actual navigation
    debugLog(state, [&] { return std::string("Entering navigator::status."); });

    // (a) Pre-stepping call from propgator
    if (not state.navigation.startVolume) {
      // Initialize and return
      initialize(state);
      return;
    }

    // Navigator status always starts without current surface
    state.navigation.currentSurface = nullptr;

    // (b) Status call within propagation loop
    // Try finding status of surfaces
    if (status(state,
               state.navigation.navSurfaces,
               state.navigation.navSurfaceIter)) {
      debugLog(state,
               [&] { return std::string("Status: in surface handling."); });
      if (state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("On surface: switch forward or release.");
        });
        if (++state.navigation.navSurfaceIter
            == state.navigation.navSurfaces.end()) {
          ++state.navigation.navLayerIter;
        }
      }
      // Set the navigation stage to surface target
      state.navigation.navigationStage = Stage::surfaceTarget;
      // Try finding status of layer
    } else if (status(state,
                      state.navigation.navLayers,
                      state.navigation.navLayerIter)) {
      debugLog(state,
               [&] { return std::string("Status: in layer handling."); });
      if (state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("On layer: update layer information.");
        });
        // Get a  navigation corrector associated to the stepper
        auto navCorr = state.stepping.corrector();
        if (resolveSurfaces(state, navCorr)) {
          // Set the navigation stage back to surface handling
          state.navigation.navigationStage = Stage::surfaceTarget;
          return;
        }
      } else {
        // Set the navigation stage to layer target
        state.navigation.navigationStage = Stage::layerTarget;
      }
      // Try finding status of boundaries
    } else if (status(state,
                      state.navigation.navBoundaries,
                      state.navigation.navBoundaryIter)) {
      debugLog(state,
               [&] { return std::string("Stauts: in boundary handling."); });

      // Are we on the boundary - then overwrite the stage
      if (state.navigation.currentSurface) {
        // Set the navigation stage back to surface handling
        debugLog(state, [&] {
          return std::string("On boundary: update volume information.");
        });
        // We are on a boundary, reset all information
        state.navigation.navSurfaces.clear();
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
        state.navigation.navLayers.clear();
        state.navigation.navLayerIter = state.navigation.navLayers.end();
        // Update volume information
        // get the attached volume information
        auto boundary = state.navigation.navBoundaryIter->object;
        state.navigation.currentVolume
            = boundary->attachedVolume(state.stepping.position(),
                                       state.stepping.direction(),
                                       state.stepping.navDir);
        // No volume anymore : end of known world
        if (!state.navigation.currentVolume) {
          debugLog(state, [&] {
            return std::string(
                "No more volume to progress to, stopping navigation.");
          });
          // Navigation break & release navigation stepping
          state.navigation.navigationBreak = true;
          state.stepping.stepSize.release(Cstep::actor);
          return;
        } else {
          debugLog(state, [&] { return std::string("Volume updated."); });
          // Forget the bounday information
          state.navigation.navBoundaries.clear();
          state.navigation.navBoundaryIter
              = state.navigation.navBoundaries.end();
        }
      } else {
        // Set the navigation stage back to boundary target
        state.navigation.navigationStage = Stage::boundaryTarget;
      }
    } else if (state.navigation.currentVolume
               == state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string("No further navigation action, proceed to target.");
      });
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      state.stepping.stepSize.release(Cstep::actor);
    } else {
      debugLog(state, [&] {
        return std::string("Status could not be determined - good luck.");
      });
    }
    return;
  }

  /// @brief Navigator target call
  ///
  /// Call options
  /// (a) there are still surfaces to be resolved: handle those
  /// (b) there no surfaces but still layers to be resolved, handle those
  /// (c) there are no surfaces nor layers to be resolved, handle boundary
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  ///
  /// @param[in,out] state is the mutable propagator state object
  template <typename propagator_state_t>
  void
  target(propagator_state_t& state) const
  {
    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    // Call the navigation helper prior to actual navigation
    debugLog(state, [&] { return std::string("Entering navigator::target."); });

    // Get a  navigation corrector associated to the stepper
    auto navCorr = state.stepping.corrector();

    // Initialize the target and target volume
    if (state.navigation.targetSurface and not state.navigation.targetVolume) {
      // Find out about the target as much as you can
      initializeTarget(state, navCorr);
    }

    // Try targeting the surfaces - then layers - then boundaries
    if (state.navigation.navigationStage <= Stage::surfaceTarget
        and targetSurfaces(state, navCorr)) {
      debugLog(state,
               [&] { return std::string("Target set to next surface."); });
    } else if (state.navigation.navigationStage <= Stage::layerTarget
               and targetLayers(state, navCorr)) {
      debugLog(state, [&] { return std::string("Target set to next layer."); });
    } else if (targetBoundaries(state, navCorr)) {
      debugLog(state,
               [&] { return std::string("Target set to next boundary."); });
    } else {
      debugLog(state, [&] {
        return std::string("No furter navigation action, proceed to target.");
      });
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      state.stepping.stepSize.release(Cstep::actor);
    }

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;

    // Return to the propagator
    return;
  }

private:
  /// --------------------------------------------------------------------
  /// Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  ///
  /// @param[in,out] state is the propagation state object
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t>
  void
  initialize(propagator_state_t& state) const
  {
    // Call the navigation helper prior to actual navigation
    debugLog(state, [&] { return std::string("Initialization."); });

    // Set the world volume if it is not set
    if (not state.navigation.worldVolume) {
      state.navigation.worldVolume = trackingGeometry->highestTrackingVolume();
    }

    // We set the current surface to the start surface
    // for eventual post-update action, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    state.navigation.currentSurface = state.navigation.startSurface;
    if (state.navigation.currentSurface) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Current surface set to start surface ";
        dstream << state.navigation.currentSurface->geoID().toString();
        return dstream.str();
      });
    }
    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.navigation.startSurface
        && state.navigation.startSurface->associatedLayer()) {
      debugLog(state, [&] {
        return std::string("Fast start initialization through association.");
      });
      // assign the current layer and volume by association
      state.navigation.startLayer
          = state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume
          = state.navigation.startLayer->trackingVolume();
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
    } else {
      debugLog(state, [&] {
        return std::string("Slow start initialization through search.");
      });
      // current volume and layer search through global search
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Starting from position "
                << toString(state.stepping.position());
        dstream << " and direction " << toString(state.stepping.direction());
        return dstream.str();
      });
      state.navigation.startVolume
          = trackingGeometry->lowestTrackingVolume(state.stepping.position());
      state.navigation.startLayer = state.navigation.startVolume
          ? state.navigation.startVolume->associatedLayer(
                state.stepping.position())
          : nullptr;
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.startVolume) {
        debugLog(state, [&] { return std::string("Start volume resolved."); });
      }
    }
    return;
  }

  /// Status call for test surfaces (surfaces, layers, boundaries)
  /// -------------------------------------------------
  ///
  /// If there are surfaces to be handled, check if the current
  /// state is on the surface
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam navigation_surfaces_t Type of the propagagor
  /// @tparam navigation_iter_t Type of the navigation iterator
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navSurfaces is the navigation status objects
  /// @param[in] surface test surface fore the status test
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t,
            typename navigation_surfaces_t,
            typename navigation_iter_t>
  bool
  status(propagator_state_t&      state,
         navigation_surfaces_t&   navSurfaces,
         const navigation_iter_t& navIter) const
  {
    // No surfaces, status check will be done on layer
    if (navSurfaces.empty() or navIter == navSurfaces.end()) {
      return false;
    }
    // Take the current surface
    auto surface = navIter->representation;
    // Check if we are at a surface
    // If we are on the surface pointed at by the iterator, we can make
    // it the current one to pass it to the other actors
    if (surface->isOnSurface(state.stepping.position(), true)) {
      debugLog(state, [&] {
        return std::string("Status Surface successfully hit, storing it.");
      });
      // Broadcast the surface to the actors and aborters
      state.navigation.currentSurface = surface;
      if (state.navigation.currentSurface) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Current surface set to surface ";
          dstream << state.navigation.currentSurface->geoID().toString();
          return dstream.str();
        });
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
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  bool
  targetSurfaces(propagator_state_t& state, const corrector_t& navCorr) const
  {

    if (state.navigation.navigationBreak) return false;

    // Make sure resolve Surfaces is called on the start layer
    if (state.navigation.startLayer
        and not state.navigation.startLayerResolved) {
      debugLog(state,
               [&] { return std::string("Start layer to be resolved."); });
      // We provide the layer to the resolve surface method in this case
      state.navigation.startLayerResolved = true;
      bool startResolved
          = resolveSurfaces(state, navCorr, state.navigation.startLayer);
      if (not startResolved
          and state.navigation.startLayer == state.navigation.targetLayer) {
        debugLog(state, [&] {
          return std::string("Start is target layer, nothing left to do.");
        });
        // set the navigation break
        state.navigation.navigationBreak = true;
        state.stepping.stepSize.release(Cstep::actor);
      }
      return startResolved;
    }

    // The call that we are on a layer and have not yet resolved the surfaces
    // No surfaces, do not return to stepper
    if (state.navigation.navSurfaces.empty()
        or state.navigation.navSurfaceIter
            == state.navigation.navSurfaces.end()) {
      debugLog(state, [&] {
        return std::string("No surfaces present, target at layer first.");
      });
      return false;
    }

    // Create the navigaton options, the surfaces have initially be ordered
    // so, to catch overstepping anyDirection is allowed here
    NavigationOptions<Surface> navOpts(anyDirection, true);
    navOpts.pathLimit = state.stepping.stepSize.value(Cstep::aborter);

    // Loop over the navigation surfaces
    while (state.navigation.navSurfaceIter
           != state.navigation.navSurfaces.end()) {
      // Screen output how much is left to try
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << std::distance(state.navigation.navSurfaceIter,
                                 state.navigation.navSurfaces.end());
        dstream << " out of " << state.navigation.navSurfaces.size();
        dstream << " surfaces remain to try.";
        return dstream.str();
      });
      // Take the surface
      auto surface = state.navigation.navSurfaceIter->object;
      // Screen output which surface you are on
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Next surface candidate will be ";
        dstream << surface->geoID().toString();
        return dstream.str();
      });
      // Now intersect (should exclude punch-through)
      auto surfaceIntersect
          = surface->intersectionEstimate(state.stepping, navOpts, navCorr);
      double surfaceDistance = surfaceIntersect.intersection.pathLength;
      if (!surfaceIntersect) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Surface intersection at path length ";
          dstream << surfaceDistance;
          dstream << " is not valid.";
          return dstream.str();
        });
        ++state.navigation.navSurfaceIter;
        continue;
      }
      // Update the step for the stepper
      updateStep(state, navCorr, surfaceDistance, true);
      // Return to the propagator
      return true;
    }

    // Reached the end of the surface iteration
    if (state.navigation.navSurfaceIter == state.navigation.navSurfaces.end()) {
      debugLog(state, [&] {
        return std::string("Last surface on layer reached, switching layer.");
      });
      // first clear the surface cache
      state.navigation.navSurfaces.clear();
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
      // now switch to the next layer
      ++state.navigation.navLayerIter;
    }
    // Do not return to the propagator
    return false;
  }

  /// @brief Target layer candidates.
  ///
  /// We are now trying to advance to the next layer (with surfaces)
  /// Check if we are on the representing surface of the layer pointed
  /// at by navLayerIter. If so, we unpack the compatible surfaces
  /// (determined by straight line intersect), and set up the iterator
  /// so that the next call to operator() will enter the surface
  /// check mode above. If no surfaces are found, we skip the layer.
  /// If we unpack a surface, the step size is set to the path length
  /// to the first surface, as determined by straight line intersect.
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  bool
  targetLayers(propagator_state_t& state, const corrector_t& navCorr) const
  {

    if (state.navigation.navigationBreak) return false;

    // if there are no layers, go back to the navigagor (not stepper yet)
    if (state.navigation.navLayers.empty()) {
      debugLog(state, [&] {
        return std::string("No layers present, resolve volume first.");
      });
      if (resolveLayers(state, navCorr)) {
        // The layer resolving worked
        return true;
      }
    }
    // loop over the available navigation layer candiates
    while (state.navigation.navLayerIter != state.navigation.navLayers.end()) {
      // The layer surface
      auto layerSurface = state.navigation.navLayerIter->representation;
      // We are on the layer
      if (state.navigation.currentSurface == layerSurface) {
        debugLog(state, [&] {
          return std::string("We are on a layer, resolve Surfaces.");
        });
        // If you found surfaces return to the propagator
        if (resolveSurfaces(state, navCorr)) {
          return true;
        } else {
          // Try the next one
          ++state.navigation.navLayerIter;
          continue;
        }
      }
      // Otherwise try to step towards it
      NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
      auto layerIntersect = layerSurface->intersectionEstimate(
          state.stepping, navOpts, navCorr);
      // check if the intersect is invalid
      if (!layerIntersect) {
        debugLog(state, [&] {
          return std::string("Layer intersection not valid, skipping it.");
        });
        ++state.navigation.navLayerIter;
      } else {
        // update the navigation step size, release the former first
        updateStep(
            state, navCorr, layerIntersect.intersection.pathLength, true);
        return true;
      }
    }
    debugLog(state, [&] {
      std::stringstream dstream;
      dstream << "Last layer";
      if (state.navigation.currentVolume == state.navigation.targetVolume) {
        dstream << " (final volume) done, proceed to target.";
      } else {
        dstream << " done, target volume boundary.";
      }
      return dstream.str();
    });
    // Set the navigation break if necessary
    state.navigation.navigationBreak
        = (state.navigation.currentVolume == state.navigation.targetVolume);
    return false;
  }

  /// Navigation through volumes
  /// -------------------------------------------------
  ///
  /// This is the boundary check routine. If the code above set up the
  /// boundary surface iterator, we advance through them here. If we are on
  /// the boundary surface, we set the current surface to the boundary
  /// surface, and get the volume pointed at by the boundary surface.  Next
  /// we unpack the layers from that volume. If the volume contains layers
  /// we set the step size to the straight line path length to the first
  /// layer.  If we don't find a next volume, the navigationBreak
  /// indicator is set.  This ends the navigation. Finally, the boundary
  /// iterator is cleared, so that the subsequent call goes back to
  /// the layer iteration logic.
  ///
  /// If we are not on the current boundary surface, we try the next one.
  /// The iterator is advanced and the step size is set. If no straight
  /// line intersect is found, the boundary surface is skipped.
  /// If we are out of boundary surfaces, the navigation is terminated.
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  bool
  targetBoundaries(propagator_state_t& state, const corrector_t& navCorr) const
  {
    if (state.navigation.navigationBreak) return false;

    if (!state.navigation.currentVolume) {
      debugLog(state, [&] {
        return std::string("No sufficient information to resolve boundary, "
                           "stopping navigation.");
      });
      state.stepping.stepSize.release(Cstep::actor);
      return false;
    } else if (state.navigation.currentVolume
               == state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string("In target volume: no need to resolve boundary, "
                           "stopping navigation.");
        state.navigation.navigationBreak = true;
        state.stepping.stepSize.release(Cstep::actor);
      });
      return true;
    }

    // Re-initialize target, because it's good to do so
    initializeTarget(state, navCorr);

    // The navigation options
    NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
    navOpts.pathLimit = state.stepping.stepSize.value(Cstep::aborter);

    // If you came until here, and you might not have boundaries
    // per definition, this is free of the self call
    if (state.navigation.navBoundaries.empty()) {
      // create the navigation options - we could give the start surface here
      NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
      // get the navigation boundaries
      state.navigation.navBoundaries
          = state.navigation.currentVolume->compatibleBoundaries(
              state.stepping, navOpts, navCorr);
      // The number of boundary candidates
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navBoundaries.size();
        dstream << " boundary candidates found at path(s): ";
        for (auto& bc : state.navigation.navBoundaries) {
          dstream << bc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // Set the begin iterator
      state.navigation.navBoundaryIter = state.navigation.navBoundaries.begin();
    }

    // Loop over the boundary surface
    while (state.navigation.navBoundaryIter
           != state.navigation.navBoundaries.end()) {
      // That is the current boundary surface
      auto boundarySurface = state.navigation.navBoundaryIter->representation;
      // Step towards the boundary surface
      auto boundaryIntersect = boundarySurface->intersectionEstimate(
          state.stepping, navOpts, navCorr);
      // Distance
      auto distance = boundaryIntersect.intersection.pathLength;
      // Check the boundary is properly intersected: we are in target mode
      if (!boundaryIntersect
          or distance * distance
              < s_onSurfaceTolerance * s_onSurfaceTolerance) {
        debugLog(state, [&] {
          return std::string("Boundary intersection not valid, skipping it.");
        });
        // get the actual boundary for the navigation & the next volume
        auto boundary = state.navigation.navBoundaryIter->object;
        state.navigation.currentVolume
            = boundary->attachedVolume(state.stepping.position(),
                                       state.stepping.direction(),
                                       state.stepping.navDir);
        // no volume anymore : end of known world
        if (!state.navigation.currentVolume) {
          debugLog(state, [&] {
            return std::string(
                "No more volume to progress to, stopping navigation.");
          });
          return false;
        }
        // store the boundary for eventual actors to work on it
        state.navigation.currentSurface = boundarySurface;
        if (state.navigation.currentSurface) {
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to boundary surface ";
            dstream << state.navigation.currentSurface->geoID().toString();
            return dstream.str();
          });
        }
        // resolve the new layer situation
        if (resolveLayers(state, navCorr)) {
          // positive layer resolving :
          // - we can invalidate the boundary surfaces and return
          state.navigation.navBoundaries.clear();
          state.navigation.navBoundaryIter
              = state.navigation.navBoundaries.end();
          // return to the stepper
          return true;
        }
        // return
        debugLog(state, [&] {
          return std::string("No layers can be reached in the new volume.");
        });
        // create the navigation options - we could give the start surface here
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
        navOpts.startObject = boundarySurface;  // exclude the current boundary
        // re-evaluate the boundary surfaces
        state.navigation.navBoundaries
            = state.navigation.currentVolume->compatibleBoundaries(
                state.stepping, navOpts, navCorr);
        state.navigation.navBoundaryIter
            = state.navigation.navBoundaries.begin();
      }
      // (re-)evaluate the distance to the boundary
      if (state.navigation.navBoundaryIter
          != state.navigation.navBoundaries.end()) {
        // we are not on the layer
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << std::distance(state.navigation.navBoundaryIter,
                                   state.navigation.navBoundaries.end());
          dstream << " out of " << state.navigation.navBoundaries.size();
          dstream << " boundaries remain to try.";
          return dstream.str();
        });

        // update the boundary Surface (could have been switched)
        boundarySurface = state.navigation.navBoundaryIter->representation;
        // intersection with next layer
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
        auto boundaryIntersect = boundarySurface->intersectionEstimate(
            state.stepping, navOpts, navCorr);
        // check if the intersect is invalid
        auto step = boundaryIntersect.intersection.pathLength;
        if (boundaryIntersect) {
          if (step * step < s_onSurfaceTolerance * s_onSurfaceTolerance) {
            // very unlikely edge-case
            lastBoundary = boundarySurface;
            continue;
          } else {
            // Update the navigation step size
            // This is a new navigation stream, release the former first
            updateStep(state, boundaryIntersect.intersection.pathLength, true);
            return true;
          }
        } else {
          // switch to new boundary and re-enter
          ++state.navigation.navBoundaryIter;
          continue;
        }
      }
    }
    // Could not do anything
    return false;
  }

  /// --------------------------------------------------------------------
  /// Navigation (re-)initialisation for the target
  ///
  /// This is only called a few times every propagation/extrapolation
  ///
  /// ---------------------------------------------------------------------
  ///
  /// As a straight line estimate can lead you to the wrong destination
  /// Volume, this will be called at:
  /// - initialization
  /// - attempted volume switch
  /// Target findin by association will not be done again
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  void
  initializeTarget(propagator_state_t& state, const corrector_t& navCorr) const
  {
    debugLog(state, [&] {
      return std::string("We do not have any layers yet, searching.");
    });
    // check if we are in the start volume
    bool start
        = (state.navigation.currentVolume == state.navigation.startVolume);
    // we do not have layers yet, get the candidates
    auto startLayer = start ? state.navigation.startLayer : nullptr;
    auto endLayer   = nullptr;  // state.navigation.targetLayer;
    // create the navigation options - and get the compatible layers
    NavigationOptions<Layer> navOpts(state.stepping.navDir,
                                     true,
                                     resolveSensitive,
                                     resolveMaterial,
                                     resolvePassive,
                                     startLayer,
                                     endLayer);
    // get a corrector associated with the stepper
    state.navigation.navLayers
        = state.navigation.currentVolume->compatibleLayers(
            state.stepping, navOpts, navCorr);

    // the number of layer candidates
    if (!state.navigation.navLayers.empty()) {
      debugLog(state, [&] {
        return std::string("Re-initialzing cancelled as it is the first step.");
      });
      return;
    }
    // Fast Navigation initialization for target:
    if (state.navigation.targetSurface
        && state.navigation.targetSurface->associatedLayer()
        && !state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string("Fast target initialization through association.");
      });
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Target surface set to ";
        dstream << state.navigation.targetSurface->geoID().toString();
        return dstream.str();
      });
      // assign the target volume and the target surface
      state.navigation.targetLayer
          = state.navigation.targetSurface->associatedLayer();
      state.navigation.targetVolume
          = state.navigation.targetLayer->trackingVolume();
    } else if (state.navigation.targetSurface) {
      // screen output that you do a re-initialization
      if (state.navigation.targetVolume) {
        debugLog(state, [&] {
          return std::string("Re-initialization of target volume triggered.");
        });
      }
      // Slow navigation initialization for target:
      // target volume and layer search through global search
      NavigationOptions<Surface> navOpts(state.stepping.navDir,
                                         false,
                                         resolveSensitive,
                                         resolveMaterial,
                                         resolvePassive);
      // take the target intersection
      auto targetIntersection
          = state.navigation.targetSurface->intersectionEstimate(
              state.stepping, navOpts, navCorr);
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Target estimate position (";
        dstream << targetIntersection.intersection.position.x() << ", ";
        dstream << targetIntersection.intersection.position.y() << ", ";
        dstream << targetIntersection.intersection.position.z() << ")";
        return dstream.str();
      });
      /// get the target volume from the intersection
      auto tPosition = targetIntersection.intersection.position;
      state.navigation.targetVolume
          = trackingGeometry->lowestTrackingVolume(tPosition);
      state.navigation.targetLayer = state.navigation.targetVolume
          ? state.navigation.targetVolume->associatedLayer(tPosition)
          : nullptr;
      if (state.navigation.targetVolume) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Target volume estimated : ";
          dstream << state.navigation.targetVolume->volumeName();
          return dstream.str();
        });
      }
    }
  }

  /// @brief Resolve the surfaces of this layer
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  bool
  resolveSurfaces(propagator_state_t& state,
                  const corrector_t&  navCorr,
                  const Layer*        cLayer = nullptr) const
  {
    // get the layer and layer surface
    auto layerSurface = cLayer ? state.navigation.startSurface
                               : state.navigation.navLayerIter->representation;
    auto navLayer = cLayer ? cLayer : state.navigation.navLayerIter->object;
    // are we on the start layer
    bool onStart      = (navLayer == state.navigation.startLayer);
    auto startSurface = onStart ? state.navigation.startSurface : layerSurface;
    NavigationOptions<Surface> navOpts(state.stepping.navDir,
                                       true,
                                       resolveSensitive,
                                       resolveMaterial,
                                       resolvePassive,
                                       startSurface,
                                       state.navigation.targetSurface);
    // Check the limit
    navOpts.pathLimit = state.stepping.stepSize.value(Cstep::aborter);
    // get the surfaces
    state.navigation.navSurfaces
        = navLayer->compatibleSurfaces(state.stepping, navOpts, navCorr);
    // the number of layer candidates
    if (!state.navigation.navSurfaces.empty()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navSurfaces.size();
        dstream << " surface candidates found at path(s): ";
        for (auto& sfc : state.navigation.navSurfaces) {
          dstream << sfc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // set the iterator
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      // Update the navigation step size before you return to the stepper
      // This is a new navigation stream, release the former step size first
      updateStep(state,
                 navCorr,
                 state.navigation.navSurfaceIter->intersection.pathLength,
                 true);
      return true;
    }
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
    debugLog(state,
             [&] { return std::string("No surface candidates found."); });
    return false;
  }

  /// Navigation through layers
  /// -------------------------------------------------
  ///
  /// Resolve layers.
  ///
  /// This initializes the layer candidates when starting
  /// or when entering a new volume
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param[in,out] state is the propagation state object
  /// @param[in] navCorr is the navigation path corrector
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename corrector_t>
  bool
  resolveLayers(propagator_state_t& state, const corrector_t& navCorr) const
  {
    debugLog(state,
             [&] { return std::string("Searching for compatible layers."); });

    // Check if we are in the start volume
    auto startLayer
        = (state.navigation.currentVolume == state.navigation.startVolume)
        ? state.navigation.startLayer
        : nullptr;
    // Create the navigation options
    // - and get the compatible layers, start layer will be excluded
    NavigationOptions<Layer> navOpts(state.stepping.navDir,
                                     true,
                                     resolveSensitive,
                                     resolveMaterial,
                                     resolvePassive,
                                     startLayer,
                                     nullptr);
    // Set also the target surface
    navOpts.targetSurface = state.navigation.targetSurface;
    navOpts.pathLimit     = state.stepping.stepSize.value(Cstep::aborter);
    // Request the compatible layers
    state.navigation.navLayers
        = state.navigation.currentVolume->compatibleLayers(
            state.stepping, navOpts, navCorr);

    // Layer candidates have been found
    if (!state.navigation.navLayers.empty()) {
      // Screen output where they are
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navLayers.size();
        dstream << " layer candidates found at path(s): ";
        for (auto& lc : state.navigation.navLayers) {
          dstream << lc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // Set the iterator to the first
      state.navigation.navLayerIter = state.navigation.navLayers.begin();
      // Setting the step size towards first
      if (state.navigation.startLayer
          && state.navigation.navLayerIter->object
              != state.navigation.startLayer) {
        debugLog(state, [&] { return std::string("Target at layer."); });
        // update the navigation step size before you return
        updateStep(state,
                   navCorr,
                   state.navigation.navLayerIter->intersection.pathLength,
                   true);
        // Trigger the return to the propagator
        return true;
      }
    }

    // Screen output - no layer candidates found
    debugLog(state, [&] {
      return std::string("No compatible layer candidates found.");
    });
    // Update the navigation step to the target step to trigger
    // step modification when requested
    updateStep(state, navCorr, state.stepping.stepSize.value(Cstep::aborter));

    return false;
  }

  /// This method updates the constrained step size
  ///
  /// @tparam propagator_state_t is the state type
  ///
  /// @param[in,out] state The state object for the step length
  /// @param[in] step the step size
  /// @param[in] release flag steers if the step is released first
  template <typename propagator_state_t, typename corrector_t>
  void
  updateStep(propagator_state_t& state,
             const corrector_t&  navCorr,
             double              navigationStep,
             bool                release = false) const
  {  //  update the step
    state.stepping.stepSize.update(navigationStep, Cstep::actor, release);
    debugLog(state, [&] {
      std::stringstream dstream;
      std::string       releaseFlag = release ? "released and " : "";
      dstream << "Navigation stepSize " << releaseFlag << "updated to ";
      dstream << state.stepping.stepSize.toString();
      return dstream.str();
    });
    /// If we have an initial step and are configured to modify it
    if (state.stepping.pathAccumulated == 0.
        and navCorr(state.stepping.stepSize)) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Initial navigation step modified to ";
        dstream << state.stepping.stepSize.toString();
        return dstream.str();
      });
    }
  }

  /// --------------------------------------------------------------------
  /// Inactive
  ///
  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  ///
  /// @param[in,out] state is the propagation state object
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t>
  bool
  inactive(propagator_state_t& state) const
  {
    // Void behavior in case no tracking geometry is present
    if (!trackingGeometry) {
      return true;
    }
    // turn the navigator into void when you are intructed to do nothing
    if (!resolveSensitive && !resolveMaterial && !resolvePassive) {
      return true;
    }

    // --------------------------------------------------------------------
    // Navigation break handling
    // This checks if a navigation break had been triggered:
    // - If so & the target exists or was hit - it simply returns
    // - If a target exists and was not yet hit, it checks for it
    // -> return is always to the stepper
    return navigationBreak(state);
  }

  /// --------------------------------------------------------------------
  /// Navigation break handling
  ///
  /// This checks if a navigation break had been triggered
  /// - If so & the target exists or was hit - it simply returns
  /// - If a target exists and was not yet hit, it checks for it
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  ///
  /// @param[in,out] state is the propagation state object
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t>
  bool
  navigationBreak(propagator_state_t& state) const
  {
    if (state.navigation.navigationBreak) {
      // target exists and reached, or no target exists
      if (state.navigation.targetReached || !state.navigation.targetSurface) {
        return true;
      }
      // the only advance could have been to the target
      if (state.navigation.targetSurface->isOnSurface(state.stepping.position(),
                                                      true)) {
        // set the target surface
        state.navigation.currentSurface = state.navigation.targetSurface;

        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface ";
          dstream << state.navigation.currentSurface->geoID().toString();
          return dstream.str();
        });
        return true;
      }
    }
    return false;
  }

  /// The private navigation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// state.options.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param[in,out] state the propagator state for the debug flag,
  ///      prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
    if (state.options.debug) {
      std::string vName = "No Volume";
      if (state.navigation.currentVolume) {
        vName = state.navigation.currentVolume->volumeName();
      }
      std::vector<std::string> lines;
      std::string              input = logAction();
      boost::split(lines, input, boost::is_any_of("\n"));
      for (const auto& line : lines) {
        std::stringstream dstream;
        dstream << ">>>" << std::setw(state.options.debugPfxWidth) << vName
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
        state.options.debugString += dstream.str();
      }
    }
  }
};

}  // namespace Acts
