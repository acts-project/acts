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

namespace Acts {

typedef detail::ConstrainedStep cstep;

/// @breif struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions
{

  /// the navigation direction
  NavigationDirection navDir = forward;

  /// the boundary check directive
  BoundaryCheck boundaryCheck = true;

  // how to resolve the geometry
  /// always look for sensitive
  bool resolveSensitive = true;
  /// always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// object to check against: at start
  const object_t* startObject = nullptr;
  /// object to check against: at end
  const object_t* endObject = nullptr;

  double pathLimit = std::numeric_limits<double>::max();

  /// Constructor
  ///
  /// @param nDir Navigation direction prescription
  /// @param bcheck Boundary check for the navigaiton action
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
    , boundaryCheck(bcheck)
    , resolveSensitive(resolves)
    , resolveMaterial(resolvem)
    , resolvePassive(resolvep)
    , startObject(sobject)
    , endObject(eobject)
    , pathLimit(ndir * std::numeric_limits<double>::max())
  {
  }
};

typedef std::vector<SurfaceIntersection> NavigationSurfaces;
typedef NavigationSurfaces::iterator     NavigationSurfaceIter;

typedef std::vector<LayerIntersection> NavigationLayers;
typedef NavigationLayers::iterator     NavigationLayerIter;

typedef std::vector<BoundaryIntersection> NavigationBoundaries;
typedef NavigationBoundaries::iterator    NavigationBoundaryIter;

/// Navigator struct
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
struct Navigator
{

  /// Constructor with shared tracking geometry
  ///
  /// @param tGeometry The tracking geometry for the navigator
  Navigator(std::shared_ptr<const TrackingGeometry> tGeometry = nullptr)
    : trackingGeometry(tGeometry)
  {
  }

  /// Tracking Geometry for this Navigator
  std::shared_ptr<const TrackingGeometry> trackingGeometry;

  /// the tolerance used to defined "reached"
  double tolerance = s_onSurfaceTolerance;

  /// Change initial step (avoids overstepping, when no correction done)
  double initialStepFactor = 0.5;

  /// Number of target attempts for the targetVolume
  int targetAttempts = 2;

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
    /// Indicator if the target is reached
    bool targetReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;
  };

  /// Unique typedef to publish to the Propagator
  typedef State state_type;

  /// Navigation action for the ActionList of the Propagator
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  ///
  /// @param[in,out] state is the mutable propagator state object
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state) const
  {

    // void behavior in case no tracking geometry is present
    if (!trackingGeometry) return;

    // turn the navigator into void when you are intructed to do nothing
    if (!resolveSensitive && !resolveMaterial && !resolvePassive) return;

    debugLog(state, [&] { return std::string("Entering navigator."); });

    // Navigator always resets the current surface first
    state.navigation.currentSurface = nullptr;

    // --------------------------------------------------------------------
    // Navigation break handling
    // This checks if a navigation break had been triggered:
    // - If so & the target exists or was hit - it simply returns
    // - If a target exists and was not yet hit, it checks for it
    // -> return is always to the stepper
    if (navigationBreak(state)) return;

    // get a  navigation corrector associated to the stepper
    auto navCorr = state.stepping.corrector();

    // -------------------------------------------------
    // Initialization
    // This should lead to:
    // - a current volume
    // - potentially also a current layer
    // -> return is always to the stepper
    if (initialize(state, navCorr)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from initialize.");
      });
      return;
    }

    // -------------------------------------------------
    // Surfaces (if present)
    // - this can only happen after a layer has  sucessfully been resolved
    // -> return is always to the stepper
    if (handleSurfaces(state, navCorr)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from surface handling.");
      });
      return;
    }

    // -------------------------------------------------
    // Layers are present
    // - this can only happen after a volume has successfully been resolved
    // -> return is always to the stepper
    if (handleLayers(state, navCorr)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from layer handling.");
      });
      return;
    }
    // -------------------------------------------------
    // Volume be handled
    // - if you arrived
    // -> return is always to the stepper
    if (handleBoundaries(state, navCorr)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from boundary handling.");
      });
      return;
    }
    // -------------------------------------------------
    // neither surfaces, layers nor boundaries triggered a return
    // navigation broken - switch navigator off
    state.navigation.navigationBreak = true;
    debugLog(state, [&] {
      return std::string("Navigation break - no valid actions left.");
    });
    // release the navigation step size
    state.stepping.stepSize.release(cstep::actor);
    return;
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
      if (state.navigation.targetReached || !state.navigation.targetSurface)
        return true;
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

  /// --------------------------------------------------------------------
  /// Navigation initialisation
  ///
  /// This is only called once for every propagation/extrapolation
  ///
  /// ---------------------------------------------------------------------
  /// Check for navigation initialisation & do it if necessary.  This means
  /// we do not have an active volume yet (since we just started).  We get
  /// the innermost volume, and set up an ordered layer iterator and step
  /// size toward the first layer.  Return to the stepper prevent execution
  /// of the subsequent logic, i.e. we want to take a step first.
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
  initialize(propagator_state_t& state, const corrector_t& navCorr) const
  {

    // No initialisation necessary
    if (state.navigation.currentVolume) return false;

    debugLog(state, [&] { return std::string("Initializing start volume."); });

    // we set the current surface to the start surface
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
    // short-cut through object association, saves navigation in the
    // geometry and volume tree search for the lowest volume
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
      // set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.startVolume) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Start volume estimated : ";
          dstream << state.navigation.startVolume->volumeName();
          return dstream.str();
        });
      }
    }

    // Target initialization
    initializeTarget(state, navCorr);

    // A current volume exists
    if (state.navigation.startVolume) {
      // assign to the currentVolume
      state.navigation.currentVolume = state.navigation.startVolume;
      // fast exit if start and target layer are identical
      if (state.navigation.startLayer == state.navigation.targetLayer) {
        debugLog(state, [&] {
          return std::string(
              "Start and target layer identical, check surfaces.");
        });
        // resolve the surfaces only
        NavigationOptions<Surface> navOpts(state.stepping.navDir,
                                           true,
                                           resolveSensitive,
                                           resolveMaterial,
                                           resolvePassive,
                                           state.navigation.startSurface,
                                           state.navigation.targetSurface);
        state.navigation.navSurfaces
            = state.navigation.startLayer->compatibleSurfaces(
                state.stepping, navOpts, navCorr);

        // the number of surface candidates
        if (state.navigation.navSurfaces.size()) {
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
          state.navigation.navSurfaceIter
              = state.navigation.navSurfaces.begin();
          // update the navigation step size before you return
          updateStep(state,
                     state.navigation.navSurfaceIter->intersection.pathLength);
          return true;
        }
        return false;
      }
      // initialize layer - if it works go ahead with processing
      if (resolveLayers(state, navCorr)) {
        if (state.stepping.stepSize == 0.) {
          debugLog(state, [&] {
            return std::string("On current layer surface, setting it.");
          });
          state.navigation.currentSurface
              = state.navigation.navLayerIter->representation;
          if (state.navigation.currentSurface)
            debugLog(state, [&] {
              std::stringstream dstream;
              dstream << "Current surface set to approach surface ";
              dstream << state.navigation.currentSurface->geoID().toString();
              return dstream.str();
            });
          // no returning to the stepper at this stage
          return false;
        }
        // we have a step size different from 0, return to the stepper
        return true;
      }
      // no layers to be resolved, out of this, but not back to stepper (yet)
      return false;
    }
    // navigation broken - switch navigator off
    state.navigation.navigationBreak = true;
    debugLog(state, [&] {
      return std::string("Navigation broken, pure propagation.");
    });
    // release the navigation step size
    state.stepping.stepSize.release(cstep::actor);
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

    if (state.navigation.targetVolume && state.stepping.pathAccumulated == 0.) {
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
  handleBoundaries(propagator_state_t& state, const corrector_t& navCorr) const
  {

    // re-initialize target: @todo we know better when to do this
    initializeTarget(state, navCorr);

    // only handle boundaries if you are not in the target volume
    if (!state.navigation.currentVolume) {
      debugLog(state, [&] {
        return std::string("No sufficient information to resolve boundary, "
                           "stopping navigation.");
      });
      return false;
    } else if (state.navigation.currentVolume
               == state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string("In target volume: no need to resolve boundary, "
                           "stopping navigation.");
      });
      return false;
    }
    // remember if you are on a boundary - needed for self call w/o layers
    const Surface* lastBoundary = nullptr;

    // if you came until here, and you might not have boundaries
    // per definition, this is free of the self call
    if (!state.navigation.navBoundaries.size()) {
      // create the navigaiton options - we could give the start surface here
      NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
      // get the navigation boundaries
      state.navigation.navBoundaries
          = state.navigation.currentVolume->compatibleBoundaries(
              state.stepping, navOpts, navCorr);
      // the number of boundary candidates
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navBoundaries.size();
        dstream << " boundary candidates found at path(s): ";
        for (auto& bc : state.navigation.navBoundaries)
          dstream << bc.intersection.pathLength << "  ";
        return dstream.str();
      });
      // set the iterator, but avoid stepping a zero step
      if (state.navigation.navBoundaries.size()) {
        state.navigation.navBoundaryIter
            = state.navigation.navBoundaries.begin();
        auto step = state.navigation.navBoundaryIter->intersection.pathLength;
        if (step * step < s_onSurfaceTolerance * s_onSurfaceTolerance) {
          lastBoundary = state.navigation.navBoundaryIter->representation;
          debugLog(state, [&] {
            return std::string("Starting from boundary surface.");
          });
        } else {
          // update the navigation step size before you return
          updateStep(state,
                     state.navigation.navBoundaryIter->intersection.pathLength);
          return true;
        }
      } else {
        debugLog(state, [&] {
          return std::string(
              "No valid boundary surface found, stopping navigation.");
        });
        return false;
      }
    }

    // loop over rest of the boundaries
    while (state.navigation.navBoundaryIter
           != state.navigation.navBoundaries.end()) {
      auto boundarySurface = state.navigation.navBoundaryIter->representation;
      // case (v-a) : you are on the boundary surface
      //              only if you hadn't just done a volume switch
      // check if we are on already in this step
      if (lastBoundary
          || boundarySurface->isOnSurface(state.stepping.position(), true)) {
        debugLog(state, [&] {
          return std::string(
              "Boundary surface reached, prepare volume switch.");
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
        if (state.navigation.currentSurface)
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to boundary surface ";
            dstream << state.navigation.currentSurface->geoID().toString();
            return dstream.str();
          });
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
        // create the navigaiton options - we could give the start surface here
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
            // update the navigation step size
            updateStep(state, boundaryIntersect.intersection.pathLength);
            return true;
          }
        } else {
          // switch to new boundary and re-enter
          ++state.navigation.navBoundaryIter;
          continue;
        }
      }
    }
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
    debugLog(state, [&] {
      return std::string("We do not have any layers yet, searching.");
    });
    // check if we are in the start volume
    bool start
        = (state.navigation.currentVolume == state.navigation.startVolume);
    // we do not have layers yet, get the candidates
    auto startLayer = start ? state.navigation.startLayer : nullptr;
    auto endLayer   = nullptr;  // state.navigation.targetLayer;
    // create the navigaiton options - and get the compatible layers
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
    if (state.navigation.navLayers.size()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navLayers.size();
        dstream << " layer candidates found at path(s): ";
        for (auto& lc : state.navigation.navLayers)
          dstream << lc.intersection.pathLength << "  ";
        return dstream.str();
      });
      // set the iterator
      state.navigation.navLayerIter = state.navigation.navLayers.begin();

      if (state.navigation.startLayer
          && state.navigation.navLayerIter->object
              != state.navigation.startLayer) {
        debugLog(state,
                 [&] { return std::string("Stepping towards first layer."); });
        // update the navigation step size before you return
        updateStep(state,
                   state.navigation.navLayerIter->intersection.pathLength);
        return true;
      } else {
        debugLog(state, [&] {
          return std::string(
              "Start layer, avoid step to layer approach surface.");
        });
        // we can call resolveSurfaces here, but first let's set the
        debugLog(state, [&] {
          return std::string("Switch layer surface to start surface.");
        });
        // setting layer surface & representation
        state.navigation.navLayerIter->representation
            = state.navigation.startSurface;
        // if you found surfaces return to the stepper
        if (resolveSurfaces(state, navCorr)) return true;
        // increase the iterator
        ++state.navigation.navLayerIter;
        return false;
      }
    }

    // your estimate about the target may have been changed
    initializeTarget(state, navCorr);

    if (state.navigation.currentVolume != state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string("No layer candidates found, switching volume.");
      });
      return false;
    }
    debugLog(state, [&] {
      return std::string(
          "Done in final volume, release stepSize & proceed to target.");
    });
    // the step size will be set to the aborter step size
    double abortStep = state.stepping.stepSize.value(cstep::aborter);
    updateStep(state, abortStep);
    state.navigation.navigationBreak = true;
    return true;
  }

  /// @brief Loop over layer candidates.
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
  handleLayers(propagator_state_t& state, const corrector_t& navCorr) const
  {
    // if there are no layers, go back to the navigagor (not stepper yet)
    if (state.navigation.navLayers.empty()) {
      debugLog(state, [&] {
        return std::string("No layers present, resolve volume first.");
      });
      return false;
    }

    // loop over the available navigation layer candiates
    while (state.navigation.navLayerIter != state.navigation.navLayers.end()) {
      // we take the layer representation surface
      auto layer        = state.navigation.navLayerIter->object;
      auto layerSurface = state.navigation.navLayerIter->representation;
      auto layerVolume  = layer->representingVolume();
      // check if we are on the layer
      // - fastest: is start layer
      bool onLayer = (layer == state.navigation.startLayer);
      // @todo: check if another fast association can be used
      if (!onLayer) {
        // - not so fast : inside volume
        onLayer = layerVolume && layerVolume->inside(state.stepping.position());
        // - not so fast: on layer surface
        onLayer
            = (onLayer
               || layerSurface->isOnSurface(state.stepping.position(), true));
      }
      // if we are on the layer alreay at this step
      // - we can resolve the compatible surfaces
      if (onLayer) {
        debugLog(state, [&] {
          return std::string("On the layer, resolving compatible surfaces.");
        });
        // this is the case where you actually start from a surface within the
        // layer : store the current surface in the propagator state
        if (state.navigation.startSurface
            && state.navigation.currentVolume == state.navigation.startVolume
            && layer == state.navigation.startLayer) {
          debugLog(state, [&] {
            return std::string("Switch layer surface to start surface.");
          });
          // setting layer surface & representation
          state.navigation.navLayerIter->representation
              = state.navigation.startSurface;
        } else {
          state.navigation.currentSurface = layerSurface;
          if (state.navigation.currentSurface)
            debugLog(state, [&] {
              std::stringstream dstream;
              dstream << "Current surface set to layer surface ";
              dstream << state.navigation.currentSurface->geoID().toString();
              return dstream.str();
            });
        }
        // if you found surfaces return to the stepper
        if (resolveSurfaces(state, navCorr)) return true;
        // increase the iterator
        ++state.navigation.navLayerIter;
      }
      // we either re-evaluate the step to the current layer
      // or evaluate the step to the next layer
      if (state.navigation.navLayerIter != state.navigation.navLayers.end()) {
        // update in case a switch was done
        layerSurface = state.navigation.navLayerIter->representation;
        // we are not on the layer
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << std::distance(state.navigation.navLayerIter,
                                   state.navigation.navLayers.end());
          dstream << " out of " << state.navigation.navLayers.size();
          dstream << " layers remain to try.";
          return dstream.str();
        });
        /// intersection with next layer
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
          // update the navigation step size
          updateStep(state, layerIntersect.intersection.pathLength);
          return true;
        }
      }
    }
    // let's re-initialize the target
    initializeTarget(state, navCorr);
    // we are at the end of trying layers
    if (state.navigation.currentVolume == state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string(
            "Done in final volume, release stepSize & proceed to target.");
      });
      // the step size will be set to the aborter step size
      double abortStep = state.stepping.stepSize.value(cstep::aborter);
      updateStep(state, abortStep);
      state.navigation.navigationBreak = true;
      return true;
    }
    debugLog(state, [&] {
      return std::string("All layers been handled, switching volume.");
    });
    // clear the layers
    state.navigation.navLayers.clear();
    state.navigation.navLayerIter = state.navigation.navLayers.end();
    // clear the boundaries
    state.navigation.navBoundaries.clear();
    state.navigation.navBoundaryIter = state.navigation.navBoundaries.end();
    // do not return to the stepper here (yet)
    return false;
  }

  /// @brief Resolve the surfaces of this layer, if not the start layer
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
  resolveSurfaces(propagator_state_t& state, const corrector_t& navCorr) const
  {
    // get the layer and layer surface
    auto layerSurface = state.navigation.navLayerIter->representation;
    auto navLayer     = state.navigation.navLayerIter->object;
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
    // get the surfaces
    state.navigation.navSurfaces
        = navLayer->compatibleSurfaces(state.stepping, navOpts, navCorr);
    // the number of layer candidates
    if (state.navigation.navSurfaces.size()) {
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
      // update the navigation step size before you return to the stepper
      updateStep(state,
                 state.navigation.navSurfaceIter->intersection.pathLength);
      return true;
    }
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
    debugLog(state, [&] {
      return std::string("No surface candidates found, switching layer.");
    });
    return false;
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
  handleSurfaces(propagator_state_t& state, const corrector_t& navCorr) const
  {
    // no surfaces, do not return to stepper
    if (!state.navigation.navSurfaces.size()) {
      debugLog(state, [&] {
        return std::string("No surfaces present, resolve layer first.");
      });
      return false;
    }
    // loop over the navigation surfaces
    while (state.navigation.navSurfaceIter
           != state.navigation.navSurfaces.end()) {
      // take the surface
      auto surface = state.navigation.navSurfaceIter->object;
      // case (s-a): we are at a surface
      //           only ask if you hadn't just come from a valid surface
      //
      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      if (surface->isOnSurface(state.stepping.position(), true)) {
        debugLog(state, [&] {
          return std::string("Surface successfully hit, storing it.");
        });
        // the surface will only appear due to correct
        // collect(Property) flag
        state.navigation.currentSurface = surface;
        if (state.navigation.currentSurface)
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to resolved surface ";
            dstream << state.navigation.currentSurface->geoID().toString();
            return dstream.str();
          });
        // break if the surface is the target surface
        if (surface == state.navigation.targetSurface) {
          debugLog(state, [&] {
            return std::string("This was the target surface. Done.");
          });
          return true;
        }
        // switch to the next candidate
        ++state.navigation.navSurfaceIter;
      }
      // screen output how much is left to try
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << std::distance(state.navigation.navSurfaceIter,
                                 state.navigation.navSurfaces.end());
        dstream << " out of " << state.navigation.navSurfaces.size();
        dstream << " surfaces remain to try.";
        return dstream.str();
      });
      // case (s-b) : update step estimation to the (current/new) surface
      if (state.navigation.navSurfaceIter
          != state.navigation.navSurfaces.end()) {
        surface = state.navigation.navSurfaceIter->object;
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Next surface candidate will be ";
          dstream << surface->geoID().toString();
          return dstream.str();
        });
        // create the navigaton options
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
        // now intersect (should exclude punch-through)
        auto surfaceIntersect
            = surface->intersectionEstimate(state.stepping, navOpts, navCorr);
        double surfaceDistance = surfaceIntersect.intersection.pathLength;

        if (!surfaceIntersect) {
          debugLog(state, [&] {
            return std::string(
                "Surface intersection is not valid, skipping it.");
          });
          ++state.navigation.navSurfaceIter;
          continue;
        } else {
          updateStep(state, surfaceDistance);
          return true;
        }
      }
    }
    // case (s-c) : reached the end of the surface iteration
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
    // do not return to the stepper
    return false;
  }

  /// This method updates the constrained step size
  ///
  /// @tparam propagator_state_t is the state type
  ///
  /// @param[in,out] state The state object for the step length
  /// @param[in] step the step size
  template <typename propagator_state_t>
  void
  updateStep(propagator_state_t& state, double step) const
  {
    /// If we have an initial step and are configured to modify it
    if (state.stepping.pathAccumulated == 0. && initialStepFactor != 1.) {
      debugLog(state, [&] {
        return std::string("Initial step modification detected.");
      });
      step *= initialStepFactor;
    }

    //  update the step
    state.stepping.stepSize.update(step, cstep::actor);
    debugLog(state, [&] {
      std::stringstream dstream;
      dstream << "Navigation stepSize updated to ";
      dstream << state.stepping.stepSize.toString();
      return dstream.str();
    });
  }

private:
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
  debugLog(propagator_state_t&          state,
           std::function<std::string()> logAction) const
  {
    if (state.options.debug) {
      std::string vName = "No Volume";
      if (state.navigation.currentVolume)
        vName = state.navigation.currentVolume->volumeName();
      std::stringstream dstream;
      dstream << ">>>" << std::setw(state.options.debugPfxWidth) << vName
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};
}  // namespace Acts
