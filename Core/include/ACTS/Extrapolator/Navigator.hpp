// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_NAVIGATOR_H
#define ACTS_NAVIGATOR_H

#include <iterator>
#include <sstream>
#include <string>
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Propagator/detail/ConstrainedStep.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"

namespace Acts {

typedef detail::ConstrainedStep cstep;

/// Struct to mimic track parameters
/// @todo harmonize to eventual future update of
/// TrackingVolume::layerCandidatesOrdered()
/// TrackingVolume::boundarySurfacesOrdered()
struct NavigationParameters
{

  /// Position
  Vector3D pos;

  /// Direction
  Vector3D dir;

  /// Access method to satisify TrackingVolume interface
  const Vector3D&
  position() const
  {
    return pos;
  }

  /// Access method to satisify TrackingVolume interface
  const Vector3D&
  momentum() const
  {
    return dir;
  }

  /// Constructor
  NavigationParameters(const Vector3D& p, const Vector3D& d) : pos(p), dir(d) {}
};

typedef std::vector<SurfaceIntersection> NavigationSurfaces;
typedef NavigationSurfaces::iterator     NavigationSurfaceIter;

typedef std::vector<LayerIntersection<NavigationParameters>> NavigationLayers;
typedef NavigationLayers::iterator NavigationLayerIter;

typedef std::vector<BoundaryIntersection<NavigationParameters>>
                                       NavigationBoundaries;
typedef NavigationBoundaries::iterator NavigationBoundaryIter;

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
/// If a surface is found, the state.currentSurface
/// pointer is set. This  enables subsequent actors to react. Secondly, this
/// actor uses the ordered  iterators  to figure out which surface, layer or
/// volume boundary is _supposed_ to be hit next. It then sets the maximum
/// step size to the path length found out by straight line intersection. If
/// the isOnSurface call fails, it also  re-computes the step size, to make
/// sure we end up at the desired surface.
struct Navigator
{

  /// Constructor with shared tracking geometry
  Navigator(std::shared_ptr<const TrackingGeometry> tGeometry = nullptr)
    : trackingGeometry(tGeometry)
  {
  }

  /// Tracking Geometry for this Navigator
  std::shared_ptr<const TrackingGeometry> trackingGeometry;

  /// the tolerance used to defined "reached"
  double tolerance = 1e-5;

  /// step safety for inital approach to layers
  /// @note can be set to 1. for non-homogeneous field
  /// @todo can be omitted with better initial step estimation
  double initialStepFactor = 1.;

  /// the navigation level, see Layer for more information
  int navigationLevel = 2;

  /// Configuration for this Navigator
  /// stop at every sensitive surface (whether it has material or not)
  bool collectSensitive = true;
  /// stop at every material surface (whether it is passive or not)
  bool collectMaterial = true;
  /// stop at every surface regardless what it is
  bool collectPassive = false;

  /// store the debug message
  bool debug = false;

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State
  {
    // Navigation on surface level
    /// the vector of navigation surfaces to work throgh
    NavigationSurfaces navSurfaces = {};
    /// the current surface iterator of the navigation state
    NavigationSurfaceIter navSurfaceIter = navSurfaces.end();

    // Navigation on layer level
    /// the vector of navigation layers to work throgh
    NavigationLayers navLayers = {};
    /// the current layer iterator of the navigation state
    NavigationLayerIter navLayerIter = navLayers.end();

    // Navigation on volume level
    /// the vector of boundary surfaces to work throgh
    NavigationBoundaries navBoundaries = {};
    /// the current boundary iterator of the navigation state
    NavigationBoundaryIter navBoundaryIter = navBoundaries.end();

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;
    /// Navigation state: the target layer
    const Layer* targetLayer = nullptr;

    /// break the navigation
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
    // fail if you have no tracking geometry
    assert(trackingGeometry != nullptr);
    debugLog(state, [&] { return std::string("Entering navigator."); });

    // navigation parameters
    NavigationParameters navPar(state.stepping.position(),
                                state.stepping.navDir
                                    * state.stepping.direction());

    // Navigator always resets the current surface first
    state.currentSurface = nullptr;

    // --------------------------------------------------------------------
    // Navigation break handling
    // This checks if a navigation break had been triggered:
    // - If so & the target exists or was hit - it simply returns
    // - If a target exists and was not yet hit, it checks for it
    // -> return is always to the stepper
    if (navigationBreak(navPar, state)) return;

    // -------------------------------------------------
    // Initialization
    // This should lead to:
    // - a current volume
    // - potentially also a current layer
    // -> return is always to the stepper
    if (initialize(navPar, state)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from initialize.");
      });
      return;
    }

    // -------------------------------------------------
    // Surfaces (if present)
    // - this can only happen after a layer has  sucessfully been resolved
    // -> return is always to the stepper
    if (handleSurfaces(navPar, state)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from surface handling.");
      });
      return;
    }

    // -------------------------------------------------
    // Layers are present
    // - this can only happen after a volume has successfully been resolved
    // -> return is always to the stepper
    if (handleLayers(navPar, state)) {
      debugLog(state, [&] {
        return std::string("Return to stepper - from layer handling.");
      });
      return;
    }
    // -------------------------------------------------
    // Volume be handled
    // - if you arrived
    // -> return is always to the stepper
    if (handleBoundaries(navPar, state)) {
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
      return std::string("Naivgation break - no valid actions left.");
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param[in] navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t>
  bool
  navigationBreak(const NavigationParameters& navPar,
                  propagator_state_t&         state) const
  {
    if (state.navigation.navigationBreak) {
      // target exists and reached, or no target exists
      if (state.targetReached || !state.targetSurface) return true;
      // the only advance could have been to the target
      if (state.targetSurface->isOnSurface(navPar.position(), true)) {
        // set the target surface
        state.currentSurface = state.targetSurface;
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface ";
          dstream << state.currentSurface->geoID().toString();
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param[in] navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  ///
  /// @return boolean trigger if successful
  template <typename propagator_state_t>
  bool
  initialize(const NavigationParameters& navPar,
             propagator_state_t&         state) const
  {

    // no initialisation necessary
    if (state.navigation.currentVolume) return false;

    debugLog(state, [&] { return std::string("Initializing start volume."); });

    // we set the current surface to the start surface
    // for eventual post-update action, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    state.currentSurface = state.startSurface;
    if (state.currentSurface) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Current surface set to start surface ";
        dstream << state.currentSurface->geoID().toString();
        return dstream.str();
      });
    }

    // Fast Navigation initialization for start condition:
    // short-cut through object association, saves navigation in the
    // geometry and volume tree search for the lowest volume
    if (state.startSurface && state.startSurface->associatedLayer()) {
      debugLog(state, [&] {
        return std::string("Fast start initialization through association.");
      });
      // assign the current layer and volume by association
      state.navigation.startLayer = state.startSurface->associatedLayer();
      state.navigation.startVolume
          = state.navigation.startLayer->trackingVolume();
    } else {
      debugLog(state, [&] {
        return std::string("Slow start initialization through search.");
      });

      // current volume and layer search through global search
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Starting from position (" << navPar.position().x();
        dstream << ", " << navPar.position().y();
        dstream << ", " << navPar.position().z() << ")";
        return dstream.str();
      });
      state.navigation.startVolume
          = trackingGeometry->lowestTrackingVolume(navPar.position());
      state.navigation.startLayer = state.navigation.startVolume
          ? state.navigation.startVolume->associatedLayer(navPar.position())
          : nullptr;
    }
    // Fast Navigation initialization for target:
    if (state.targetSurface && state.targetSurface->associatedLayer()) {
      debugLog(state, [&] {
        return std::string("Fast target initialization through association.");
      });
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Target surface set to";
        dstream << state.targetSurface->geoID().toString();
        return dstream.str();
      });
      // assign the target volume and the target surface
      state.navigation.targetLayer = state.targetSurface->associatedLayer();
      state.navigation.targetVolume
          = state.navigation.targetLayer->trackingVolume();
    } else if (state.targetSurface) {
      // Slow navigation initialization for target:
      // target volume and layer search through global search
      auto targetIntersection = state.targetSurface->intersectionEstimate(
          navPar.position(), navPar.momentum(), true, false);
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Target estimate position (";
        dstream << targetIntersection.position.x() << ", ";
        dstream << targetIntersection.position.y() << ", ";
        dstream << targetIntersection.position.z() << ")";
        return dstream.str();
      });
      /// get the target volume from the intersection
      state.navigation.targetVolume
          = trackingGeometry->lowestTrackingVolume(targetIntersection.position);
      state.navigation.targetLayer = state.navigation.targetVolume
          ? state.navigation.targetVolume->associatedLayer(
                targetIntersection.position)
          : nullptr;
    }
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
        // resolve the surfaces only - this could be changed to a resolve()
        // method
        state.navigation.navSurfaces
            = state.navigation.startLayer->getCompatibleSurfaces(
                navPar,
                forward,
                true,
                collectSensitive,
                collectMaterial,
                collectPassive,
                navigationLevel,
                state.startSurface,
                state.targetSurface);
        // the number of layer candidates
        if (state.navigation.navSurfaces.size()) {
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << state.navigation.navSurfaces.size();
            dstream << " surface candidates found.";
            return dstream.str();
          });
          // set the iterator
          state.navigation.navSurfaceIter
              = state.navigation.navSurfaces.begin();
          // update the navigation step size before you return
          updateStep(state, state.navigation.navSurfaceIter);
          return true;
        }
        return false;
      }
      // initialize layer - if it works go ahead with processing
      if (resolveLayers(navPar, state)) {
        if (state.stepping.stepSize == 0.) {
          debugLog(state, [&] {
            return std::string("On current layer surface, setting it.");
          });
          state.currentSurface = state.navigation.navLayerIter->representation;
          if (state.currentSurface)
            debugLog(state, [&] {
              std::stringstream dstream;
              dstream << "Current surface set to approach surface";
              dstream << state.currentSurface->geoID().toString();
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param[in] navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  /// @parqam skipCurrent is a flag in order to ignore the current iterator
  ///
  /// return (bool) triggers return to the stepper
  template <typename propagator_state_t>
  bool
  handleBoundaries(const NavigationParameters& navPar,
                   propagator_state_t&         state,
                   bool                        skipCurrent = false) const
  {
    // only handle boundaries if you are not in the target volume
    if (state.navigation.currentVolume == state.navigation.targetVolume
        || !state.navigation.currentVolume)
      return false;
    // if you came until here, and you have no boundaries
    // then retrieve them
    if (!state.navigation.navBoundaries.size()) {
      // get the navigation boundaries
      state.navigation.navBoundaries
          = state.navigation.currentVolume->boundarySurfacesOrdered(
              navPar, forward, skipCurrent);
      // the number of boundary candidates
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navBoundaries.size();
        dstream << " boundary surface candidates found.";
        return dstream.str();
      });
      // set the iterator - if we have boundary surfaces
      if (state.navigation.navBoundaries.size()) {
        state.navigation.navBoundaryIter
            = state.navigation.navBoundaries.begin();
        // update the navigation step size before you return
        updateStep(state, state.navigation.navBoundaryIter);
        return true;
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
      if (boundarySurface->isOnSurface(navPar.position(), true)) {
        debugLog(state, [&] {
          return std::string(
              "Boundary surface reached, prepare volume switch.");
        });
        // get the actual boundary for the navigation & the next volume
        auto boundary = state.navigation.navBoundaryIter->object;
        state.navigation.currentVolume = boundary->attachedVolume(
            navPar.position(), navPar.momentum(), forward);
        // no volume anymore : end of known world
        if (!state.navigation.currentVolume) {
          debugLog(state, [&] {
            return std::string(
                "No more volume to progress to, stopping navigation.");
          });
          return false;
        }
        // store the boundary for eventual actors to work on it
        state.currentSurface = boundarySurface;
        if (state.currentSurface)
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to boundary surface";
            dstream << state.currentSurface->geoID().toString();
            return dstream.str();
          });
        // and we can invalidate the boundary surfaces and return
        state.navigation.navBoundaries.clear();
        state.navigation.navBoundaryIter = state.navigation.navBoundaries.end();
        // resolve the new layer situation
        if (resolveLayers(navPar, state)) return true;
        // return
        debugLog(state, [&] {
          return std::string("No layers can be reached in the new volume.");
        });
        // self call for new boundaries
        return handleBoundaries(navPar, state, true);
      }
      ++state.navigation.navBoundaryIter;
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param[in] navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  ///
  /// @return indicates to return back to stepper
  template <typename propagator_state_t>
  bool
  resolveLayers(const NavigationParameters& navPar,
                propagator_state_t&         state) const
  {
    debugLog(state, [&] {
      return std::string("We do not have any layers yet, searching.");
    });
    // check if we are in the start volume
    bool start
        = (state.navigation.currentVolume == state.navigation.startVolume);
    // we do not have layers yet, get the candidates
    state.navigation.navLayers = state.navigation.currentVolume->decompose(
        (start ? state.navigation.startLayer : nullptr),
        state.navigation.targetLayer,
        navPar,
        true,
        state.stepping.stepSize.value(cstep::aborter),
        collectSensitive,
        collectMaterial,
        collectPassive);
    // the number of layer candidates
    if (state.navigation.navLayers.size()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navLayers.size();
        dstream << " layer candidates found.";
        return dstream.str();
      });
      // set the iterator
      state.navigation.navLayerIter = state.navigation.navLayers.begin();
      if (state.navigation.navLayerIter->object
          != state.navigation.startLayer) {
        debugLog(state,
                 [&] { return std::string("Stepping towards first layer."); });
        // update the navigation step size before you return
        updateStep(state, state.navigation.navLayerIter);
        return true;
      } else {
        debugLog(state, [&] {
          return std::string(
              "Start layer, avoid step to layer approach surface.");
        });
        return false;
      }
    }
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
    state.stepping.stepSize.update(
        state.stepping.stepSize.value(cstep::aborter), cstep::actor);
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  ///
  /// return (bool) triggers return to the stepper
  template <typename propagator_state_t>
  bool
  handleLayers(const NavigationParameters& navPar,
               propagator_state_t&         state) const
  {
    // if there are no layers, go back to the navigagor (not stepper yet)
    if (state.navigation.navLayers.empty()) return false;

    // loop over the available navigation layer candiates
    while (state.navigation.navLayerIter != state.navigation.navLayers.end()) {
      // we take the layer representation surface
      auto layer        = state.navigation.navLayerIter->object;
      auto layerSurface = state.navigation.navLayerIter->representation;
      auto layerVolume  = layer->representingVolume();
      // check if we are on the layer
      bool onLayer
          = state.navigation.navLayerIter->intersection.pathLength == 0;
      onLayer = onLayer || layerSurface->isOnSurface(navPar.position(), true);
      if (!onLayer
          && state.navigation.startLayer
              == state.navigation.navLayerIter->object) {
        onLayer = (layerVolume && layerVolume->inside(navPar.position()));
      }
      // check if we are on the layer
      if (onLayer) {
        // store the current surface in the propagator state
        if (state.startSurface
            && state.navigation.currentVolume == state.navigation.startVolume
            && layer == state.navigation.startLayer) {
          debugLog(state, [&] {
            return std::string("Switch layer surface to start surface.");
          });
          // setting layer surface & representation
          state.navigation.navLayerIter->representation = state.startSurface;
        } else {
          state.currentSurface = layerSurface;
          if (state.currentSurface)
            debugLog(state, [&] {
              std::stringstream dstream;
              dstream << "Current surface set to layer surface";
              dstream << state.currentSurface->geoID().toString();
              return dstream.str();
            });
        }
        // if you found surfaces return to the stepper
        if (resolveSurfaces(navPar, state)) return true;
        // increase the iterator
        ++state.navigation.navLayerIter;
      }
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
        // intersection with next layer
        auto layerIntersect = layerSurface->intersectionEstimate(
            navPar.position(), navPar.momentum(), true, false);
        // check if the intersect is invalid
        if (!layerIntersect) {
          debugLog(state, [&] {
            return std::string("Layer intersection not valid, skipping it.");
          });
          ++state.navigation.navLayerIter;
        } else {
          // update the navigation step size
          state.stepping.stepSize.update(
              state.stepping.navDir * layerIntersect.pathLength, cstep::actor);
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Navigation stepSize towards layer updated to ";
            dstream << state.stepping.stepSize.toString();
            return dstream.str();
          });
          return true;
        }
      }
    }
    // we are at the end of trying layers
    if (state.navigation.currentVolume == state.navigation.targetVolume) {
      debugLog(state, [&] {
        return std::string(
            "Done in final volume, release stepSize & proceed to target.");
      });
      // the step size will be set to the aborter step size
      state.stepping.stepSize.update(
          state.stepping.stepSize.value(cstep::aborter), cstep::actor);
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param navPar are the current navigation parameters
  /// @param state is the propagation state object
  ///
  /// @return whether you found surfaces or not
  template <typename propagator_state_t>
  bool
  resolveSurfaces(const NavigationParameters& navPar,
                  propagator_state_t&         state) const
  {
    // get the layer and layer surface
    auto layerSurface = state.navigation.navLayerIter->representation;
    auto navLayer     = state.navigation.navLayerIter->object;
    // are we on the start layer
    bool onStart      = (navLayer == state.navigation.startLayer);
    auto startSurface = onStart ? state.startSurface : layerSurface;
    // get the surfaces
    // @todo: could symmetrise with decompose() method
    state.navigation.navSurfaces
        = navLayer->getCompatibleSurfaces(navPar,
                                          forward,
                                          true,
                                          collectSensitive,
                                          collectMaterial,
                                          collectPassive,
                                          navigationLevel,
                                          startSurface,
                                          state.targetSurface);
    // the number of layer candidates
    if (state.navigation.navSurfaces.size()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navSurfaces.size();
        dstream << " surface candidates found.";
        return dstream.str();
      });
      // set the iterator
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      // update the navigation step size before you return to the stepper
      updateStep(state, state.navigation.navSurfaceIter);
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
  /// @tparam propagator_state_t is the state type of the propagagor
  ///
  /// @param[in] navPar are the current navigation parameters
  /// @param[in,out] state is the propagation state object
  ///
  /// return (bool) triggers a return to the stepper
  template <typename propagator_state_t>
  bool
  handleSurfaces(const NavigationParameters& navPar,
                 propagator_state_t&         state) const
  {
    // no surfaces, do not return to stepper
    if (!state.navigation.navSurfaces.size()) return false;
    // loop over the navigation surfaces
    while (state.navigation.navSurfaceIter
           != state.navigation.navSurfaces.end()) {
      // take the surface
      auto surface = state.navigation.navSurfaceIter->object;
      // case (s-a): we are at a surface
      //           only ask if you hadn't just come from a valid surface
      //
      // @todo : we should have a good idea whether this check should be done
      // @todo : add tolerance
      //
      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      if (surface->isOnSurface(navPar.position(), true)) {
        debugLog(state, [&] {
          return std::string("Surface successfully hit, storing it.");
        });
        // the surface will only appear due to correct
        // collect(Property) flag
        state.currentSurface = surface;
        if (state.currentSurface)
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to resolved surface";
            dstream << state.currentSurface->geoID().toString();
            return dstream.str();
          });
        // break if the surface is the target surface
        if (surface == state.targetSurface) {
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
      // case (s-b) : update step estimation to the new surface
      if (state.navigation.navSurfaceIter
          != state.navigation.navSurfaces.end()) {
        surface               = state.navigation.navSurfaceIter->object;
        auto surfaceIntersect = surface->intersectionEstimate(
            navPar.position(), navPar.momentum(), true, false);
        double surfaceDistance = surfaceIntersect.pathLength;
        if (!surfaceIntersect) {
          debugLog(state, [&] {
            return std::string(
                "Surface intersection is not valid, skipping it.");
          });
          ++state.navigation.navSurfaceIter;
          continue;
        } else {
          state.stepping.stepSize.update(
              state.stepping.navDir * surfaceDistance, cstep::actor);
          debugLog(state, [&] {
            std::stringstream dstream;
            dstream << "Navigation stepSize towards surface updated to ";
            dstream << state.stepping.stepSize.toString();
            return dstream.str();
          });
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
  /// @tparam propagator_state_t is the state type
  /// @tparam type_t is the intersection type (surface, layer, etc.)
  template <typename propagator_state_t, typename type_t>
  void
  updateStep(propagator_state_t& state, type_t& type) const
  {
    //  update the step
    double ustep = initialStepFactor * state.stepping.navDir
        * type->intersection.pathLength;
    state.stepping.stepSize.update(ustep, cstep::actor);
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
    if (debug) {
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
}

#endif  // end of namespace Acts
