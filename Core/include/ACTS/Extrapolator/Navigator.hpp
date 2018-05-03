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

// typedef std::vector<FullSurfaceIntersection<NavigationParameters>
// NavigationSurfaces;
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
/// The current navigation stage is cached in the result / and updated when
/// ever necessary. If any surface in the extrapolation / flow is hit, it is
/// set to the cache, such that other actors can / deal wit it.  / This actor
/// always needs to run first!  / It does two things: it figures out the order
/// of volumes, layers and / surfaces.  / For each propagation step, the
/// operator() runs, which checks if the / current / surface (or layer or
/// volume boundary) is reached via isOnSurface.  / The current target surface
/// is the surface pointed to by of the iterators / for the surfaces, layers
/// or
/// volume boundaries.  / If a surface is found, the pCache.currentSurface
/// pointer is set. This / enables subsequent actors to react. Secondly, this
/// actor uses the ordered / iterators / to figure out which surface, layer or
/// volume boundary is _supposed_ to be / hit / next. It then sets the maximum
/// step size to the path length found out by / straight line intersection. If
/// the isOnSurface call fails, it also / re-computes / the step size, to make
/// sure we end up at the desired surface.
struct Navigator
{
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

  /// Simple result struct to be returned
  /// It mainly acts as an interal state cache which is
  /// created for every propagation/extrapolation step
  struct this_state
  {
    /// Navigation on surface level
    /// the vector of navigation surfaces to work throgh
    NavigationSurfaces navSurfaces = {};
    /// the surface iterator
    NavigationSurfaceIter navSurfaceIter = navSurfaces.end();

    /// Navigation on layer level
    /// the vector of navigation layers to work throgh
    NavigationLayers    navLayers    = {};
    NavigationLayerIter navLayerIter = navLayers.end();

    /// Navigation on volume level
    // the vector of boundary surfaces to work throgh
    NavigationBoundaries   navBoundaries   = {};
    NavigationBoundaryIter navBoundaryIter = navBoundaries.end();

    /// Navigation cache: the start volume and layer
    const TrackingVolume* startVolume = nullptr;
    const Layer*          startLayer  = nullptr;
    /// Navigation cache: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation cache: the target volume and target layer
    const Layer*          targetLayer  = nullptr;
    const TrackingVolume* targetVolume = nullptr;

    /// break the navigation
    bool navigationBreak = false;
  };

  typedef this_state result_type;

  /// Navigation action for the ActionList of the Propagator
  ///
  /// @tparam propagator_cache_t is the type of Propagatgor cache
  /// @tparam stepper_cache_t is the type of Stepper cache
  ///
  /// @param pCache is the mutable stepper cache object
  /// @param sCache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& pCache,
             stepper_cache_t&    sCache,
             result_type&        result) const
  {
    // fail if you have no tracking geometry
    assert(trackingGeometry != nullptr);
    debugLog(
        pCache, result, [&] { return std::string("Entering navigator."); });

    // navigation parameters
    NavigationParameters navPar(sCache.position(),
                                sCache.navDir * sCache.direction());

    // Navigator always resets teh current surface first
    pCache.currentSurface = nullptr;

    // --------------------------------------------------------------------
    // Navigation break handling
    // This checks if a navigation break had been triggered:
    // - If so & the target exists or was hit - it simply returns
    // - If a target exists and was not yet hit, it checks for it
    // -> return is always to the stepper
    if (navigationBreak(navPar, pCache, result)) return;

    // -------------------------------------------------
    // Initialization
    // This should lead to:
    // - a current volume
    // - potentially also a current layer
    // -> return is always to the stepper
    if (initialize(navPar, pCache, sCache, result)) {
      debugLog(pCache, result, [&] {
        return std::string("Return to stepper - from initialize.");
      });
      return;
    }

    // -------------------------------------------------
    // Surfaces (if present)
    // - this can only happen after a layer has  sucessfully been resolved
    // -> return is always to the stepper
    if (handeSurfaces(navPar, pCache, sCache, result)) {
      debugLog(pCache, result, [&] {
        return std::string("Return to stepper - from surface handling.");
      });
      return;
    }

    // -------------------------------------------------
    // Layers are present
    // - this can only happen after a volume has successfully been resolved
    // -> return is always to the stepper
    if (handleLayers(navPar, pCache, sCache, result)) {
      debugLog(pCache, result, [&] {
        return std::string("Return to stepper - from layer handling.");
      });
      return;
    }
    // -------------------------------------------------
    // Volume be handled
    // - if you arrived
    // -> return is always to the stepper
    if (handleBoundaries(navPar, pCache, sCache, result)) {
      debugLog(pCache, result, [&] {
        return std::string("Return to stepper - from boundary handling.");
      });
      return;
    }
    // -------------------------------------------------
    // neither surfaces, layers nor boundaries triggered a return
    // navigation broken - switch navigator off
    result.navigationBreak = true;
    debugLog(pCache, result, [&] {
      return std::string("Naivgation break - no valid actions left.");
    });
    // release the navigation step size
    sCache.stepSize.release(cstep::actor);
    return;
  }

  /// Pure observer interface
  /// - this does not apply to the navigator
  ///
  /// @tparam propagator_cache_t is the type of Propagatgor cache
  /// @tparam stepper_cache_t is the type of Stepper cache
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& /*pCache*/, stepper_cache_t& /*sCache*/) const
  {
  }

  /// --------------------------------------------------------------------
  /// Navigation break handling
  ///
  /// This checks if a navigation break had been triggered
  /// - If so & the target exists or was hit - it simply returns
  /// - If a target exists and was not yet hit, it checks for it
  ///
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam result_t is the cache type
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param result is the result object
  ///
  /// boolean return triggers exit to stepper
  template <typename propagation_cache_t, typename result_type>
  bool
  navigationBreak(const NavigationParameters& navPar,
                  propagation_cache_t&        pCache,
                  result_type&                result) const
  {
    if (result.navigationBreak) {
      // target exists and reached, or no target exists
      if (pCache.targetReached || !pCache.targetSurface) return true;
      // the only advande could have been to the target
      if (pCache.targetSurface->isOnSurface(navPar.position(), true)) {
        // set the target surface
        pCache.currentSurface = pCache.targetSurface;
        debugLog(pCache, result, [&] {
          std::stringstream dstream;
          dstream << "Current surface set to target surface ";
          dstream << pCache.currentSurface->geoID().toString();
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
  /// size toward the first layer.  The return prevent execution of the
  /// subsequent logic, we want to make a step first.
  ///
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  ///
  /// @return boolean trigger if successful
  template <typename propagtor_cache_t,
            typename stepper_cache_t,
            typename result_type>
  bool
  initialize(const NavigationParameters& navPar,
             propagtor_cache_t&          pCache,
             stepper_cache_t&            sCache,
             result_type&                result) const
  {

    // no initialisation necessary
    if (result.currentVolume) return false;

    debugLog(pCache, result, [&] {
      return std::string("Initializing start volume.");
    });

    // we set the current surface to the start surface
    // for eventual post-update actio, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    pCache.currentSurface = pCache.startSurface;
    if (pCache.currentSurface)
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << "Current surface set to start surface ";
        dstream << pCache.currentSurface->geoID().toString();
        return dstream.str();
      });

    // Fast Navigation initialization for start condition:
    // short-cut through object association, saves navigation in the
    // geometry and volume tree search for the lowest volume
    if (pCache.startSurface && pCache.startSurface->associatedLayer()) {
      debugLog(pCache, result, [&] {
        return std::string("Fast start initialization through association.");
      });
      // assign the current layer and volume by association
      result.startLayer  = pCache.startSurface->associatedLayer();
      result.startVolume = result.startLayer->trackingVolume();
    } else {
      debugLog(pCache, result, [&] {
        return std::string("Slow start initialization through search.");
      });

      // current volume and layer search through global search
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << "Starting from position (" << navPar.position().x();
        dstream << ", " << navPar.position().y();
        dstream << ", " << navPar.position().z() << ")";
        return dstream.str();
      });
      result.startVolume
          = trackingGeometry->lowestTrackingVolume(navPar.position());
      result.startLayer = result.startVolume
          ? result.startVolume->associatedLayer(navPar.position())
          : nullptr;
    }
    // Fast Navigation initialization for target:
    if (pCache.targetSurface && pCache.targetSurface->associatedLayer()) {
      debugLog(pCache, result, [&] {
        return std::string("Fast target initialization through association.");
      });
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << "Target surface set to";
        dstream << pCache.targetSurface->geoID().toString();
        return dstream.str();
      });
      // assign the target volume and the target surface
      result.targetLayer  = pCache.targetSurface->associatedLayer();
      result.targetVolume = result.targetLayer->trackingVolume();
    } else if (pCache.targetSurface) {
      // Slow navigation initialization for target:
      // target volume and layer search through global search
      auto targetIntersection = pCache.targetSurface->intersectionEstimate(
          navPar.position(), navPar.momentum(), true, false);
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << "Target estimate position (";
        dstream << targetIntersection.position.x() << ", ";
        dstream << targetIntersection.position.y() << ", ";
        dstream << targetIntersection.position.z() << ")";
        return dstream.str();
      });
      /// get the target volume from the intersection
      result.targetVolume
          = trackingGeometry->lowestTrackingVolume(targetIntersection.position);
      result.targetLayer = result.targetVolume
          ? result.targetVolume->associatedLayer(targetIntersection.position)
          : nullptr;
    }
    // A current volume exists
    if (result.startVolume) {
      // assign to the currentVolume
      result.currentVolume = result.startVolume;
      // fast exit if start and target layer are identical
      if (result.startLayer == result.targetLayer) {
        debugLog(pCache, result, [&] {
          return std::string(
              "Start and target layer identical, check surfaces.");
        });
        // resolve the surfaces only - this could be changed to a resolve()
        // method
        result.navSurfaces
            = result.startLayer->getCompatibleSurfaces(navPar,
                                                       forward,
                                                       true,
                                                       collectSensitive,
                                                       collectMaterial,
                                                       collectPassive,
                                                       navigationLevel,
                                                       pCache.startSurface,
                                                       pCache.targetSurface);
        // the number of layer candidates
        if (result.navSurfaces.size()) {
          debugLog(pCache, result, [&] {
            std::stringstream dstream;
            dstream << result.navSurfaces.size();
            dstream << " surface candidates found.";
            return dstream.str();
          });
          // set the iterator
          result.navSurfaceIter = result.navSurfaces.begin();
          // update the navigation step size before you return
          updateStep(pCache, sCache, result, result.navSurfaceIter);
          return true;
        }
        return false;
      }
      // initialize layer - if it works go ahead
      if (resolveLayers(navPar, pCache, sCache, result)) {
        if (sCache.stepSize == 0.) {
          debugLog(pCache, result, [&] {
            return std::string("On current layer surface, setting it.");
          });
          pCache.currentSurface = result.navLayerIter->representation;
          if (pCache.currentSurface)
            debugLog(pCache, result, [&] {
              std::stringstream dstream;
              dstream << "Current surface set to approach surface";
              dstream << pCache.currentSurface->geoID().toString();
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
    result.navigationBreak = true;
    debugLog(pCache, result, [&] {
      return std::string("Navigation broken, pure propagation.");
    });
    // release the navigation step size
    sCache.stepSize.release(cstep::actor);
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
  /// line intersect if found, the boundary surface is skipped.
  /// If we are out of boundary surfaces, the navigation is terminated.
  ///
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  ///
  /// return (bool) triggers return to the stepper
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t>
  bool
  handleBoundaries(const NavigationParameters& navPar,
                   propagator_cache_t&         pCache,
                   stepper_cache_t&            sCache,
                   result_t&                   result,
                   bool                        skipCurrent = false) const
  {
    // only handle boundaries if you are not in the target volume
    if (result.currentVolume == result.targetVolume || !result.currentVolume)
      return false;
    // if you came until here, and you have no boundaries
    // then retrieve them
    if (!result.navBoundaries.size()) {
      // get the navigation boundaries
      result.navBoundaries = result.currentVolume->boundarySurfacesOrdered(
          navPar, forward, skipCurrent);
      // the number of boundary candidates
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << result.navBoundaries.size();
        dstream << " boundary surface candidates found.";
        return dstream.str();
      });
      // set the iterator - if we have boundary surfaces
      if (result.navBoundaries.size()) {
        result.navBoundaryIter = result.navBoundaries.begin();
        // update the navigation step size before you return
        updateStep(pCache, sCache, result, result.navBoundaryIter);
        return true;
      } else {
        debugLog(pCache, result, [&] {
          return std::string(
              "No valid boundary surface found, stopping navigation.");
        });
        return false;
      }
    }
    // loop over rest of the boundaries
    while (result.navBoundaryIter != result.navBoundaries.end()) {
      auto boundarySurface = result.navBoundaryIter->representation;
      // case (v-a) : you are on the boundary surface
      //              only if you hadn't just done a volume switch
      // check if we are on already in this step
      if (boundarySurface->isOnSurface(navPar.position(), true)) {
        debugLog(pCache, result, [&] {
          return std::string(
              "Boundary surface reached, prepare volume switch.");
        });
        // get the actual boundary for the navigation & the next volume
        auto boundary        = result.navBoundaryIter->object;
        result.currentVolume = boundary->attachedVolume(
            navPar.position(), navPar.momentum(), forward);
        // no volume anymore : end of known world
        if (!result.currentVolume) {
          debugLog(pCache, result, [&] {
            return std::string(
                "No more volume to progress to, stopping navigation.");
          });
          return false;
        }
        // store the boundary for eventual actors to work on it
        pCache.currentSurface = boundarySurface;
        if (pCache.currentSurface)
          debugLog(pCache, result, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to boundary surface";
            dstream << pCache.currentSurface->geoID().toString();
            return dstream.str();
          });
        // and we can invalidate the boundary surfaces and return
        result.navBoundaries.clear();
        result.navBoundaryIter = result.navBoundaries.end();
        // resolve the new layer situation
        if (resolveLayers(navPar, pCache, sCache, result)) return true;
        // return
        debugLog(pCache, result, [&] {
          return std::string("No layers can be reached in the new volume.");
        });
        // self call for new boundaries
        return handleBoundaries(navPar, pCache, sCache, result, true);
      }
      ++result.navBoundaryIter;
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
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  /// @param result is the result object
  ///
  /// @return indicates to return back to stepper
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t>
  bool
  resolveLayers(const NavigationParameters& navPar,
                propagator_cache_t&         pCache,
                stepper_cache_t&            sCache,
                result_t&                   result) const
  {
    debugLog(pCache, result, [&] {
      return std::string("We do not have any layers yet, searching.");
    });
    // check if we are in the start volume
    bool start = (result.currentVolume == result.startVolume);
    // we do not have layers yet, get the candidates
    result.navLayers
        = result.currentVolume->decompose((start ? result.startLayer : nullptr),
                                          result.targetLayer,
                                          navPar,
                                          true,
                                          sCache.stepSize.value(cstep::aborter),
                                          collectSensitive,
                                          collectMaterial,
                                          collectPassive);
    // the number of layer candidates
    if (result.navLayers.size()) {
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << result.navLayers.size();
        dstream << " layer candidates found.";
        return dstream.str();
      });
      // set the iterator
      result.navLayerIter = result.navLayers.begin();
      if (result.navLayerIter->object != result.startLayer) {
        debugLog(pCache, result, [&] {
          return std::string("Stepping towards first layer.");
        });
        // update the navigation step size before you return
        updateStep(pCache, sCache, result, result.navLayerIter);
        return true;
      } else {
        debugLog(pCache, result, [&] {
          return std::string(
              "Start layer, avoid step to layer approach surface.");
        });
        return false;
      }
    }
    if (result.currentVolume != result.targetVolume) {
      debugLog(pCache, result, [&] {
        return std::string("No layer candidates found, switching volume.");
      });
      return false;
    }
    debugLog(pCache, result, [&] {
      return std::string(
          "Done in final volume, release stepSize & proceed to target.");
    });
    // the step size will be set to the aborter step size
    sCache.stepSize.update(sCache.stepSize.value(cstep::aborter), cstep::actor);
    result.navigationBreak = true;
    return true;
  }

  // Loop over layer candidates.
  //
  // We are now trying to advance to the next layer (with surfaces)
  // Check if we are on the representing surface of the layer pointed
  // at by navLayerIter. If so, we unpack the compatible surfaces
  // (determined by straight line intersect), and set up the iterator
  // so that the next call to operator() will enter the surface
  // check mode above. If no surfaces are found, we skip the layer.
  // If we unpack a surface, the step size is set to the path length
  // to the first surface, as determined by straight line intersect.
  //
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  /// @param result is the result object
  //
  // return (bool) triggers return to the stepper
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t>
  bool
  handleLayers(const NavigationParameters& navPar,
               propagator_cache_t&         pCache,
               stepper_cache_t&            sCache,
               result_t&                   result) const
  {
    // of course only
    if (result.navLayers.size()) {
      while (result.navLayerIter != result.navLayers.end()) {
        // we take the layer representation surface
        auto layer        = result.navLayerIter->object;
        auto layerSurface = result.navLayerIter->representation;
        auto layerVolume  = layer->representingVolume();
        // check if we are on the layer
        bool onLayer = result.navLayerIter->intersection.pathLength == 0;
        onLayer = onLayer || layerSurface->isOnSurface(navPar.position(), true);
        if (!onLayer && result.startLayer == result.navLayerIter->object) {
          onLayer = (layerVolume && layerVolume->inside(navPar.position()));
        }
        // check if we are on the layer
        if (onLayer) {
          // store the current surface in the cache
          if (pCache.startSurface && result.currentVolume == result.startVolume
              && layer == result.startLayer) {
            debugLog(pCache, result, [&] {
              return std::string("Switch layer surface to start surface.");
            });
            // setting layer surface & representation
            result.navLayerIter->representation = pCache.startSurface;
          } else {
            pCache.currentSurface = layerSurface;
            if (pCache.currentSurface)
              debugLog(pCache, result, [&] {
                std::stringstream dstream;
                dstream << "Current surface set to layer surface";
                dstream << pCache.currentSurface->geoID().toString();
                return dstream.str();
              });
          }
          // if you found surfaces return to the stepper
          if (resolveSurfaces(navPar, pCache, sCache, result)) return true;
          // increase the iterator
          ++result.navLayerIter;
        }
        if (result.navLayerIter != result.navLayers.end()) {
          // update in case a switch was done
          layerSurface = result.navLayerIter->representation;
          // we are not on the layer
          debugLog(pCache, result, [&] {
            std::stringstream dstream;
            dstream << std::distance(result.navLayerIter,
                                     result.navLayers.end());
            dstream << " out of " << result.navLayers.size();
            dstream << " layers remain to try.";
            return dstream.str();
          });
          // intersection with next layer
          auto layerIntersect = layerSurface->intersectionEstimate(
              navPar.position(), navPar.momentum(), true, false);
          // check if the intersect is invalid
          if (!layerIntersect) {
            debugLog(pCache, result, [&] {
              return std::string("Layer intersection not valid, skipping it.");
            });
            ++result.navLayerIter;
          } else {
            // update the navigation step size
            sCache.stepSize.update(sCache.navDir * layerIntersect.pathLength,
                                   cstep::actor);
            debugLog(pCache, result, [&] {
              std::stringstream dstream;
              dstream << "Navigation stepSize towards layer updated to ";
              dstream << sCache.stepSize.toString();
              return dstream.str();
            });
            return true;
          }
        }
      }
      // we are at the end of trying layers
      if (result.currentVolume == result.targetVolume) {
        debugLog(pCache, result, [&] {
          return std::string(
              "Done in final volume, release stepSize & proceed to target.");
        });
        // the step size will be set to the aborter step size
        sCache.stepSize.update(sCache.stepSize.value(cstep::aborter),
                               cstep::actor);
        result.navigationBreak = true;
        return true;
      }
      debugLog(pCache, result, [&] {
        return std::string("All layers been handled, switching volume.");
      });
      // clear the layers
      result.navLayers.clear();
      result.navLayerIter = result.navLayers.end();
      // clear the boundaries
      result.navBoundaries.clear();
      result.navBoundaryIter = result.navBoundaries.end();
    }
    // do not return to the stepper here
    return false;
  }

  /// Resolve the surfaces of this layer, if not the start layer
  ///
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  /// @param result is the result object
  //
  /// @return whether you found surfaces or not
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t>
  bool
  resolveSurfaces(const NavigationParameters& navPar,
                  propagator_cache_t&         pCache,
                  stepper_cache_t&            sCache,
                  result_t&                   result) const
  {
    // get the layer and layer surface
    auto layerSurface = result.navLayerIter->representation;
    auto navLayer     = result.navLayerIter->object;
    // are we on the start layer
    bool onStart      = (navLayer == result.startLayer);
    auto startSurface = onStart ? pCache.startSurface : layerSurface;
    // get the surfaces
    // @todo: could symmetrise with decompose() method
    result.navSurfaces = navLayer->getCompatibleSurfaces(navPar,
                                                         forward,
                                                         true,
                                                         collectSensitive,
                                                         collectMaterial,
                                                         collectPassive,
                                                         navigationLevel,
                                                         startSurface,
                                                         pCache.targetSurface);
    // the number of layer candidates
    if (result.navSurfaces.size()) {
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << result.navSurfaces.size();
        dstream << " surface candidates found.";
        return dstream.str();
      });
      // set the iterator
      result.navSurfaceIter = result.navSurfaces.begin();
      // update the navigation step size before you return to the stepper
      updateStep(pCache, sCache, result, result.navSurfaceIter);
      return true;
    }
    result.navSurfaceIter = result.navSurfaces.end();
    debugLog(pCache, result, [&] {
      return std::string("No surface candidates found, switching layer.");
    });
    return false;
  }

  /// Loop over surface candidates here:
  ///  - if an intersect is  valid but not yet reached
  ///    then return with updated step size
  ///  - if an intersect is not valid, switch to next
  ///
  /// @tparam propagator_cache_t is the cache type of the propagagor
  /// @tparam stepper_cache_t is the cache type of the stepper
  /// @tparam result_t is the cache type
  ///
  /// @param navPar are the current navigation parameters
  /// @param pCache is the propagation cache object
  /// @param sCache is the stepper cache object
  ///
  /// return (bool) triggers a return to the stepper
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t>
  bool
  handeSurfaces(const NavigationParameters& navPar,
                propagator_cache_t&         pCache,
                stepper_cache_t&            sCache,
                result_t&                   result) const
  {
    // no surfaces, do not return to stepper
    if (!result.navSurfaces.size()) return false;
    // loop over the navigation surfaces
    while (result.navSurfaceIter != result.navSurfaces.end()) {
      // take the surface
      auto surface = result.navSurfaceIter->object;
      // case (s-a): we are at a surface
      //           only ask if you hadn't just come from a valid surface
      //
      // @todo : we should have a good idea whether this check should be done
      // @todo : add tolerance
      //
      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      if (surface->isOnSurface(navPar.position(), true)) {
        debugLog(pCache, result, [&] {
          return std::string("Surface successfully hit, storing it.");
        });
        // the surface will only appear due to correct
        // collect(Property) flag
        pCache.currentSurface = surface;
        if (pCache.currentSurface)
          debugLog(pCache, result, [&] {
            std::stringstream dstream;
            dstream << "Current surface set to resolved surface";
            dstream << pCache.currentSurface->geoID().toString();
            return dstream.str();
          });
        // break if the surface is the target surface
        if (surface == pCache.targetSurface) {
          debugLog(pCache, result, [&] {
            return std::string("This was the target surface. Done.");
          });
          return true;
        }
        // switch to the next candidate
        ++result.navSurfaceIter;
      }
      // screen output how much is left to try
      debugLog(pCache, result, [&] {
        std::stringstream dstream;
        dstream << std::distance(result.navSurfaceIter,
                                 result.navSurfaces.end());
        dstream << " out of " << result.navSurfaces.size();
        dstream << " surfaces remain to try.";
        return dstream.str();
      });
      // case (s-b) : update step estimation to the new surface
      if (result.navSurfaceIter != result.navSurfaces.end()) {
        surface               = result.navSurfaceIter->object;
        auto surfaceIntersect = surface->intersectionEstimate(
            navPar.position(), navPar.momentum(), true, false);
        double surfaceDistance = surfaceIntersect.pathLength;
        if (!surfaceIntersect) {
          debugLog(pCache, result, [&] {
            return std::string(
                "Surface intersection is not valid, skipping it.");
          });
          ++result.navSurfaceIter;
          continue;
        } else {
          sCache.stepSize.update(sCache.navDir * surfaceDistance, cstep::actor);
          debugLog(pCache, result, [&] {
            std::stringstream dstream;
            dstream << "Navigation stepSize towards surface updated to ";
            dstream << sCache.stepSize.toString();
            return dstream.str();
          });
          return true;
        }
      }
    }
    // case (s-c) : reached the end of the surface iteration
    if (result.navSurfaceIter == result.navSurfaces.end()) {
      debugLog(pCache, result, [&] {
        return std::string("Last surface on layer reached, switching layer.");
      });
      // first clear the surface cache
      result.navSurfaces.clear();
      result.navSurfaceIter = result.navSurfaces.end();
      // now switch to the next layer
      ++result.navLayerIter;
    }
    // do not return to the stepper
    return false;
  }

  /// This method updates the constrained step size
  /// @tparam stepper_cache_t is the cache type
  /// @tparam result_t is the cache type
  /// @tparam type_t is the cache type
  template <typename propagator_cache_t,
            typename stepper_cache_t,
            typename result_t,
            typename type_t>
  void
  updateStep(propagator_cache_t& pCache,
             stepper_cache_t&    sCache,
             result_t&           result,
             type_t&             type) const
  {
    //  update the step
    double ustep
        = initialStepFactor * sCache.navDir * type->intersection.pathLength;
    sCache.stepSize.update(ustep, cstep::actor);
    debugLog(pCache, result, [&] {
      std::stringstream dstream;
      dstream << "Navigation stepSize updated to ";
      dstream << sCache.stepSize.toString();
      return dstream.str();
    });
  }

private:
  /// The private navigation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the pCache.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_cache_t Type of the propagator cache
  /// @result_t Type of the reulst
  ///
  /// @param pCache the propagator cache for the debug flag, prefix and length
  /// @param result the result object of the navigator
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_cache_t, typename result_t>
  void
  debugLog(propagator_cache_t&          pCache,
           result_t&                    result,
           std::function<std::string()> logAction) const
  {
    if (debug) {
      std::string vName               = "No Volume";
      if (result.currentVolume) vName = result.currentVolume->volumeName();
      std::stringstream dstream;
      dstream << ">>>" << std::setw(pCache.debugPfxWidth) << vName << " | ";
      dstream << std::setw(pCache.debugMsgWidth) << logAction() << '\n';
      pCache.debugString += dstream.str();
    }
  }
};
}

#endif  // end of namespace Acts
