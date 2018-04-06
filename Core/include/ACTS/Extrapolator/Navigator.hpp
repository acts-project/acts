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
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Propagator/detail/constrained_step.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"

#ifndef NAVIGATOR_DEBUG_OUTPUTS
#define NAVIGATOR_DEBUG_OUTPUTS
#define NAVLOG(cache, result, message)                                         \
  if (debug) {                                                                 \
    std::string vName                = "No Volume";                            \
    if (result.current_volume) vName = result.current_volume->volumeName();    \
    std::stringstream dstream;                                                 \
    dstream << ">>>" << std::setw(cache.debug_pfx_width) << vName << " | ";    \
    dstream << std::setw(cache.debug_msg_width) << message << '\n';            \
    cache.debug_string += dstream.str();                                       \
  }
#endif

namespace Acts {

typedef detail::constrained_step cstep;

/// Struct to mimmick track parameters
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
/// volume boundaries.  / If a surface is found, the cache.current_surface
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
  bool collectPassive = true;

  /// store the debug message
  bool debug = false;

  /// Simple result struct to be returned
  /// It mainly acts as an interal state cache which is
  /// created for every propagation/extrapolation step
  struct this_state
  {
    /// Navigation on surface level
    /// the vector of navigation surfaces to work throgh
    NavigationSurfaces nav_surfaces = {};
    /// the surface iterator
    NavigationSurfaceIter nav_surface_iter = nav_surfaces.end();

    /// Navigation on layer level
    /// the vector of navigation layers to work throgh
    NavigationLayers    nav_layers     = {};
    NavigationLayerIter nav_layer_iter = nav_layers.end();

    /// Navigation on volume level
    // the vector of boundary surfaces to work throgh
    NavigationBoundaries   nav_boundaries    = {};
    NavigationBoundaryIter nav_boundary_iter = nav_boundaries.end();

    /// Navigation cache: the start volume and layer
    const TrackingVolume* start_volume = nullptr;
    const Layer*          start_layer  = nullptr;
    /// Navigation cache: the current volume
    const TrackingVolume* current_volume = nullptr;
    /// Navigation cache: the target volume and target layer
    const Layer*          target_layer  = nullptr;
    const TrackingVolume* target_volume = nullptr;

    /// check if we stay on the layer
    bool stay_on_layer = false;

    /// break the navigation
    bool navigation_break = false;
  };

  typedef this_state result_type;

  /// Navigation action for the ActionList of the Propagator
  ///
  /// @tparam cache_t is the type of Stepper cache
  ///
  /// @param cache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename cache_t>
  void
  operator()(cache_t& cache, result_type& result) const
  {

    // reset navigation if stream was broken
    if (result.navigation_break) return;

    // clean-up the navigator
    if (cache.target_reached) {
      result.current_volume   = nullptr;
      result.navigation_break = false;
      result.stay_on_layer    = false;
      return;
    }

    // fail if you have no tracking geometry
    assert(trackingGeometry != nullptr);

    // the navigator erases the current surfaces
    cache.current_surface = nullptr;

    // navigation parameters
    NavigationParameters nav_par(cache.position(),
                                 cache.nav_dir * cache.direction());
    // -------------------------------------------------
    // Initialization
    // - this should lead to
    // - a current volume
    // - potentially also a current layer
    // -> return is always to the stepper
    if (initialize(nav_par, cache, result)) return;

    // -------------------------------------------------
    // Surfaces (if present)
    // - this can only happen after a layer has  sucessfully been resolved
    // -> return is always to the stepper
    if (handle_surfaces(nav_par, cache, result)) return;

    // -------------------------------------------------
    // Layers are present
    // - this can only happen after a volume has successfully been resolved
    // -> return is always to the stepper
    if (handle_layers(nav_par, cache, result)) return;

    // -------------------------------------------------
    // Volume be handled
    // - if you arrived
    // -> return is always to the stepper
    if (handle_boundaries(nav_par, cache, result)) return;

    // -------------------------------------------------
    // neither surfaces, layers nor boundaries triggered a return
    // navigation broken - switch navigator off
    result.navigation_break = true;
    NAVLOG(cache, result, "No more navigation actions, pure propagation.");
    return;
  }

  /// Pure observer interface
  /// This does not apply to the navigator
  template <typename cache_t>
  void
  operator()(cache_t& cache) const
  {
    (void)cache;
  }

  // --------------------------------------------------------------------
  // Navigation initialisation
  //
  // This is only called once for every propagation/extrapolation
  //
  // ---------------------------------------------------------------------
  // Check for navigation initialisation & do it if necessary.  This means
  // we do not have an active volume yet (since we just started).  We get
  // the innermost volume, and set up an ordered layer iterator and step
  // size toward the first layer.  The return prevent execution of the
  // subsequent logic, we want to make a step first.
  //
  //
  template <typename cache_t, typename result_type>
  bool
  initialize(const NavigationParameters& nav_par,
             cache_t&                    cache,
             result_type&                result) const
  {

    // no initialisation necessary
    if (result.current_volume) return false;

    NAVLOG(cache, result, "Initializing start volume.");
    // we set the current surface to the start surface
    // for eventual post-update actio, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    cache.current_surface = cache.start_surface;

    // Fast Navigation initialization for start condition:
    // short-cut through object association, saves navigation in the
    // geometry and volume tree search for the lowest volume
    if (cache.start_surface && cache.start_surface->associatedLayer()) {
      NAVLOG(cache, result, "Fast start initialization through association.")
      // assign the current layer and volume by association
      result.start_layer  = cache.start_surface->associatedLayer();
      result.start_volume = result.start_layer->trackingVolume();
    } else {
      NAVLOG(cache, result, "Slow start initialization through search.")
      // Slow navigation initialization for start condition:
      // current volume and layer search through global search
      result.start_volume
          = trackingGeometry->lowestTrackingVolume(nav_par.position());
      result.start_layer = result.current_volume
          ? result.current_volume->associatedLayer(nav_par.position())
          : nullptr;
    }
    // Fast Navigation initialization for target:
    if (cache.target_surface && cache.target_surface->associatedLayer()) {
      // assign the target volume and the target surface
      result.target_layer  = cache.target_surface->associatedLayer();
      result.target_volume = result.target_layer->trackingVolume();
    } else if (cache.target_surface) {
      // Slow navigation initialization for target:
      // target volume and layer search through global search
      auto target_intersection = cache.target_surface->intersectionEstimate(
          nav_par.position(), nav_par.momentum(), true, false);
      /// get the target volume from the intersection
      result.target_volume = trackingGeometry->lowestTrackingVolume(
          target_intersection.position);
      result.target_layer = result.target_volume
          ? result.target_volume->associatedLayer(target_intersection.position)
          : nullptr;
    }
    // A current volume exists
    if (result.start_volume) {
      // assign to the current_volume
      result.current_volume = result.start_volume;
      // initialize layer - if it works go ahead
      if (resolve_layers(nav_par, cache, result)) {
        if (cache.step_size == 0.) {
          NAVLOG(cache, result, "On the current layer surface, setting it.");
          cache.current_surface = result.nav_layer_iter->representation;
          // no returning to the stepper here
          return false;
        }
        return true;
      }
      return false;
    }
    // navigation broken - switch navigator off
    result.navigation_break = true;
    NAVLOG(cache, result, "Navigation broken, pure propagation.");
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
  /// layer.  If we don't find a next volume, the navigation_break
  /// indicator is set.  This ends the navigation. Finally, the boundary
  /// iterator is cleared, so that the subsequent call goes back to
  /// the layer iteration logic.
  ///
  /// If we are not on the current boundary surface, we try the next one.
  /// The iterator is advanced and the step size is set. If no straight
  /// line intersect if found, the boundary surface is skipped.
  /// If we are out of
  /// boundary surfaces, the navigation is terminated.
  /// return (bool) triggers return to the stepper
  template <typename cache_t, typename result_t>
  bool
  handle_boundaries(const NavigationParameters& nav_par,
                    cache_t&                    cache,
                    result_t&                   result) const
  {
    // if you came until here, and you have no boundaries
    // then retrieve them
    if (!result.nav_boundaries.size()) {
      // get the navigation boundaries
      result.nav_boundaries
          = result.current_volume->boundarySurfacesOrdered(nav_par, forward);
      // the number of boundary candidates
      NAVLOG(cache,
             result,
             result.nav_boundaries.size()
                 << " boundary surface candidates found.");
      // set the iterator
      result.nav_boundary_iter = result.nav_boundaries.begin();
      // update the navigation step size before you return
      update_step(cache, result, result.nav_boundary_iter);
      return true;
    }
    // loop over rest of the boundaries
    while (result.nav_boundary_iter != result.nav_boundaries.end()) {
      auto boundary_surface = result.nav_boundary_iter->representation;
      // case (v-a) : you are on the boundary surface
      //              only if you hadn't just done a volume switch
      // check if we are on already in this step
      if (boundary_surface->isOnSurface(nav_par.position(), true)) {
        NAVLOG(
            cache, result, "Boundary surface reached, prepare volume switch.");
        // get the actual boundary for the navigation & the next volume
        auto boundary         = result.nav_boundary_iter->object;
        result.current_volume = boundary->attachedVolume(
            nav_par.position(), nav_par.momentum(), forward);
        // store the boundary for eventual actors to work on it
        cache.current_surface = boundary_surface;
        // and we can invalidate the boundary surfaces and return
        result.nav_boundaries.clear();
        result.nav_boundary_iter = result.nav_boundaries.end();
        // resolve the new layer situation
        if (resolve_layers(nav_par, cache, result)) return true;
        // return
        NAVLOG(cache, result, "No layers can be reached in the new volume.");
        // self call for new boundaries
        return handle_boundaries(nav_par, cache, result);
      }
      ++result.nav_boundary_iter;
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
  /// @tparam cache_t is the cache type
  /// @tparam result_t is the cache type
  ///
  /// @return indicates to return back to stepper
  template <typename cache_t, typename result_t>
  bool
  resolve_layers(const NavigationParameters& nav_par,
                 cache_t&                    cache,
                 result_t&                   result) const
  {
    NAVLOG(cache, result, "We do not have any layers yet, searching.");
    // check if we are in the start volume
    bool start = (result.current_volume == result.start_volume);
    // we do not have layers yet, get the candidates
    result.nav_layers = result.current_volume->decompose(
        (start ? result.start_layer : nullptr),
        result.target_layer,
        nav_par,
        true,
        cache.step_size.value(cstep::aborter),
        collectSensitive,
        collectMaterial,
        collectPassive);
    // the number of layer candidates
    if (result.nav_layers.size()) {
      NAVLOG(cache,
             result,
             result.nav_layers.size() << " layer candidates found.");
      // set the iterator
      result.nav_layer_iter = result.nav_layers.begin();
      if (result.nav_layer_iter->object != result.start_layer) {
        NAVLOG(cache, result, "Stepping towards first layer.");
        // update the navigation step size before you return
        update_step(cache, result, result.nav_layer_iter);
        return true;
      } else {
        NAVLOG(cache,
               result,
               "Start layer, avoid step to layer approach surface.");
        return false;
      }
    }
    NAVLOG(cache, result, "No layer candidates found, switching volume.");
    return false;
  }

  // Loop over layer candidates.
  //
  // We are now trying to advance to the next layer (with surfaces)
  // Check if we are on the representing surface of the layer pointed
  // at by nav_layer_iter. If so, we unpack the compatible surfaces
  // (determined by straight line intersect), and set up the iterator
  // so that the next call to operator() will enter the surface
  // check mode above. If no surfaces are found, we skip the layer.
  // If we unpack a surface, the step size is set to the path length
  // to the first surface, as determined by straight line intersect.
  //
  // return (bool) triggers return to the stepper
  template <typename cache_t, typename result_t>
  bool
  handle_layers(const NavigationParameters& nav_par,
                cache_t&                    cache,
                result_t&                   result) const
  {
    // of course only
    if (result.nav_layers.size()) {
      while (result.nav_layer_iter != result.nav_layers.end()) {
        // we take the layer representation surface
        auto layer_surface = result.nav_layer_iter->representation;
        auto layer_volume = result.nav_layer_iter->object->representingVolume();
        // check if we are on the layer
        bool on_layer = result.nav_layer_iter->intersection.pathLength == 0;
        on_layer
            = on_layer || layer_surface->isOnSurface(nav_par.position(), true);
        if (!on_layer && result.start_layer == result.nav_layer_iter->object) {
          on_layer = (layer_volume && layer_volume->inside(nav_par.position()));
        }
        // check if we are on the layer
        if (on_layer) {
          NAVLOG(cache, result, "On layer, storing & resolve.");
          // store the current surface in the cache
          cache.current_surface = layer_surface;
          // if you found surfaces return to the stepper
          if (resolve_surfaces(nav_par, cache, result)) return true;
          // increase the iterator
          ++result.nav_layer_iter;
        }
        if (result.nav_layer_iter != result.nav_layers.end()) {
          // update in case a switch was done
          layer_surface = result.nav_layer_iter->representation;
          // we are not on the layer
          NAVLOG(cache,
                 result,
                 std::distance(result.nav_layer_iter, result.nav_layers.end())
                     << " out of "
                     << result.nav_layers.size()
                     << " layers remain to try.");
          // intersection with next layer
          auto layer_intersect = layer_surface->intersectionEstimate(
              nav_par.position(), nav_par.momentum(), true, false);
          // check if the intersect is invalid
          if (!layer_intersect) {
            NAVLOG(cache, result, "Layer intersection not valid, skipping it.");
            ++result.nav_layer_iter;
          } else {
            // update the navigation step size
            cache.step_size.update(cache.nav_dir * layer_intersect.pathLength,
                                   cstep::actor);
            NAVLOG(cache,
                   result,
                   "Navigation step_size towards layer updated to "
                       << cache.step_size.toString());
            return true;
          }
        }
      }
      // we are at the end of trying layers
      NAVLOG(cache, result, "All layers been handled, switching volume.");
      // clear the layers
      result.nav_layers.clear();
      result.nav_layer_iter = result.nav_layers.end();
      // clear the boundaries
      result.nav_boundaries.clear();
      result.nav_boundary_iter = result.nav_boundaries.end();
    }
    // do not return to the stepper here
    return false;
  }

  // Navigation through surfaces
  // -----------------------------------------------
  ///
  // If we are currently iterating over surfaces, we get the current one
  // pointed at by the iterator here.  This is the main routine that
  // advances through the surfaces in order.
  ///
  /// -----------------------------------------------
  /* result.nav_surfaces
     = result.start_layer->getCompatibleSurfaces(nav_par,
                                                   forward,
                                                   true,
                                                   collectSensitive,
                                                   collectMaterial,
                                                   collectPassive,
                                                   navigationLevel,
                                                   cache.start_surface,
                                                   cache.target_surface);
 if (result.nav_surfaces.size()) {
   NAVLOG(cache,
          result,
          result.nav_surfaces.size() << " surface candidates found.");
   // set the iterator
   result.nav_surface_iter = result.nav_surfaces.begin();
   // potentially switch the first one
   if (onStartSurfaceCheck
       && std::abs(result.nav_surface_iter->intersection.pathLength)
         < s_onSurfaceTolerance  ) {
       //
       NAVLOG(cache,
              result,
              "This is on the start surface, skipping it.");
       // increase the iterator, but bail out if it's the last one
       ++result.nav_surface_iter;
       if (result.nav_surface_iter == result.nav_surfaces.end()){
           result.nav_surfaces.clear();
           result.nav_surface_iter = result.nav_surfaces.end();
           return false;
       }
   }
   // update the navigation step size before you return
   update_step(cache, result, result.nav_surface_iter);
   // we still return
   return true;
  }
  // no returing to the stepper
  return false;
 }
  */

  /// Resolve the surfaces of this layer, if not the start layer
  ///
  /// @tparam cache_t is the cache type
  /// @tparam result_t is the cache type
  ///
  /// @return whether you found surfaces or not
  template <typename cache_t, typename result_t>
  bool
  resolve_surfaces(const NavigationParameters& nav_par,
                   cache_t&                    cache,
                   result_t&                   result) const
  {
    // get the layer and layer surface
    auto layer_surface = result.nav_layer_iter->representation;
    auto nav_layer     = result.nav_layer_iter->object;
    // are we on the start layer
    bool on_start      = (nav_layer == result.start_layer);
    auto start_surface = on_start ? cache.start_surface : layer_surface;
    // get the surfaces
    // @todo: could symmetrise with decompose() method
    result.nav_surfaces
        = nav_layer->getCompatibleSurfaces(nav_par,
                                           forward,
                                           true,
                                           collectSensitive,
                                           collectMaterial,
                                           collectPassive,
                                           navigationLevel,
                                           start_surface,
                                           cache.target_surface);
    // the number of layer candidates
    if (result.nav_surfaces.size()) {
      NAVLOG(cache,
             result,
             result.nav_surfaces.size() << " surface candidates found.");
      // set the iterator
      result.nav_surface_iter = result.nav_surfaces.begin();
      // update the navigation step size before you return
      update_step(cache, result, result.nav_surface_iter);
      return true;
    }
    result.nav_surface_iter = result.nav_surfaces.end();
    NAVLOG(cache, result, "No surface candidates found, switching layer.");
    return false;
  }

  // Loop over surface candidates here:
  //  - if an intersect is  valid but not yet reached
  //    then return with updated step size
  //  - if an intersect is not valid, switch to next
  //
  // return (bool) triggers a return to the stepper
  template <typename cache_t, typename result_t>
  bool
  handle_surfaces(const NavigationParameters& nav_par,
                  cache_t&                    cache,
                  result_t&                   result) const
  {
    // no surfaces, do not return to stepper
    if (!result.nav_surfaces.size()) return false;
    // loop over the navigation surfaces
    while (result.nav_surface_iter != result.nav_surfaces.end()) {
      // take the surface
      auto surface = result.nav_surface_iter->object;
      // case (s-a): we are at a surface
      //           only ask if you hadn't just come from a valid surface
      //
      // @todo : we should have a good idea whether this check should be done
      // @todo : add tolerance
      //
      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      if (surface->isOnSurface(nav_par.position(), true)) {
        NAVLOG(cache, result, "Surface successfully hit, storing it.");
        // the surface will only appear due to correct
        // collect(Property) flag
        cache.current_surface = surface;
        // break if the surface is the target surface
        if (surface == cache.target_surface) {
          NAVLOG(cache, result, "This was the target surface. Done.");
          return true;
        }
        // switch to the next candidate
        ++result.nav_surface_iter;
      }
      // screen output how much is left to try
      NAVLOG(cache,
             result,
             std::distance(result.nav_surface_iter, result.nav_surfaces.end())
                 << " out of "
                 << result.nav_surfaces.size()
                 << " surfaces remain to try.");
      // case (s-b) : update step estimation to the new surface
      if (result.nav_surface_iter != result.nav_surfaces.end()) {
        surface                = result.nav_surface_iter->object;
        auto surface_intersect = surface->intersectionEstimate(
            nav_par.position(), nav_par.momentum(), true, false);
        double surface_distance = surface_intersect.pathLength;
        if (!surface_intersect) {
          NAVLOG(
              cache, result, "Surface intersection is not valid, skipping it.");
          ++result.nav_surface_iter;
          continue;
        } else {
          cache.step_size.update(cache.nav_dir * surface_distance,
                                 cstep::actor);
          NAVLOG(cache,
                 result,
                 "Navigation step_size towards surface updated to "
                     << cache.step_size.toString());
          return true;
        }
      }
    }
    // case (s-c) : reached the end of the surface iteration
    if (result.nav_surface_iter == result.nav_surfaces.end()) {
      NAVLOG(cache, result, "Last surface on layer reached, switching layer.");
      // first clear the surface cache
      result.nav_surfaces.clear();
      result.nav_surface_iter = result.nav_surfaces.end();
      // now switch to the next layer
      ++result.nav_layer_iter;
    }
    // do not return to the stepper
    return false;
  }

  /// This method updates the constrained step size
  /// @tparam cache_t is the cache type
  /// @tparam result_t is the cache type
  /// @tparam type_t is the cache type
  template <typename cache_t, typename result_t, typename type_t>
  void
  update_step(cache_t& cache, result_t& result, type_t& type) const
  {
    //  update the step
    double ustep
        = initialStepFactor * cache.nav_dir * type->intersection.pathLength;
    cache.step_size.update(ustep, cstep::actor);
    NAVLOG(cache,
           result,
           "Navigation step_size updated to " << cache.step_size.toString());
  }
};
}
#endif  // end of namespace Acts
