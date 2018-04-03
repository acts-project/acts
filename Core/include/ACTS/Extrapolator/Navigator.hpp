// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_NAVIGATOR_H
#define ACTS_NAVIGATOR_H

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

typedef std::vector<SurfaceIntersection>                     NavigationSurfaces;
typedef std::vector<LayerIntersection<NavigationParameters>> NavigationLayers;
typedef std::vector<BoundaryIntersection<NavigationParameters>>
    NavigationBoundaries;

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
    /// the vector of navigation surfaces to
    NavigationSurfaces nav_surfaces = {};
    /// the surface iterator
    NavigationSurfaces::const_iterator nav_surface_iter = nav_surfaces.end();

    /// Navigation on layer level
    /// the vector of navigation layer to
    NavigationLayers                 nav_layers     = {};
    NavigationLayers::const_iterator nav_layer_iter = nav_layers.end();

    /// Navigation on volume level
    /// Navigation cache: the current volume and current layer
    const Layer*          current_layer  = nullptr;
    const TrackingVolume* current_volume = nullptr;
    /// Navigation cache: the target volume and target layer
    const Layer*          target_layer  = nullptr;
    const TrackingVolume* target_volume = nullptr;

    NavigationBoundaries                 nav_boundaries = {};
    NavigationBoundaries::const_iterator nav_boundary_iter
        = nav_boundaries.end();

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

    /// do nothing if the navigation stream was broken
    if (result.navigation_break) return;

    // fail if you have no tracking geometry
    assert(trackingGeometry != nullptr);

    // the navigator erases the current surfaces
    cache.current_surface = nullptr;

    // navigation parameters
    NavigationParameters nav_par(cache.position(), cache.direction());

    // Navigation initialisation
    //
    // This is only called once for every propagation/extrapolation
    //
    // -------------------------------------------------
    // Check for navigation initialisation & do it if necessary.  This means
    // we do not have an active volume yet (since we just started).  We get
    // the innermost volume, and set up an ordered layer iterator and step
    // size toward the first layer.  The return prevent execution of the
    // subsequent logic, we want to make a step first.
    if (!result.current_volume) {
      NAVLOG(cache, result, "Initializing start volume.");
      // Fast Navigation initialization:
      // short-cut through object association, saves navigation in the
      // geometry and volume tree search for the lowser volume
      if (cache.start_surface && cache.start_surface->associatedLayer()) {
        // assign the current layer and volume by association
        result.current_layer  = cache.start_surface->associatedLayer();
        result.current_volume = result.current_layer->trackingVolume();
      }
      // Fast Navigation initialization for target:
      if (cache.target_surface && cache.target_surface->associatedLayer()) {
        // assign the target volume and the target surface
        result.target_layer  = cache.target_surface->associatedLayer();
        result.target_volume = result.target_layer->trackingVolume();
      }
      // flag for on-layer propagation
      bool stay_on_layer = result.current_layer
          && (result.current_layer == result.target_layer);
      // check and (if necessary) spend the time to navigate to
      // the lowest trackign volume - not neccesary if you stay on a layer
      if (!stay_on_layer) {
        result.current_volume = result.current_volume
            ? result.current_volume
            : trackingGeometry->lowestTrackingVolume(cache.pos);
        // safety check for misconfigured geometry
        // skip if start layer is target layer
        if (result.current_volume) {
          // we use start and (if current_volume == target_volume)
          // als the target layer
          const Layer* tLayer = (result.current_volume == result.target_volume)
              ? result.target_layer
              : nullptr;
          // now request the navigation layers
          result.nav_layers = result.current_volume->layerCandidatesOrdered(
              result.current_layer,
              tLayer,
              nav_par,
              cache.nav_dir,
              true,
              collectSensitive,
              collectMaterial,
              collectPassive);
          NAVLOG(cache,
                 result,
                 "Initial Volume with " << result.nav_layers.size()
                                        << " layer candidates");
          result.nav_layer_iter = result.nav_layers.begin();
          if (result.nav_layers.size()) {
            // update the step size to the first layer
            double ustep = initialStepFactor * cache.nav_dir
                * result.nav_layer_iter->intersection.pathLength;
            cache.step_size.update(ustep, cstep::actor);
            NAVLOG(cache,
                   result,
                   "Navigation step_size towards first layer estimated as "
                       << ustep);
            return;
          }
        } else {
          // this is the do-nothing case
          NAVLOG(cache,
                 result,
                 "Navigation could not be intialised. Pure progation.");
          return;
        }
      } else {
        // only one layer is given to the navigation
        // all yo uneed is resolve start and end surface on the layer
        if (result.current_layer->resolve(
                collectSensitive, collectMaterial, collectPassive)) {
          result.nav_surfaces = result.current_layer->getCompatibleSurfaces(
              nav_par,
              cache.nav_dir,
              true,
              collectSensitive,
              collectMaterial,
              collectPassive,
              navigationLevel,
              cache.start_surface,
              cache.target_surface);
          // set the start iterator
          result.nav_surface_iter = result.nav_surfaces.begin();
        }
      }
    }

    bool surfaceSwitch = false;
    bool layerSwitch   = false;
    bool volumeSwitch  = false;

    // Navigation through surfaces
    // -----------------------------------------------
    // If we are currently iterating over surfaces, we get the current one
    // pointed at by the iterator here.  This is the main routine that
    // advances through the surfaces in order.
    //
    // Loop over surface candidates here:
    //  - if an intersect is  valid but not yet reached
    //    then return with updated step size
    //  - if an intersect is not valid, switch to next
    while (result.nav_surface_iter != result.nav_surfaces.end()) {

      auto surface = result.nav_surface_iter->object;
      // case (s-a): we are at a surface
      //           only ask if you hadn't just come from a valid surface
      //
      // @todo : we should have a good idea whether this check should be done
      // @todo : add tolerance
      //
      // If we are on the surface pointed at by the iterator, we can make
      // it the current one to pass it to the other actors
      if (!surfaceSwitch && surface->isOnSurface(cache.pos, true)) {
        NAVLOG(cache, result, "Surface successfully hit, storing it.");
        // the surface will only appear due to correct
        // collect(Property) flag
        cache.current_surface = surface;
        // switch to the next candidate
        ++result.nav_surface_iter;
        // remember that you have done a surface switch
        surfaceSwitch = true;
        // call continue
        continue;
      }
      // case (s-b) : update step estimation to the new surface
      // @todo: in straight line case, we can re-use
      surface                = result.nav_surface_iter->object;
      auto surface_intersect = surface->intersectionEstimate(
          cache.pos, cache.nav_dir * cache.dir, true, false);
      double surface_distance = surface_intersect.pathLength;
      if (!surface_intersect) {
        NAVLOG(
            cache, result, "Surface intersection is not valid, skipping it.");
        ++result.nav_surface_iter;
      } else {
        double ustep = cache.nav_dir * surface_distance;
        cache.step_size.update(ustep, cstep::actor);
        NAVLOG(cache,
               result,
               "Navigation step_size towards surface updated to "
                   << cache.step_size);
        return;
      }
    }

    // case (s-c) : reached the end of the surface iteration
    if (result.nav_surfaces.size()
        && result.nav_surface_iter == result.nav_surfaces.end()) {
      NAVLOG(cache, result, "Last surface on layer reached, switch to next.");
      // first clear the surface cache
      result.nav_surfaces.clear();
      result.nav_surface_iter = result.nav_surfaces.end();
      // now switch to the next layer
      ++result.nav_layer_iter;
      // remember that you have done a layer switch
      layerSwitch = true;
    }

    // Navigation through layers
    // -------------------------------------------------
    //
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
    while (result.nav_layer_iter != result.nav_layers.end()) {
      // the representing layer surface, e.g. the approach surface
      auto layer_surface = result.nav_layer_iter->representation;
      // case (l-a): we are at the representing layer surface
      //           only ask if you hadn't just done a layer switch
      //
      if (!layerSwitch && layer_surface->isOnSurface(cache.pos, true)) {
        // store the current surface in the cache
        cache.current_surface = layer_surface;
        // now get the surfaces from the layer
        auto nav_layer = result.nav_layer_iter->object;
        // check if this layer has to be resolved
        // and retrieve the compatible sdurfaces if available
        if (nav_layer->resolve(
                collectSensitive, collectMaterial, collectPassive)) {
          result.nav_surfaces
              = nav_layer->getCompatibleSurfaces(nav_par,
                                                 cache.nav_dir,
                                                 true,
                                                 collectSensitive,
                                                 collectMaterial,
                                                 collectPassive,
                                                 navigationLevel,
                                                 layer_surface);
          result.nav_surface_iter = result.nav_surfaces.begin();
          // no compatible surface means switch to next layer
          if (!result.nav_surfaces.size()) {
            NAVLOG(
                cache,
                result,
                "No compatible surfaces on this layer, switch to next layer.");
            ++result.nav_layer_iter;
            // break if you run into the last layer
            if (result.nav_layer_iter == result.nav_layers.end()) break;
          } else {
            NAVLOG(cache,
                   result,
                   result.nav_surfaces.size()
                       << " compatible surfaces found to try.");
            // update the step size towards the first one
            result.nav_surface_iter = result.nav_surfaces.begin();
            double ustep            = initialStepFactor * cache.nav_dir
                * result.nav_surface_iter->intersection.pathLength;
            cache.step_size.update(ustep, cstep::actor);
            NAVLOG(cache,
                   result,
                   "Navigation step_size towards surface estimated as "
                       << ustep);
            // @todo - here we could store that a surface search has
            // been done
            return;
          }
        }
      }
      // case (l-b) :  update step estimation to the layer representing surface
      auto layer_intersect = layer_surface->intersectionEstimate(
          cache.pos, cache.nav_dir * cache.dir, true, false);
      double layer_distance = layer_intersect.pathLength;
      // check if the intersect is invalid
      if (!layer_intersect) {
        NAVLOG(cache, result, "Layer intersection not valid, skipping it.");
        ++result.nav_layer_iter;
      } else {
        // update the navigation step size
        double ustep = cache.nav_dir * layer_distance;
        cache.step_size.update(ustep, cstep::actor);
        NAVLOG(cache,
               result,
               "Navigation step_size towards layer updated to " << ustep);
        return;
      }
    }

    // case (l-c) :  We reached the last layer
    //
    // If we skipped above and the layer pointed at by the iterator
    // is the last one, we are done with this volume. We clear the
    // layer and surface iterators, and get the current volumes boundaries
    // to determine which one is the next volume.
    // We set up the boundary iterator and initialize the step size
    // to the straigh line path length to the first boundary surface
    if (result.nav_layers.size()
        && result.nav_layer_iter == result.nav_layers.end()) {
      // clear the navigation layer cache
      result.nav_layers.clear();
      result.nav_layer_iter = result.nav_layers.end();
      // clear the surface layer cache (should not be necessary, check)
      result.nav_surfaces.clear();
      result.nav_surface_iter = result.nav_surfaces.end();
      NAVLOG(cache, result, "Last layer reached in this volume, switching.");
      // get the boundary surfaces and set the iterator
      result.nav_boundaries = result.current_volume->boundarySurfacesOrdered(
          nav_par, cache.nav_dir);
      result.nav_boundary_iter = result.nav_boundaries.begin();
      NAVLOG(cache,
             result,
             result.nav_boundaries.size() << " boundaries provided.");
      // remember that you have done a volume switch
      volumeSwitch = true;
    }

    // Navigation through volumes
    // -------------------------------------------------
    // This is the boundary check routine. If the code above set up the
    // boundary surface iterator, we advance through them here. If we are on
    // the boundary surface, we set the current surface to the boundary
    // surface, and get the volume pointed at by the boundary surface.  Next
    // we unpack the layers from that volume. If the volume contains layers
    // we set the step size to the straight line path length to the first
    // layer.  If we don't find a next volume, the navigation_break
    // indicator
    // is set.  This ends the navigation. Finally, the boundary iterator is
    // cleared, so that the subsequent call goes back to the layer iteration
    // logic.
    //
    // If we are not on the current boundary surface, we try the next one.
    // The iterator is advanced and the step size is set. If no straight
    // line
    // intersect if found, the boundary surface is skipped.  If we are out
    // of
    // boundary surfaces, the navigation is terminated.
    while (result.nav_boundary_iter != result.nav_boundaries.end()) {
      auto boundary_surface = result.nav_boundary_iter->representation;
      // case (v-a) : you are on the boundary surface
      //              only if you hadn't just done a volume switch
      // check if we are on already in this step
      // @todo add tolerance
      if (!volumeSwitch && boundary_surface->isOnSurface(cache.pos, true)) {
        NAVLOG(
            cache, result, "Boundary surface reached, prepare volume switch.");
        // get the actual boundary for the navigation & the next volume
        auto boundary = result.nav_boundary_iter->object;
        result.current_volume
            = boundary->attachedVolume(cache.pos, cache.dir, cache.nav_dir);
        // store the boundary for eventual actors to work on it
        cache.current_surface = boundary_surface;
        // We still have a volume to deal with
        if (result.current_volume) {
          // get the layer candidates
          result.nav_layers
              = result.current_volume->layerCandidatesOrdered(nullptr,
                                                              nullptr,
                                                              nav_par,
                                                              cache.nav_dir,
                                                              true,
                                                              collectSensitive,
                                                              collectMaterial,
                                                              collectPassive);
          NAVLOG(cache,
                 result,
                 "Next Volume with " << result.nav_layers.size()
                                     << " layer candidates");
          result.nav_layer_iter = result.nav_layers.begin();
          if (result.nav_layers.size()) {
            // update the navigation step size to the first layer
            double ustep = initialStepFactor * cache.nav_dir
                * result.nav_layer_iter->intersection.pathLength;
            cache.step_size.update(ustep, cstep::actor);
            NAVLOG(cache,
                   result,
                   "Navigation step_size towards first layer estimated as "
                       << ustep);
          }
        } else {
          NAVLOG(cache, result, "No next volume found. Navigation break.");
          result.navigation_break = true;
          cache.step_size.release(cstep::actor);
          return;
        }
        // and we can invalidate the boundary surfaces and return
        result.nav_boundaries.clear();
        result.nav_boundary_iter = result.nav_boundaries.end();
        return;
      }
      // case (v-b) : update the step size towards the boundary
      // by straignt line intersection
      auto boundary_intersect = boundary_surface->intersectionEstimate(
          cache.pos, cache.nav_dir * cache.dir, true, false);
      double boundary_distance = boundary_intersect.pathLength;
      // test the next boundary if this one does not work
      if (!boundary_intersect) {
        NAVLOG(cache, result, "Boundary intersection not valid, skipping it.");
        ++result.nav_boundary_iter;
        // if there is no more boundary to leave, we're done
        if (result.nav_boundary_iter == result.nav_boundaries.end()) {
          NAVLOG(cache,
                 result,
                 "No more boundary surfaces to leave this volume. "
                 "Navigation break.");
          result.navigation_break = true;
          return;
        }
      } else {
        double ustep = cache.nav_dir * boundary_distance;
        cache.step_size.update(ustep, cstep::actor);
        NAVLOG(cache,
               result,
               "Navigation step_size towards boundary updated to " << ustep);
        return;
      }
    }
  }

  /// Pure observer interface
  /// This does not apply to the navigator
  template <typename cache_t>
  void
  operator()(cache_t& cache) const
  {
    (void)cache;
  }
};
}

#endif
