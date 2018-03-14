// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_NAVIGATOR_H
#define ACTS_NAVIGATOR_H

#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"

#include <sstream>

#ifndef NAVIGATOROUTPUTS
#define NAVIGATOROUTPUTS
#define VLOG(result, dump)                                                     \
  if (debug) {                                                                 \
    std::string vName                = "No Volume";                            \
    if (result.current_volume) vName = result.current_volume->volumeName();    \
    std::stringstream dstream;                                                 \
    dstream << "[ " << std::setw(30) << vName << " ] ";                        \
    dstream << std::setw(50) << dump << '\n';                                  \
    std::cout << dstream.str();                                                \
    result.debug_string += dstream.str();                                      \
  }
#endif

namespace Acts {

namespace extrapolation {

  /// Struct to mimmick track parameters
  /// @todo harmonize to eventual future update of
  /// TrackingVolume::layerCandidatesOrdered()
  /// TrackingVolume::boundarySurfacesOrdered()
  struct NavParameters
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
    NavParameters(const Vector3D& p, const Vector3D& d) : pos(p), dir(d) {}
  };

  typedef std::vector<SurfaceIntersection>                 NavSurfaces;
  typedef std::vector<LayerIntersection<NavParameters>>    NavLayers;
  typedef std::vector<BoundaryIntersection<NavParameters>> NavBoundaries;

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

    // the navigation level, see Layer for more information
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
    struct this_state
    {
      /// Navigation on surface level
      /// the vector of navigation surfaces to
      NavSurfaces nav_surfaces = {};
      /// the surface iterator
      NavSurfaces::const_iterator nav_surface_iter = nav_surfaces.end();

      /// Navigation on layer level
      /// the vector of navigation layer to
      NavLayers                 nav_layers     = {};
      NavLayers::const_iterator nav_layer_iter = nav_layers.end();

      /// Navigation on volume level
      /// Navigation cache: the current volume
      const TrackingVolume* current_volume = nullptr;
      /// Navigation cache: the target volume
      const TrackingVolume* target_volume = nullptr;

      NavBoundaries                 nav_boundaries    = {};
      NavBoundaries::const_iterator nav_boundary_iter = nav_boundaries.end();

      /// Debug string buffer
      std::string debug_string;

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

      int direction = 1;

      // fail if you have no tracking geometry
      assert(trackingGeometry != nullptr);

      // the navigator erases the current surfaces
      cache.current_surface = nullptr;

      // navigation parameters
      NavParameters nav_par(cache.position(), cache.direction());
      auto nav_dir = PropDirection(1);  // @todo !needs to come from the cache

      // Navigation initialisation
      // -------------------------------------------------
      // Vheck for navigation initialisation & do it if necessary.  This means
      // we do not have an active volume yet (since we just started).  We get
      // the innermost volume, and set up an ordered layer iterator and step
      // size toward the first layer.  The return prevent execution of the
      // subsequent logic, we want to make a step first.
      if (!result.current_volume) {
        result.current_volume
            = trackingGeometry->lowestTrackingVolume(cache.pos);
        if (result.current_volume) {
          result.nav_layers
              = result.current_volume->layerCandidatesOrdered(nullptr,
                                                              nullptr,
                                                              nav_par,
                                                              nav_dir,
                                                              true,
                                                              collectSensitive,
                                                              collectMaterial,
                                                              collectPassive);
          VLOG(result,
               "Initialised with " << result.nav_layers.size()
                                   << " layer candidates");
          result.nav_layer_iter = result.nav_layers.begin();
          if (result.nav_layers.size()) {
            // update the step size to the first layer
            cache.step_size = result.nav_layer_iter->intersection.pathLength;
            VLOG(result,
                 "Initial step size towards layer updated to "
                     << cache.step_size);
            return;
          }
        }
      }

      // Navigation through surfaces
      // -----------------------------------------------
      // If we are currently iterating over surfaces, we get the current one
      // pointed at by the iterator here.  This is the main routine that
      // advances through the surfaces in order.
      if (result.nav_surface_iter != result.nav_surfaces.end()) {
        auto surface = result.nav_surface_iter->object;
        // @todo : we should have a good idea if this check is already to be
        // done
        // @todo : add tolerance
        // If we are on the surface pointed at by the iterator, we can make
        // it the current one to pass it to the other actors.
        // Also make iterator point to the next one.
        if (surface->isOnSurface(cache.pos, true)) {
          VLOG(result, "Surface successfully, storing it.");
          // the surface will only appear due to correct
          // collect(Property) flag
          cache.current_surface = surface;
          // swith to the next candidate
          ++result.nav_surface_iter;
        }

        // Check if we still have a candidate here. If we previously found a
        // surface and advanced the iterator, we want to check here if there is
        // another one, and if so update the step size accordingly. If there is
        // a next surface but we don't intersect, we skip it.  If we intersect,
        // we terminate this call by returning.
        if (result.nav_surface_iter != result.nav_surfaces.end()) {
          // update to the new surface
          /// @todo: in straight line case, we can re-use
          surface                = result.nav_surface_iter->object;
          auto surface_intersect = surface->intersectionEstimate(
              cache.pos, direction * cache.dir, true, false);
          double surface_distance = surface_intersect.pathLength;
          if (!surface_intersect) {
            VLOG(result, "Surface intersection is not valid, skipping it.");
            ++result.nav_surface_iter;
          } else if (std::abs(cache.step_size) > std::abs(surface_distance)) {
            cache.step_size = surface_distance;
            VLOG(result,
                 "Step size towards surface updated to " << cache.step_size);
            return;
          }
        }

        // the surface iterator may have been updated
        // If we skipped and are now out of surfaces, we
        // switch to the next layer, and set the step size accordingly
        if (result.nav_surface_iter == result.nav_surfaces.end()) {
          result.nav_surfaces.clear();
          result.nav_surface_iter = result.nav_surfaces.end();
          VLOG(result, "Last surface hit, switching layer.");
          ++result.nav_layer_iter;
          if (result.nav_layer_iter != result.nav_layers.end()) {
            // adjust the next steo size and return
            auto layer_surface   = result.nav_layer_iter->representation;
            auto layer_intersect = layer_surface->intersectionEstimate(
                cache.pos, direction * cache.dir, true, false);
            double layer_distance = layer_intersect.pathLength;
            cache.step_size       = layer_distance;
            VLOG(result,
                 "Initial step size towards layer updated to "
                     << cache.step_size);
            return;
          }
        }
      }

      // Navigation through layers
      // -------------------------------------------------
      // We are now trying to advance to the next layer (with surfaces)
      // Check if we are on the representing surface of the layer pointed
      // at by nav_layer_iter. If so, we unpack the compatible surfaces
      // (determined by straight line intersect), and set up the iterator
      // so that the next call to operator() will enter the surface
      // check mode above. If no surfaces are found, we skip the layer.
      // If we unpack a surface, the step size is set to the path length
      // to the first surface, as determined by straight line intersect.
      if (result.nav_layer_iter != result.nav_layers.end()) {
        auto layer_surface = result.nav_layer_iter->representation;
        // check if we are on already on surface: we should have a good idea
        // when to ask
        // @todo add tolerance
        if (layer_surface->isOnSurface(cache.pos, true)) {
          VLOG(result, "Layer reached, prepare surfaces to be processed.");
          // collect if configured to do so
          // if we don't set this here, the rest of the actors
          // will not process this surface
          if ((layer_surface->associatedMaterial() && collectMaterial)
              || collectPassive)
            cache.current_surface = layer_surface;
          // now get the surfaces from the layer
          auto nav_layer = result.nav_layer_iter->object;
          // check if this layer has to be resolved
          if (nav_layer->resolve(
                  collectSensitive, collectMaterial, collectPassive)) {
            result.nav_surfaces
                = nav_layer->getCompatibleSurfaces(nav_par,
                                                   nav_dir,
                                                   true,
                                                   collectSensitive,
                                                   collectMaterial,
                                                   collectPassive,
                                                   navigationLevel,
                                                   layer_surface);
            result.nav_surface_iter = result.nav_surfaces.begin();
            // no compatible surface means switch to next layer
            if (!result.nav_surfaces.size()) {
              VLOG(result,
                   "No compatible surfaces on this layer, skipping it.");
              ++result.nav_layer_iter;
            } else {
              VLOG(result,
                   result.nav_surfaces.size()
                       << " compatible surfaces to try found.");
              // update the step size towards the first one
              result.nav_surface_iter = result.nav_surfaces.begin();
              cache.step_size
                  = result.nav_surface_iter->intersection.pathLength;
              VLOG(result,
                   "Initial step size towards surface updated to "
                       << cache.step_size);
              return;
            }
          }
        }

        // If we skipped the layer above, and we have a next one,
        // we check if straight line intersects, and if so, set the
        // max step size. If not, we skip the layer.
        if (result.nav_layer_iter != result.nav_layers.end()) {
          layer_surface        = result.nav_layer_iter->representation;
          auto layer_intersect = layer_surface->intersectionEstimate(
              cache.pos, direction * cache.dir, true, false);
          double layer_distance = layer_intersect.pathLength;
          // check if the intersect is invalid
          if (!layer_intersect) {
            VLOG(result, "Layer intersection not valid, skipping it.");
            ++result.nav_layer_iter;
          } else if (std::abs(cache.step_size) > std::abs(layer_distance)) {
            cache.step_size = layer_distance;
            VLOG(result,
                 "Step size towards layer updated to " << cache.step_size);
          }
        }

        // If we skipped above and the layer pointed at by the iterator
        // is the last one, we are done with this volume. We clear the
        // layer and surface iterators, and get the current volumes boundaries
        // to determine which one is the next volume.
        // We set up the boundary iterator and initialize the step size
        // to the straigh line path length to the first boundary surface
        if (result.nav_layer_iter == result.nav_layers.end()) {
          // clear the navigation layers of the last volume
          result.nav_layers.clear();
          result.nav_layer_iter = result.nav_layers.end();
          result.nav_surfaces.clear();
          result.nav_surface_iter = result.nav_surfaces.end();
          VLOG(result, "Last layer reached in this volume, switching.");
          result.nav_boundaries
              = result.current_volume->boundarySurfacesOrdered(nav_par,
                                                               nav_dir);
          result.nav_boundary_iter = result.nav_boundaries.begin();
          VLOG(result, result.nav_boundaries.size() << " boundaries provided.");
          // we can update the cache size here and return
          cache.step_size = result.nav_boundary_iter->intersection.pathLength;
          VLOG(result,
               "Step size towards boundary updated to " << cache.step_size);
          return;
        }
      }

      // Navigation through volumes
      // -------------------------------------------------
      // This is the boundary check routine. If the code above set up the
      // boundary surface iterator, we advance through them here. If we are on
      // the boundary surface, we set the current surface to the boundary
      // surface, and get the volume pointed at by the boundary surface.  Next
      // we unpack the layers from that volume. If the volume contains layers
      // we set the step size to the straight line path length to the first
      // layer.  If we don't find a next volume, the navigation_break indicator
      // is set.  This ends the navigation. Finally, the boundary iterator is
      // cleared, so that the subsequent call goes back to the layer iteration
      // logic.
      //
      // If we are not on the current boundary surface, we try the next one.
      // The iterator is advanced and the step size is set. If no straight line
      // intersect if found, the boundary surface is skipped.  If we are out of
      // boundary surfaces, the navigation is terminated.
      if (result.nav_boundary_iter != result.nav_boundaries.end()) {
        auto boundary_surface = result.nav_boundary_iter->representation;
        // check if we are on already in this step @todo add tolerance
        if (boundary_surface->isOnSurface(cache.pos, true)) {
          VLOG(result, "Boundary surface reached, prepare volume switch.");
          // get the actual boundary for the navigation & the next volume
          auto boundary = result.nav_boundary_iter->object;
          result.current_volume
              = boundary->attachedVolume(cache.pos, cache.dir, nav_dir);
          // store the boundary if configured to do so
          // If we don't set it here, it wont be processed by the other actors.
          if ((boundary_surface->associatedMaterial() && collectMaterial)
              || collectPassive)
            cache.current_surface = boundary_surface;
          // We still have a volume to deal with
          if (result.current_volume) {
            // get the layer candidates
            result.nav_layers = result.current_volume->layerCandidatesOrdered(
                nullptr,
                nullptr,
                nav_par,
                nav_dir,
                true,
                collectSensitive,
                collectMaterial,
                collectPassive);
            VLOG(result,
                 "Initialised with " << result.nav_layers.size()
                                     << " layer candidates");
            result.nav_layer_iter = result.nav_layers.begin();
            if (result.nav_layers.size()) {
              // update the step size to the first layer
              cache.step_size = result.nav_layer_iter->intersection.pathLength;
              VLOG(result,
                   "Initial step size towards layer updated to "
                       << cache.step_size);
            }
          } else {
            VLOG(result, "No next volume found. Navigation break.");
            result.navigation_break = true;
            return;
          }
          // and we can invalidate the boundary surfaces and return
          result.nav_boundaries.clear();
          result.nav_boundary_iter = result.nav_boundaries.end();
          return;
        }
        // intersect the boundary
        auto boundary_intersect = boundary_surface->intersectionEstimate(
            cache.pos, direction * cache.dir, true, false);
        double boundary_distance = boundary_intersect.pathLength;
        // test the next boundary if this one does not work
        if (!boundary_intersect) {
          VLOG(result,
               "Boundary intersection not valid, skipping it."
                   << cache.step_size);
          ++result.nav_boundary_iter;
          // if there is no more boundary to leave, we're done
          if (result.nav_boundary_iter == result.nav_boundaries.end()) {
            VLOG(result,
                 "No more boundary surfaces to leave this volume. "
                 "Navigation break.");
            result.navigation_break = true;
            return;
          }
        } else if (std::abs(cache.step_size) > std::abs(boundary_distance)) {
          cache.step_size = boundary_distance;
          VLOG(result,
               "Stepsize towards boundary updated to " << cache.step_size);
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
}

#endif
