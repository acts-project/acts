// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_NAVIGATOR_H
#define ACTS_NAVIGATOR_H

#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"

#include <sstream>

#ifndef NAVIGATOROUTPUTS
#define NAVIGATOROUTPUTS
#define VLOG(result,dump)                                \
  if (debug){                                            \
    std::string volumeName = result.current_volume ?     \
      result.current_volume->volumeName() : "No Volume"; \
    std::stringstream dstream;                           \
    dstream << "[ " << std::setw(30)                     \
            << result.current_volume->volumeName()       \
            << " ] " <<  std::setw(50)                   \
            << dump << '\n';                             \
    std::cout << dstream.str();                          \
    result.debug_string += dstream.str();}                                                       
#endif

namespace Acts {
      
  namespace extrapolation {

  /// Struct to mimmick track parameters
  /// @todo harmonize to eventual future update of 
  /// TrackingVolume::layerCandidatesOrdered() 
  /// TrackingVolume::boundarySurfacesOrdered()
  struct NavParameters {
    
    /// Position
    Vector3D pos;
    /// Direction
    Vector3D dir;
    
    /// Access method to satisify TrackingVolume interface
    const Vector3D& position() const 
      {return pos;}
    
    /// Access method to satisify TrackingVolume interface
    const Vector3D& momentum() const 
      {return dir;}
    
    /// Constructor 
    NavParameters(const Vector3D& p, const Vector3D& d)
       : pos(p), dir(d) {} 
    
  };

  typedef std::vector<SurfaceIntersection>  NavSurfaces;
  typedef std::vector<LayerIntersection<NavParameters> > NavLayers;
  typedef std::vector< BoundaryIntersection<NavParameters> > NavBoundaries;  

  /// Navigator struct
  ///
  /// This is an Actor to be added to the ActorList in order to navigate
  /// through the static tracking geometry setup.
  ///      
  /// The current navigation stage is cached in the result
  /// and updated when ever necessary. If any surface in the extrapolation
  /// flow is hit, it is set to the cache, such that other actors can 
  /// deal wit it. 
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
    bool collectMaterial  = true;
    /// stop at every surface regardless what it is 
    bool collectPassive   = true;
    
    /// store the debug message
    bool debug            = false;
      
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
      NavLayers nav_layers = {};
      NavLayers::const_iterator nav_layer_iter = nav_layers.end();
      
      /// Navigation on volume level
      /// Navigation cache: the current volume
      const TrackingVolume* current_volume = nullptr;
      /// Navigation cache: the target volume
      const TrackingVolume* target_volume  = nullptr; 
      
      NavBoundaries nav_boundaries = {}; 
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
      NavParameters nav_par(cache.position(),cache.direction());
      auto nav_dir = PropDirection(1); // @todo !needs to come from the cache
      
      // Navigation initialisation -------------------------------------------------
      // check for navigation initialisation & do it if necessary
      if (!result.current_volume){
          result.current_volume = trackingGeometry->lowestTrackingVolume(cache.pos);
          if (result.current_volume){
            result.nav_layers = result.current_volume->layerCandidatesOrdered(nullptr, 
                                                                              nullptr,
                                                                              nav_par,
                                                                              nav_dir,
                                                                              true,
                                                                              collectSensitive,
                                                                              collectMaterial,
                                                                              collectPassive);
            VLOG(result,"Initialised with " << result.nav_layers.size() << " layer candidates");                                                                  
            result.nav_layer_iter = result.nav_layers.begin();
            if (result.nav_layers.size()){
              // update the step size to the first layer
              cache.step_size = result.nav_layer_iter->intersection.pathLength;
              VLOG(result, "Initial step size towards layer updated to " << cache.step_size);
              return;
           }                                                                          
        }                                                                      
      }
       
      // Navigation through surfaces -----------------------------------------------
      // check if there are surfaces to be processed
      if (result.nav_surface_iter != result.nav_surfaces.end()){
        auto surface = result.nav_surface_iter->object;
        // @todo : we should have a good idea if this check is already to be done
        // @todo : add tolerance
        if (surface->isOnSurface(cache.pos, true)){
          VLOG(result, "Surface successfully, storing it.");
          // the surface will only appear due to correct
          // collect(Property) flag
          cache.current_surface = surface;
          // swith to the next candidate
          ++result.nav_surface_iter;
        }
        // check if we still have a candidate here
        if (result.nav_surface_iter != result.nav_surfaces.end()){
          // update to the new surface
          /// @todo: in straight line case, we can re-use
          surface = result.nav_surface_iter->object;
          auto surface_intersect = surface->intersectionEstimate(cache.pos,
                                                                 direction * cache.dir,
                                                                 true,
                                                                 false);
          double surface_distance = surface_intersect.pathLength;
          if (!surface_intersect){
            VLOG(result, "Surface intersection is not valid, skipping it.");
            ++result.nav_surface_iter;
          } else if (std::abs(cache.step_size) > std::abs(surface_distance)) {
            cache.step_size = surface_distance;     
            VLOG(result, "Step size towards surface updated to " << cache.step_size);
            return;
          }
        }
        
        // the surface iterator may have been updated
        if (result.nav_surface_iter == result.nav_surfaces.end()){
          result.nav_surfaces.clear();
          result.nav_surface_iter = result.nav_surfaces.end();
          VLOG(result, "Last surface hit, switching layer.");
          ++result.nav_layer_iter;
          if (result.nav_layer_iter != result.nav_layers.end()){
            // adjust the next steo size and return
            auto layer_surface = result.nav_layer_iter->representation;
            auto layer_intersect
                 = layer_surface->intersectionEstimate(cache.pos,
                                                 direction * cache.dir,
                                                 true,
                                                 false);
            double layer_distance = layer_intersect.pathLength;   
            cache.step_size = layer_distance;
            VLOG(result, "Initial step size towards layer updated to " << cache.step_size);  
            return;
          }
        }
      } 
              
      // Navigation through layers -------------------------------------------------
      // check if there are layers to be processed
      if (result.nav_layer_iter != result.nav_layers.end()){
        auto layer_surface = result.nav_layer_iter->representation;
        // check if we are on already on surface: we should have a good idea when to ask
        // @todo add tolerance
        if (layer_surface->isOnSurface(cache.pos, true)){
            VLOG(result, "Layer reached, prepare surfaces to be processed.");
            // collect if configured to do so
            if ( (layer_surface->associatedMaterial() && collectMaterial) || collectPassive )
              cache.current_surface = layer_surface;
            // now get the surfaces from the layer
            auto nav_layer = result.nav_layer_iter->object;
            // check if this layer has to be resolved
            if (nav_layer->resolve(collectSensitive,collectMaterial,collectPassive)){
              result.nav_surfaces = nav_layer->getCompatibleSurfaces(nav_par,
                                                                     nav_dir,
                                                                     true,
                                                                     collectSensitive,
                                                                     collectMaterial,
                                                                     collectPassive,
                                                                     navigationLevel,
                                                                     layer_surface);
             result.nav_surface_iter = result.nav_surfaces.begin();
             // no compatible surface means switch to next layer
             if (!result.nav_surfaces.size()){
               VLOG(result, "No compatible surfaces on this layer, skipping it.");
               ++result.nav_layer_iter;
             } else {
               VLOG(result, result.nav_surfaces.size() << " compatible surfaces to try found.");
               // update the step size towards the first one
               result.nav_surface_iter = result.nav_surfaces.begin();
               cache.step_size = result.nav_surface_iter->intersection.pathLength;
               VLOG(result, "Initial step size towards surface updated to " << cache.step_size);  
               return;
             }
          }
        }
        // run the step estimation for the (next) layer 
        if (result.nav_layer_iter != result.nav_layers.end()){
          layer_surface = result.nav_layer_iter->representation;
          auto layer_intersect
               = layer_surface->intersectionEstimate(cache.pos,
                                               direction * cache.dir,
                                               true,
                                               false);
          double layer_distance = layer_intersect.pathLength;   
          // check if the intersect is invalid
          if (!layer_intersect) {
            VLOG(result, "Layer intersection not valid, skipping it.");
            ++result.nav_layer_iter;
          } else if (std::abs(cache.step_size) > std::abs(layer_distance)) {
            cache.step_size = layer_distance;
            VLOG(result, "Step size towards layer updated to " << cache.step_size);
          } 
        }
        // the layer may have chagned, check if we are at the end of it
        if (result.nav_layer_iter == result.nav_layers.end()){
           // clear the navigation layers of the last volume
           result.nav_layers.clear();
           result.nav_layer_iter = result.nav_layers.end();
           result.nav_surfaces.clear();
           result.nav_surface_iter = result.nav_surfaces.end();
           VLOG(result, "Last layer reached in this volume, switching.");
           result.nav_boundaries
             = result.current_volume->boundarySurfacesOrdered(nav_par,nav_dir);
           result.nav_boundary_iter = result.nav_boundaries.begin();
           VLOG(result, result.nav_boundaries.size() << " boundaries provided.");
           // we can update the cache size here and return
           cache.step_size = result.nav_boundary_iter->intersection.pathLength;
           VLOG(result, "Step size towards boundary updated to " << cache.step_size);
           return;
        }      
      }                                  
            
      // Navigation through volumes ---------------------------------------
      // navigation boundaries to work off 
      if (result.nav_boundary_iter != result.nav_boundaries.end()){
        auto boundary_surface = result.nav_boundary_iter->representation;
        // check if we are on already in this step @todo add tolerance
        if (boundary_surface->isOnSurface(cache.pos, true)){
            VLOG(result, "Boundary surface reached, prepare volume switch.");
            // get the actual boundary for the navigation & the next volume
            auto boundary = result.nav_boundary_iter->object;
            result.current_volume = 
              boundary->attachedVolume(cache.pos, cache.dir, nav_dir);
            // store the boundary if configured to do so
            if ((boundary_surface->associatedMaterial() && collectMaterial)
              || collectPassive ) cache.current_surface = boundary_surface;
            // We still have a volume to deal with 
            if (result.current_volume){
            // get the layer candidates
              result.nav_layers
                = result.current_volume->layerCandidatesOrdered(nullptr, 
                                                                nullptr,
                                                                nav_par,
                                                                nav_dir,
                                                                true,
                                                                collectSensitive,
                                                                collectMaterial,
                                                                collectPassive);
              VLOG(result,"Initialised with " << result.nav_layers.size() << " layer candidates");  
              result.nav_layer_iter = result.nav_layers.begin();
              if (result.nav_layers.size()){
                // update the step size to the first layer
                cache.step_size = result.nav_layer_iter->intersection.pathLength;
                VLOG(result, "Initial step size towards layer updated to " << cache.step_size);
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
        auto boundary_intersect
            = boundary_surface->intersectionEstimate(cache.pos,
                                            direction * cache.dir,
                                            true,
                                            false);
        double boundary_distance = boundary_intersect.pathLength;    
        // test the next boundary if this one does not work
        if (!boundary_intersect) {
          VLOG(result, "Boundary intersection not valid, skipping it." << cache.step_size);                 
          ++result.nav_boundary_iter;
          // if there is no more boundary to leave, we're done
          if (result.nav_boundary_iter == result.nav_boundaries.end()){
            VLOG(result, "No more boundary surfaces to leave this volume. Navigation break."); 
            result.navigation_break = true;
            return;
          }
        } else if (std::abs(cache.step_size) > std::abs(boundary_distance)){
          cache.step_size = boundary_distance; 
          VLOG(result, "Stepsize towards boundary updated to " << cache.step_size);                 
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
