
// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MATERIALCOLLECTOR_H
#define ACTS_MATERIALCOLLECTOR_H

#include <sstream>
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"

namespace Acts {

/// The information to written out per hit surface
struct MaterialHit
{
  const Surface* surface = nullptr;
  Vector3D       position;
  Vector3D       direction;
  Material       material;
  double         pathLength;
};

/// A Material Collector struct struct
struct MaterialCollector
{
  
  /// In the detailed collection mode the material
  /// per surface is collected, otherwise only the total
  /// pathlength in X0 or L0 are recorded
  bool detailedCollection = false;

  /// Simple result struct to be returned
  /// It collects the indivdual 
  struct this_result
  {
    std::vector<const MaterialHit> collected = {};
    double materialInX0 = 0.;
    double materialInL0 = 0.;
  };

  typedef this_result result_type;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the cache has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam cache_t is the type of Stepper cache
  ///
  /// @param cache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename cache_t>
  void
  operator()(cache_t& cache, result_type& result) const
  {
    // a current surface has been already assigned by the navigator
    if (cache.current_surface && cache.current_surface->associatedMaterial()){
      // get the material propertices and only continue
      const MaterialProperties* mProperties =
      cache.current_surface->associatedMaterial()->material(cache.position()); 
      if (mProperties){
        // the path correction from the surface intersection
        double pCorrection = 
          cache.current_surface->pathCorrection(cache.position(),
                                                cache.direction());
        // the full material
        materialInX0 += pCorrection*thicknessInX0();
        materialInL0 += pCorrection*thicknessInL0();
        // if configured, record the individual material hits
        if (detailedCollection){
          // create for recording
          MaterialHit material_hit;
          material_hit.surface   = cache.current_surface;
          material_hit.position  = cache.position();
          material_hit.direction = cache.direction();
          // get the material & path length
          material_hit.material = mProperties->material();
          material_hit.pathLength = pCorrection * mProperties->thickness();
          // save if in the result
          result.collected.push_back(material_hit);
        }
      }
    }
  }

  /// Pure observer interface
  /// - this does not apply to the surface collector
  template <typename cache_t>
  void
  operator()(cache_t& cache) const
  {
    (void)cache;
  }
};
}

#endif
