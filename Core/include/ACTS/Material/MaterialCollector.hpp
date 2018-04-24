// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MATERIALCOLLECTOR_H
#define ACTS_MATERIALCOLLECTOR_H

#include <sstream>
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/Surface.hpp"

#ifndef MATCOLLECTOR_DEBUG_OUTPUTS
#define MATCOLLECTOR_DEBUG_OUTPUTS
#define MATCLOG(cache, result, message)                                        \
  if (debug) {                                                                 \
    std::stringstream dstream;                                                 \
    dstream << "   " << std::setw(cache.debugPfxWidth);                        \
    dstream << "material collection"                                           \
            << " | ";                                                          \
    dstream << std::setw(cache.debugMsgWidth) << message << '\n';              \
    cache.debugString += dstream.str();                                        \
  }
#endif

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

  /// Screen output steering
  bool debug = false;

  /// Simple result struct to be returned
  /// It collects the indivdual
  struct this_result
  {
    std::vector<const MaterialHit> collected;
    double                         materialInX0 = 0.;
    double                         materialInL0 = 0.;
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
    // if we are on target, everything should have been done
    if (cache.targetReached) return;

    if (cache.currentSurface)
      MATCLOG(cache,
              result,
              "Material check on surface "
                  << cache.currentSurface->geoID().toString());

    // a current surface has been already assigned by the navigator
    if (cache.currentSurface && cache.currentSurface->associatedMaterial()) {

      // get the material propertices and only continue
      const MaterialProperties* mProperties
          = cache.currentSurface->associatedMaterial()->material(
              cache.position());
      if (mProperties) {
        // check if you have a factor for pre/post/full update to do
        double prepofu = 1.;
        if (cache.startSurface == cache.currentSurface) {
          MATCLOG(cache, result, "Update on start surface: post-update mode.");
          prepofu = cache.currentSurface->associatedMaterial()->factor(
              cache.navDir, postUpdate);
        } else if (cache.targetSurface == cache.currentSurface) {
          MATCLOG(cache, result, "Update on target surface: pre-update mode.");
          prepofu = cache.currentSurface->associatedMaterial()->factor(
              cache.navDir, preUpdate);
        } else
          MATCLOG(cache, result, "Update while pass through: full mode.");

        if (prepofu == 0.) {
          MATCLOG(cache, result, "Pre/Post factor set material to zero.");
          return;
        }

        MATCLOG(cache, result, "Material properties found for this surface.");
        // the path correction from the surface intersection
        double pCorrection = prepofu
            * cache.currentSurface->pathCorrection(cache.position(),
                                                   cache.direction());
        // the full material
        result.materialInX0 += pCorrection * mProperties->thicknessInX0();
        result.materialInL0 += pCorrection * mProperties->thicknessInL0();

        MATCLOG(cache, result, "t/X0 increased to " << result.materialInX0);
        MATCLOG(cache, result, "t/L0 increased to " << result.materialInL0);

        // if configured, record the individual material hits
        if (detailedCollection) {
          // create for recording
          MaterialHit material_hit;
          material_hit.surface   = cache.currentSurface;
          material_hit.position  = cache.position();
          material_hit.direction = cache.direction();
          // get the material & path length
          material_hit.material   = mProperties->material();
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
