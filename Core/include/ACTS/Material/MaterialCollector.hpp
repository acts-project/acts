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

namespace Acts {

/// The information to be writtern out per hit surface
struct MaterialHit
{
  const Surface* surface = nullptr;
  Vector3D       position;
  Vector3D       direction;
  Material       material;
  double         pathLength;
};

/// A Material Collector struct
struct MaterialCollector
{

  /// In the detailed collection mode the material
  /// per surface is collected, otherwise only the total
  /// pathlength in X0 or L0 are recorded
  bool detailedCollection = false;

  /// Enable debug output writing
  bool debug = false;

  /// Simple result struct to be returned
  ///
  /// Result of the material collection process
  /// It collects the overall X0 and L0 path lengths,
  /// and optionally a detailed per-material breakdown
  struct this_result
  {
    std::vector<MaterialHit> collected;
    double                   materialInX0 = 0.;
    double                   materialInL0 = 0.;
  };

  typedef this_result result_type;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the state has a current surface,
  /// in which case the action is performed:
  /// - it records the surface given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_state_t is the type of Stepper state
  ///
  /// @param propState is the mutable propagator state object
  /// @param stepState is the mutable stepper state object
  /// @param result is the result object to be filled
  template <typename propagator_state_t, typename stepper_state_t>
  void
  operator()(propagator_state_t& propState,
             stepper_state_t&    stepState,
             result_type&        result) const
  {
    // if we are on target, everything should have been done
    if (propState.targetReached) return;

    if (propState.currentSurface) {
      debugLog(propState, [&] {
        std::stringstream dstream;
        dstream << "Material check on surface ";
        dstream << propState.currentSurface->geoID().toString();
        return dstream.str();
      });
    }

    // a current surface has been already assigned by the navigator
    if (propState.currentSurface
        && propState.currentSurface->associatedMaterial()) {

      // get the material propertices and only continue
      const MaterialProperties* mProperties
          = propState.currentSurface->associatedMaterial()->material(
              stepState.position());
      if (mProperties) {
        // pre/post/full update
        double prepofu = 1.;
        if (propState.startSurface == propState.currentSurface) {
          debugLog(propState, [&] {
            return std::string("Update on start surface: post-update mode.");
          });
          prepofu = propState.currentSurface->associatedMaterial()->factor(
              stepState.navDir, postUpdate);
        } else if (propState.targetSurface == propState.currentSurface) {
          debugLog(propState, [&] {
            return std::string("Update on target surface: pre-update mode");
          });
          prepofu = propState.currentSurface->associatedMaterial()->factor(
              stepState.navDir, preUpdate);
        } else {
          debugLog(propState, [&] {
            return std::string("Update while pass through: full mode.");
          });
        }

        // the pre/post factor has been applied
        // now check if there's still something to do
        if (prepofu == 0.) {
          debugLog(propState, [&] {
            return std::string("Pre/Post factor set material to zero.");
          });
          return;
        }
        // more debugging output to the screen
        debugLog(propState, [&] {
          return std::string("Material properties found for this surface.");
        });

        // the path correction from the surface intersection
        double pCorrection = prepofu
            * propState.currentSurface->pathCorrection(stepState.position(),
                                                       stepState.direction());
        // the full material
        result.materialInX0 += pCorrection * mProperties->thicknessInX0();
        result.materialInL0 += pCorrection * mProperties->thicknessInL0();

        debugLog(propState, [&] {
          std::stringstream dstream;
          dstream << "t/X0 (t/L0) increased to ";
          dstream << result.materialInX0 << " (";
          dstream << result.materialInL0 << " )";
          return dstream.str();
        });

        // if configured, record the individual material hits
        if (detailedCollection) {
          // create for recording
          MaterialHit material_hit;
          material_hit.surface   = propState.currentSurface;
          material_hit.position  = stepState.position();
          material_hit.direction = stepState.direction();
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
  template <typename propagator_state_t, typename stepper_state_t>
  void
  operator()(propagator_state_t&, stepper_state_t&) const
  {
  }

private:
  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param propState the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&          propState,
           std::function<std::string()> logAction) const
  {
    if (debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(propState.options.debugPfxWidth);
      dstream << "material collector"
              << " | ";
      dstream << std::setw(propState.options.debugMsgWidth) << logAction()
              << '\n';
      propState.options.debugString += dstream.str();
    }
  }
};
}

#endif
