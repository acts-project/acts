// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"

namespace Acts {

  namespace detail {
  
  using cstep = detail::ConstrainedStep;  
    
  /// Update surface status - Single component
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param surface [in] The surface provided 
  /// @param bcheck [in] The boundary check for this status update
    template <typename stepper_t>
    Acts::Intersection::Status 
    updateSurfaceStatus_sc(const stepper_t& stepper,
                          typename stepper_t::State& state, 
                          const Surface& surface,
                          const BoundaryCheck& bcheck){
  
      // Now intersect (should exclude punch-through)
      auto sIntersection = surface.intersect(state.geoContext, 
                                             stepper.position(state),
                                             state.navDir*stepper.direction(state),
                                             bcheck);
                                             
      // The intersection is on surface already                                           
      if (sIntersection.intersection.status == Intersection::Status::onSurface){
        // Release navigation step size
        state.stepSize.release(cstep::actor);
        return Intersection::Status::onSurface)
      } else if (sIntersection.intersection or sIntersection.alternative){
        // Path and overstep limit checking
        double pLimit = state.stepSize.value(cstep::aborter);
        double oLimit = overstepLimit(state.stepping);
        auto checkIntersection = [&](const Intersection& intersection) -> bool {
          double cLimit = intersection.pathLength;
          bool accept = (cLimit > oLimit and cLimit*cLimit < pLimit*pLimit); 
          if (accept){
            double distance = state.navDir*sIntersection.intersection.pathLength;
            // Update the step size 
            state.stepSize.update(distance, cstep::actor, true);
          }
          return true;
        };
        // If either of the two intersections are viable return reachable
        if (checkIntersection(sIntersection.intersection)
            or (sIntersection.alternative 
                  and checkIntersection(sIntersection.alternative))){
              return Intersection::Status::reachable;
        }                                                                                    
      }                                            
      return Intersection::Status::unreachable;    
    }
  }
}