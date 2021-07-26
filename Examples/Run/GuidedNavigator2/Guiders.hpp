#pragma once

#include "Acts/Surfaces/Surface.hpp"

#include <vector>

namespace Acts {

struct DirectGuider {
  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t &state, const stepper_t &) const {
    ++state.navigation.navSurfaceIter;
  }
};

struct TrialAndErrorGuider {
  std::vector<const Surface *> possibleSurfaces;

  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t &state,
                       const stepper_t &stepper) const {
    const auto &logger = state.options.logger;

    // Important parameters
    const FreeVector &params = state.stepping.pars;
    const GeometryContext &gctx = state.geoContext;
    const double oLimit = stepper.overstepLimit(state.stepping);
    const double pLimit = state.options.pathLimit;

    ACTS_VERBOSE("updateOnSurface at "
                 << params.segment<3>(eFreePos0).transpose() << " with dir "
                 << params.segment<3>(eFreeDir0).transpose());

    // Make intersection objects
    std::vector<SurfaceIntersection> sfis(possibleSurfaces.size());
    std::transform(possibleSurfaces.begin(), possibleSurfaces.end(),
                   sfis.begin(), [&](auto s) {
                     return s->intersect(
                         gctx, params.segment<3>(eFreePos0),
                         state.stepping.navDir * params.segment<3>(eFreeDir0),
                         true);
                   });

    ACTS_VERBOSE("have " << sfis.size() << " candidate intersections");

    // Filter out intersections which are not valid
    auto isIntersectionValid = [&](const SurfaceIntersection &sfi) {
      const double sfiPath = sfi.intersection.pathLength;
      return sfi.intersection.status == Intersection3D::Status::reachable &&
             sfiPath > oLimit && sfiPath * sfiPath <= pLimit * pLimit;
    };

    sfis.erase(std::remove_if(sfis.begin(), sfis.end(),
                              [&](const auto &sfi) {
                                return not isIntersectionValid(sfi);
                              }),
               sfis.end());

    ACTS_VERBOSE("after remove " << sfis.size() << " candidates remain");

    std::sort(sfis.begin(), sfis.end());
    for(auto &sfi : sfis) {
      sfi.intersection.pathLength *= std::copysign(1., state.stepping.navDir);
    }

    // Transform back to vector of surface pointers
    std::vector<const Surface *> finalSurfaces;
    finalSurfaces.reserve(sfis.size());
    for (const auto &si : sfis) {
      finalSurfaces.push_back(si.representation);
    }

    state.navigation.navSurfaces = std::move(finalSurfaces);
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
  }
};

}  // namespace Acts
