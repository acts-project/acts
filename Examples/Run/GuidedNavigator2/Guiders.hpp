#pragma once

#include "Acts/Surfaces/Surface.hpp"

#include <vector>



namespace Acts {

struct DirectGuider {
  template <typename propagator_state_t, typename stepper_t>
  void updateUnreachable(propagator_state_t &state, const stepper_t &) const {
    ++state.navigation.navSurfaceIter;
  }

  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t &state, const stepper_t &) const {
    ++state.navigation.navSurfaceIter;
  }
};

struct TrialAndErrorGuider {
  std::vector<const Surface *> m_possibleSurfaces;

  template <typename propagator_state_t, typename stepper_t>
  void updateUnreachable(propagator_state_t &state, const stepper_t &) const {
    ++state.navigation.navSurfaceIter;
  }

  template <typename propagator_state_t, typename stepper_t>
  void updateOnSurface(propagator_state_t &state, const stepper_t &) const {
    const auto &logger = state.options.logger;

    // Important parameters
    const FreeVector &params = state.stepping.pars;
    const GeometryContext &gctx = state.geoContext;
    //     const double olimit = stepper.overstepLimit(state.stepping);

    ACTS_VERBOSE("updateOnSurface at "
                 << params.segment<3>(eFreePos0).transpose() << " with dir "
                 << params.segment<3>(eFreeDir0).transpose());

    // Make intersection objects
    std::vector<SurfaceIntersection> sfis(m_possibleSurfaces.size());
    std::transform(m_possibleSurfaces.begin(), m_possibleSurfaces.end(),
                   sfis.begin(), [&](auto s) {
                     return s->intersect(
                         gctx, params.segment<3>(eFreePos0),
                         state.stepping.navDir * params.segment<3>(eFreeDir0),
                         true);
                   });

    ACTS_VERBOSE("have " << sfis.size() << " candidate intersections");

    // TODO overstepLimit etc. missing here
    sfis.erase(
        std::remove_if(sfis.begin(), sfis.end(),
                       [](auto sfi) {
                         return sfi.intersection.status !=
                                    Acts::Intersection3D::Status::reachable ||
                                sfi.intersection.pathLength < 0;
                       }),
        sfis.end());

    ACTS_VERBOSE("after remove " << sfis.size() << " candidates remain");

    for (const auto &s : sfis)
      std::cout << s.intersection.pathLength
                << "\tsurface: " << s.object->center(gctx).transpose() << "\n";

    if (state.stepping.navDir == forward) {
      std::sort(sfis.begin(), sfis.end());
    } else {
      std::sort(sfis.begin(), sfis.end(), std::greater<>());
    }

    // Transform back to vector of surface pointers
    std::vector<const Surface *> finalSurfaces;
    finalSurfaces.reserve(sfis.size());
    for (const auto &si : sfis) {
      finalSurfaces.push_back(si.representation);
    }

    std::cout << "SURFACE CANDIDATES:\n";
    for (const auto s : finalSurfaces)
      std::cout << s->geometryId() << "\t" << s->center(gctx).transpose()
                << "\n";

    state.navigation.navSurfaces = std::move(finalSurfaces);
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
  }
};

}  // namespace Acts
