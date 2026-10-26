// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/interface/IAssignmentFinder.hpp"
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/DefaultActor.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <map>
#include <utility>
#include <vector>

namespace Acts {

/// @brief selector for finding surface
///
/// This is to be used with a SurfaceCollector<>
struct MaterialSurfaceIdentifier {
  /// check if the surface has material
  /// @param sf Surface to check for material
  /// @return True if the surface has material assigned to it
  bool operator()(const Surface& sf) const { return sf.hasMaterial(); }
};

/// An Interaction volume collector with unique counting
struct InteractionVolumeCollector {
  /// Simple state struct that holds the collected volumes
  ///
  /// @note the map is to avoid double counting as this is
  /// called in an action list
  struct this_result {
    /// Map of collected volume assignments by geometry identifier
    std::map<GeometryIdentifier, IAssignmentFinder::VolumeAssignment> collected;
  };

  /// Type alias for the volume collection state
  using state_type = this_result;

  /// Collector action for the ActionList of the Propagator
  /// It checks if the propagator state has a current volume,
  /// in which case the action is performed:
  /// - it records the volume given the configuration
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper used for the propagation
  /// @tparam navigator_t Type of the navigator used for the propagation
  ///
  /// @param [in,out] state is the mutable stepper state object
  /// @param [in] stepper The stepper in use
  /// @param [in] navigator The navigator in use
  /// @param [in,out] actorState is the mutable state of the actor
  /// @return Result object indicating success or failure
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, state_type& actorState,
                   const Logger& /*logger*/) const {
    // Retrieve the current volume
    auto currentVolume = navigator.currentVolume(state.navigation);

    // The current volume has been assigned by the navigator
    if (currentVolume != nullptr) {
      auto collIt = actorState.collected.find(currentVolume->geometryId());
      // Check if the volume has been collected and if it has material
      if (collIt == actorState.collected.end() &&
          currentVolume->hasMaterial()) {
        Vector3 entryPosition = stepper.position(state.stepping);
        Vector3 exitPosition = entryPosition;
        IAssignmentFinder::VolumeAssignment vAssignment{
            InteractionVolume(currentVolume), entryPosition, exitPosition};
        actorState.collected.emplace(currentVolume->geometryId(), vAssignment);
      } else if (collIt != actorState.collected.end()) {
        // Update the exit position
        (collIt->second).exit = stepper.position(state.stepping);
      }
    }
    return Result<void>::success();
  }
};

/// @class PropagatorMaterialAssigner
///
/// A Propagator based material assigner that uses the navigation and
/// transport of the propagator to assign the material to the surface
/// or the volume.
///
/// @note eventual navigation problems would affect he material mapping
///
/// @ingroup material_mapping
template <typename propagator_t>
class PropagatorMaterialAssigner /*final*/ : public IAssignmentFinder {
 public:
  /// @brief  Construct with propagator
  /// @param propagator
  explicit PropagatorMaterialAssigner(propagator_t propagator)
      : m_propagator(std::move(propagator)) {}

  /// @brief Method for generating assignment candidates for the
  /// material interaction assignment to surfaces or volumes
  ///
  /// @param gctx is the geometry context
  /// @param mctx is the magnetic field context
  /// @param position is the position of the initial ray
  /// @param direction is the direction of initial ray
  ///
  /// @return a vector of Surface Assignments and Volume Assignments
  std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
            std::vector<IAssignmentFinder::VolumeAssignment>>
  assignmentCandidates(const GeometryContext& gctx,
                       const MagneticFieldContext& mctx,
                       const Vector3& position,
                       const Vector3& direction) const final {
    // Return container
    std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
              std::vector<IAssignmentFinder::VolumeAssignment>>
        candidates;

    using VectorHelpers::makeVector4;
    // Neutral curvilinear parameters
    BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
        makeVector4(position, 0), direction, 1, std::nullopt,
        ParticleHypothesis::geantino());

    // Collect the material surfaces and volumes up to the end of the world
    using MaterialSurfaceCollector =
        SurfaceCollector<MaterialSurfaceIdentifier>;
    using ActorList = ActorList<DefaultActor, MaterialSurfaceCollector,
                                InteractionVolumeCollector>;
    using PropagatorOptions =
        typename propagator_t::template Options<ActorList>;

    PropagatorOptions options(gctx, mctx);
    options.stopAtEndOfWorld = true;

    auto state = m_propagator.makeState(options);
    m_propagator.propagate(state, start, nullptr, false).value();

    // The surface collection results
    const auto& scResult =
        state.actorState.template get<MaterialSurfaceCollector::state_type>();

    for (const auto& mSurface : scResult.collected) {
      candidates.first.push_back(IAssignmentFinder::SurfaceAssignment{
          mSurface.surface, mSurface.position, direction});
    }

    // The volume collection results
    const auto& vcResult =
        state.actorState.template get<InteractionVolumeCollector::state_type>();
    for (const auto& [geoId, vIntersection] : vcResult.collected) {
      candidates.second.push_back(vIntersection);
    }

    // Return the candidates
    return candidates;
  }

 private:
  propagator_t m_propagator;
};

}  // namespace Acts
