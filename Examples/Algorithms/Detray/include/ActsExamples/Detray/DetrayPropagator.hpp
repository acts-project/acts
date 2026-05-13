// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/PropagationSummary.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <covfie/core/field.hpp>
#include <detray/navigation/caching_navigator.hpp>
#include <detray/propagator/actor_chain.hpp>
#include <detray/propagator/propagation_config.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/test/common/bfield.hpp>
#include <detray/test/utils/inspectors.hpp>
#include <detray/test/validation/material_validation_utils.hpp>

namespace ActsExamples {

template <typename stepper_t, typename detector_t, typename field_t = bool>
class DetraySterilePropagator : public PropagatorInterface {
 public:
  /// Create a DetrayPropagator
  ///
  /// @param dGeometry configuration struct
  /// @param sterile whether to run the propagator in sterile mode
  /// @param logger The logger instance
  explicit DetraySterilePropagator(
      detector_t&& dGeometry,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "DetraySterilePropagator", Acts::Logging::INFO))
      : PropagatorInterface(),
        m_detrayGeometry(std::move(dGeometry)),
        m_logger(std::move(logger)) {}

  ///@brief  Execute a propagation for charged particle parameters
  ///
  ///@param context The algorithm context
  ///@param cfg  The propagation algorithm configuration
  ///@param logger A logger wrapper instance
  ///@param startParameters The start parameters
  ///@return PropagationOutput
  Acts::Result<PropagationOutput> execute(
      const AlgorithmContext& context,
      [[maybe_unused]] const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const final {
    // Get the geometry context form the algorithm context
    const auto& geoContext = context.geoContext;
    // Get the track information
    const Acts::Vector3 position = startParameters.position(geoContext);
    const Acts::Vector3 direction = startParameters.momentum().normalized();

    ACTS_VERBOSE("Starting propagation at " << position.transpose()
                                            << " with direction "
                                            << direction.transpose());

    // Return material
    RecordedMaterial recordedMaterial;
    PropagationSummary summary(startParameters);

    /**
// Now follow that ray with the same track and check, if we find
// the same volumes and distances along the way
detray::free_track_parameters<DetrayAlgebraType> track(
{position.x(), position.y(), position.z()}, 0.f,
{direction.x(), direction.y(), direction.z()},
startParameters.charge());



if (!m_sterile) {
/// Aggregation of multiple inspectors
using DetrayInspector = detray::aggregate_inspector<DetrayObjectTracer>;

// Navigation with inspection
using DetrayNavigator =
detray::caching_navigator<detector_t,
detray::navigation::default_cache_size,
DetrayInspector, DetrayIntersection>;

// Propagator with empty actor chain (for the moment)
using Propagator =
detray::propagator<stepper_t, DetrayNavigator,
detray::actor_chain<DetrayMaterialTracer>>;

detray::propagation::config prop_cfg{};
auto dCtx = prop_cfg.context;

// With and without field dispatch at compile time
if constexpr (std::is_same<field_t, bool>::value == true) {
typename Propagator::state propagation(track, m_detrayGeometry, dCtx);
propagateAndRecord<Propagator>(propagation, summary, recordedMaterial);
} else {
typename Propagator::state propagation(track, m_cfg.field,
       m_detrayGeometry, dCtx);
propagateAndRecord<Propagator>(propagation, summary, recordedMaterial);
}

} else {
// Navigation with inspection
using DetrayNavigator = detray::caching_navigator<detector_t>;

// Propagator with empty actor chain (for the moment)
using Propagator =
detray::propagator<stepper_t, DetrayNavigator, detray::actor_chain<>>;

detray::propagation::config prop_cfg{};
auto dCtx = prop_cfg.context;

using DetrayConfig = detray::propagation::config;
DetrayConfig dCfg{};
Propagator propagator(dCfg);

// With and without field dispatch at compile time
if constexpr (std::is_same<field_t, bool>::value == true) {
typename Propagator::state propagation(track, m_detrayGeometry, dCtx);
// Run the detray propagation
propagator.propagate(propagation);
} else {
typename Propagator::state propagation(track, m_cfg.field,
       m_detrayGeometry, dCtx);
// Run the detray propagation
propagator.propagate(propagation);
}
}
**/
    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  /**
  /// Propagate and record the steps and material
  /// @tparam propagator_type The type of the propagator
  /// @param propagation The propagation state
  /// @param summary The propagation summary to fill
  /// @param recordedMaterial The recorded material to fill
  template <typename propagator_type>
  void propagateAndRecord(typename propagator_type::state& propagation,
                          PropagationSummary& summary,
                          RecordedMaterial& recordedMaterial) const {
    using DetrayConfig = detray::propagation::config;
    DetrayConfig dCfg{};

    propagator_type propagator(dCfg);

    DetrayMaterialTracer::state materialTracerState{
        *m_cfg.detrayStore->memoryResource};

    auto actorStates = detray::tie(materialTracerState);

    // Run the detray propagation
    propagator.propagate(propagation, actorStates);

    // Retrieve navigation information
    auto& inspector = propagation.navigation().inspector();
    auto& objectTracer = inspector.template get<DetrayObjectTracer>();

    // Translate the objects into the steps
    for (const auto& object : objectTracer.object_trace) {
      // Get the position of the object
      const auto& dposition = object.pos;
      const auto& sfDesription = object.intersection.surface();
      const auto sf = detray::tracking_surface{m_detrayGeometry, sfDesription};
      Acts::GeometryIdentifier geoID(sf.source());
      // Create a step from the object
      Acts::detail::Step step;
      step.position = Acts::Vector3(dposition[0], dposition[1], dposition[2]);
      step.geoID = geoID;
      step.navDir = object.intersection.is_along()
                        ? Acts::Direction::Forward()
                        : Acts::Direction::Backward();
      summary.steps.emplace_back(step);
    }

    // Retrieve the material information
    const auto& detrayMaterial = materialTracerState.get_material_record();
    recordedMaterial.materialInX0 = detrayMaterial.sX0;
    recordedMaterial.materialInL0 = detrayMaterial.sL0;
  }
  */

  /// The detray detector store and memory resource
  detector_t m_detrayGeometry;

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
