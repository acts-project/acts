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
#include "ActsExamples/Traccc/DetrayStore.hpp"

#include <detray/navigation/navigator.hpp>
#include <detray/propagator/actor_chain.hpp>
#include <detray/propagator/propagation_config.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/test/utils/inspectors.hpp>
#include <detray/test/validation/material_validation_utils.hpp>

namespace ActsExamples {

/// Define the algebra type
using DetrayAlgebraType =
    typename ActsPlugins::DetrayHostDetector::algebra_type;

/// Type that holds the intersection information
using DetrayIntersection = detray::intersection2D<
    typename ActsPlugins::DetrayHostDetector::surface_type, DetrayAlgebraType,
    false>;

using DetrayMaterialTracer =
    detray::material_validator::material_tracer<double, vecmem::vector>;

/// Inspector that records all encountered surfaces
using DetrayObjectTracer =
    detray::navigation::object_tracer<DetrayIntersection, detray::dvector,
                                      detray::navigation::status::e_on_module,
                                      detray::navigation::status::e_on_portal>;

template <typename stepper_t, typename detray_store_t, typename field_t = bool>
class DetrayPropagator : public PropagatorInterface {
 public:
  /// Configuration struct
  struct Config {
    /// The detray store
    std::shared_ptr<const detray_store_t> detrayStore = nullptr;
    /// Switch to sterile
    bool sterile = false;
    /// The field (if any)
    field_t field = field_t();
  };

  /// Create a DetrayPropagator
  ///
  /// @param cfg configuration struct
  /// @param logger The logger instance
  explicit DetrayPropagator(const Config& cfg,
                            std::unique_ptr<const Acts::Logger> logger =
                                Acts::getDefaultLogger("DetrayPropagator",
                                                       Acts::Logging::INFO))
      : PropagatorInterface(), m_cfg(cfg), m_logger(std::move(logger)) {}

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

    // Now follow that ray with the same track and check, if we find
    // the same volumes and distances along the way
    detray::free_track_parameters<DetrayAlgebraType> track(
        {position.x(), position.y(), position.z()}, 0.f,
        {direction.x(), direction.y(), direction.z()},
        startParameters.charge());

    // Return material
    RecordedMaterial recordedMaterial;
    PropagationSummary summary(startParameters);

    if (!m_cfg.sterile) {
      /// Aggregation of multiple inspectors
      using DetrayInspector = detray::aggregate_inspector<DetrayObjectTracer>;

      // Navigation with inspection
      using DetrayNavigator =
          detray::navigator<ActsPlugins::DetrayHostDetector,
                            detray::navigation::default_cache_size,
                            DetrayInspector>;

      // Propagator with empty actor chain (for the moment)
      using Propagator =
          detray::propagator<stepper_t, DetrayNavigator,
                             detray::actor_chain<DetrayMaterialTracer>>;

      using DetrayContext = typename Propagator::state::context_type;
      DetrayContext dCtx{};

      // With and without field dispatch at compile time
      if constexpr (std::is_same<field_t, bool>::value == true) {
        typename Propagator::state propagation(
            track, m_cfg.detrayStore->detector, dCtx);
        propagateAndRecord<Propagator>(propagation, summary, recordedMaterial);
      } else {
        typename Propagator::state propagation(
            track, m_cfg.field, m_cfg.detrayStore->detector, dCtx);
        propagateAndRecord<Propagator>(propagation, summary, recordedMaterial);
      }

    } else {
      // Navigation with inspection
      using DetrayNavigator =
          detray::navigator<ActsPlugins::DetrayHostDetector,
                            detray::navigation::default_cache_size>;

      // Propagator with empty actor chain (for the moment)
      using Propagator =
          detray::propagator<stepper_t, DetrayNavigator, detray::actor_chain<>>;

      using DetrayContext = typename Propagator::state::context_type;
      DetrayContext dCtx{};

      using DetrayConfig = detray::propagation::config;
      DetrayConfig dCfg{};
      Propagator propagator(dCfg);

      // With and without field dispatch at compile time
      if constexpr (std::is_same<field_t, bool>::value == true) {
        typename Propagator::state propagation(
            track, m_cfg.detrayStore->detector, dCtx);
        // Run the detray propagation
        propagator.propagate(propagation);
      } else {
        typename Propagator::state propagation(
            track, m_cfg.field, m_cfg.detrayStore->detector, dCtx);
        // Run the detray propagation
        propagator.propagate(propagation);
      }
    }

    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
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
    auto& inspector = propagation._navigation.inspector();
    auto& objectTracer = inspector.template get<DetrayObjectTracer>();

    // Translate the objects into the steps
    for (const auto& object : objectTracer.object_trace) {
      // Get the position of the object
      const auto& dposition = object.pos;
      const auto& sfDesription = object.intersection.sf_desc;
      const auto sf =
          detray::tracking_surface{m_cfg.detrayStore->detector, sfDesription};
      Acts::GeometryIdentifier geoID{sf.source()};
      // Create a step from the object
      Acts::detail::Step step;
      step.position = Acts::Vector3(dposition[0], dposition[1], dposition[2]);
      step.geoID = geoID;
      step.navDir = object.intersection.direction ? Acts::Direction::Forward()
                                                  : Acts::Direction::Backward();
      summary.steps.emplace_back(step);
    }

    // Retrieve the material information
    const auto& detrayMaterial = materialTracerState.get_material_record();
    recordedMaterial.materialInX0 = detrayMaterial.sX0;
    recordedMaterial.materialInL0 = detrayMaterial.sL0;
  }

  /// The detray detector store and memory resource
  Config m_cfg;

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
