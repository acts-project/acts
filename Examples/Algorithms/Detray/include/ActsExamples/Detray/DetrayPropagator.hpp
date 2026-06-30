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

#include <array>

#include <covfie/core/field.hpp>
#include <detray/navigation/caching_navigator.hpp>
#include <detray/propagator/actor_chain.hpp>
#include <detray/propagator/propagation_config.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/test/common/bfield.hpp>
#include <detray/test/utils/inspectors.hpp>
#include <detray/test/validation/material_validation_utils.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsExamples {

template <typename stepper_t, typename detector_t, typename field_t,
          bool kSTERILE>
class DetrayConcretePropagator : public PropagatorInterface {
 public:
  /// Create a DetrayConcretePropagator - this is for testing purposes only,
  /// it runs without actor chain.
  ///
  /// @param dGeometry the detray geometry
  /// @param memoryResource the memory resource to use for the detray propagation
  /// @param logger The logger instance
  /// @param searchWindow the grid acceleration search window (per axis).
  ///        Defaults to {0, 0}, i.e. detray only looks at the current bin and
  ///        does no neighbor lookup. Increasing this lets a client probe
  ///        whether the surface-grid filling needs to merge neighbor cells:
  ///        if results change with a larger window, the converted grids are
  ///        missing neighbor entries and the filling should be optimized.
  explicit DetrayConcretePropagator(
      std::shared_ptr<detector_t> dGeometry,
      vecmem::memory_resource& memoryResource,
      std::unique_ptr<const Acts::Logger> logger =
          Acts::getDefaultLogger("DetrayPropagator", Acts::Logging::INFO),
      std::array<unsigned int, 2> searchWindow = {0u, 0u})
      : PropagatorInterface(),
        m_detrayGeometry(std::move(dGeometry)),
        m_memoryResource(&memoryResource),
        m_searchWindow(searchWindow),
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

    // Now follow that ray with the same track and check, if we find
    // the same volumes and distances along the way
    using DAlgebra = typename detector_t::algebra_type;

    detray::free_track_parameters<DAlgebra> track(
        {static_cast<float>(position.x()), static_cast<float>(position.y()),
         static_cast<float>(position.z())},
        0.f,
        {static_cast<float>(direction.x()), static_cast<float>(direction.y()),
         static_cast<float>(direction.z())},
        static_cast<float>(startParameters.charge()));

    // Return material
    RecordedMaterial recordedMaterial;
    PropagationSummary summary(startParameters);

    detray::propagation::config dCfg{};

    // Configure the grid acceleration search window. {0, 0} (the default) means
    // detray only inspects the current bin and performs no neighbor lookup.
    // A client can widen this to test whether the converted surface grids are
    // missing neighbor-cell entries (see ctor documentation).
    dCfg.navigation.search_window = {m_searchWindow[0], m_searchWindow[1]};

    // Sterile propagation without actor chain
    if constexpr (kSTERILE) {
      using DNavigator = detray::caching_navigator<detector_t>;
      using DPropagator =
          detray::propagator<stepper_t, DNavigator, detray::actor_chain<>>;

      DPropagator propagator(dCfg);
      // With and without field dispatch at compile time
      if constexpr (std::is_same<field_t, bool>::value == true) {
        typename DPropagator::state propagation(track, *m_detrayGeometry,
                                                dCfg.context);
        // Run the detray propagation
        propagator.propagate(propagation);
      } else {
        typename DPropagator::state propagation(
            track, m_field, *m_detrayGeometry, dCfg.context);
        // Run the detray propagation
        propagator.propagate(propagation);
      }
    } else {
      /// Types for propagation tracing and inspection
      using DIntersection =
          detray::intersection2D<typename detector_t::surface_type, DAlgebra>;
      using DMaterialTracer =
          detray::material_validator::material_tracer<float, vecmem::vector>;
      using DObjectTracer = detray::navigation::object_tracer<
          detector_t, detray::dvector, detray::navigation::status::e_on_object,
          detray::navigation::status::e_on_portal>;
      using DInspector = detray::aggregate_inspector<DObjectTracer>;

      // Navigation with inspection
      using DNavigator =
          detray::caching_navigator<detector_t,
                                    detray::navigation::default_cache_size,
                                    DInspector, DIntersection>;

      // Propagator with empty actor chain (for the moment)
      using DPropagator =
          detray::propagator<stepper_t, DNavigator,
                             detray::actor_chain<DMaterialTracer>>;

      // With and without field dispatch at compile time
      if constexpr (std::is_same<field_t, bool>::value == true) {
        typename DPropagator::state propagation(track, *m_detrayGeometry,
                                                dCfg.context);
        DPropagator propagator(dCfg);

        propagateAndRecord<DMaterialTracer, DObjectTracer>(
            propagation, propagator, summary, recordedMaterial);
      } else {
        typename DPropagator::state propagation(
            track, m_field, *m_detrayGeometry, dCfg.context);
        DPropagator propagator(dCfg);
        propagateAndRecord<DMaterialTracer, DObjectTracer>(
            propagation, propagator, summary, recordedMaterial);
      }
    }

    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  /// Propagate and record the steps and material
  /// @tparam propagator_type The type of the propagator
  /// @param propagation The propagation state
  /// @param propagator The propagator to run
  /// @param summary The propagation summary to fill
  /// @param recordedMaterial The recorded material to fill
  template <typename material_tracer_t, typename object_tracer_t,
            typename propagator_t>
  void propagateAndRecord(typename propagator_t::state& propagation,
                          const propagator_t& propagator,
                          PropagationSummary& summary,
                          RecordedMaterial& recordedMaterial) const {
    typename material_tracer_t::state materialTracerState{*m_memoryResource};

    auto actorStates = detray::tie(materialTracerState);

    // Run the detray propagation
    propagator.propagate(propagation, actorStates);

    // Retrieve navigation information
    auto& inspector = propagation.navigation().inspector();
    auto& objectTracer = inspector.template get<object_tracer_t>();

    // Translate the objects into the steps
    for (const auto& object : objectTracer.object_trace) {
      // Get the position of the object
      const auto& dposition = object.pos;
      const auto& sfDesription = object.intersection.surface();
      const auto sf = detray::tracking_surface{*m_detrayGeometry, sfDesription};
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

  /// The detray detector store and memory resource
  std::shared_ptr<detector_t> m_detrayGeometry = nullptr;
  vecmem::memory_resource* m_memoryResource = nullptr;

  /// Grid acceleration search window (per axis) applied to the navigation
  /// config. {0, 0} disables neighbor lookup.
  std::array<unsigned int, 2> m_searchWindow = {0u, 0u};

  field_t m_field = field_t();

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

template <typename stepper_t, typename detector_t, typename field_t = bool>
using DetraySterilePropagator =
    DetrayConcretePropagator<stepper_t, detector_t, field_t, true>;

template <typename stepper_t, typename detector_t, typename field_t = bool>
using DetrayPropagator =
    DetrayConcretePropagator<stepper_t, detector_t, field_t, false>;

}  // namespace ActsExamples
