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
  /// @param logger The logger instance
  explicit DetraySterilePropagator(
      std::shared_ptr<detector_t> dGeometry,
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

    // Now follow that ray with the same track and check, if we find
    // the same volumes and distances along the way
    using DetrayAlgebraType = typename detector_t::algebra_type;

    detray::free_track_parameters<DetrayAlgebraType> track(
        {static_cast<float>(position.x()), static_cast<float>(position.y()),
         static_cast<float>(position.z())},
        0.f,
        {static_cast<float>(direction.x()), static_cast<float>(direction.y()),
         static_cast<float>(direction.z())},
        static_cast<float>(startParameters.charge()));

    // Return material
    RecordedMaterial recordedMaterial;
    PropagationSummary summary(startParameters);

    using navigator_t = detray::caching_navigator<detector_t>;
    using propagator_t =
        detray::propagator<stepper_t, navigator_t, detray::actor_chain<>>;

    using DetrayConfig = detray::propagation::config;
    DetrayConfig dCfg{};
    propagator_t propagator(dCfg);

    // With and without field dispatch at compile time
    if constexpr (std::is_same<field_t, bool>::value == true) {
      typename propagator_t::state propagation(track, *m_detrayGeometry,
                                               dCfg.context);
      // Run the detray propagation
      propagator.propagate(propagation);
    } else {
      typename propagator_t::state propagation(track, m_field, m_detrayGeometry,
                                               dCfg.context);
      // Run the detray propagation
      propagator.propagate(propagation);
    }
    // Return empty results
    return std::pair{std::move(summary), std::move(recordedMaterial)};
  }

 private:
  /// The detray detector store and memory resource
  std::shared_ptr<detector_t> m_detrayGeometry;

  field_t m_field = field_t();

  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger = nullptr;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
