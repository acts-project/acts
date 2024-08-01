// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"

namespace ActsExamples {

///@brief Propagator wrapper
///
/// This class encapsulates a propagator instance and exposes
/// propagation functions for neutral and charged particles.
/// @note This interface exists to decouple the concrete propagators from
///       the @c PropagationAlgorithm
class PropagatorInterface {
 public:
  virtual ~PropagatorInterface() = default;

  ///@brief  Execute a propagation for charged particle parameters
  ///
  ///@param context The algorithm context
  ///@param cfg  The propagation algorithm configuration
  ///@param logger A logger wrapper instance
  ///@param startParameters The start parameters
  ///@return PropagationOutput
  virtual PropagationOutput execute(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const = 0;
};

///@brief Concrete instance of a propagator
///
/// This class implements a wrapped version of a concrete propagator.
/// It can be instantiated to comply with @c PropagatorInterface.
/// Setting up a @c PropagationAlgorithm with this looks like:
/// ```cpp
/// config.propagatorImpl =
///     std::make_shared<ActsExamples::ConcretePropagator<Propagator>>(
///         std::move(propagator));
/// ```
///@tparam propagator_t The concrete propagator to instantiate with
template <typename propagator_t>
class ConcretePropagator : public PropagatorInterface {
 public:
  ConcretePropagator(propagator_t propagator)
      : m_propagator{std::move(propagator)} {}

  PropagationOutput execute(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger,
      const Acts::BoundTrackParameters& startParameters) const override {
    return executeTest(context, cfg, logger, startParameters);
  }

 private:
  /// Templated execute test method for
  /// charged and neutral particles
  /// @param [in] context is the contextual data of this event
  /// @param [in] startParameters the start parameters
  /// @param [in] pathLength the maximal path length to go
  template <typename parameters_t>
  PropagationOutput executeTest(
      const AlgorithmContext& context, const PropagationAlgorithm::Config& cfg,
      const Acts::Logger& logger, const parameters_t& startParameters,
      double pathLength = std::numeric_limits<double>::max()) const {
    ACTS_DEBUG("Test propagation/extrapolation starts");

    PropagationOutput pOutput;

    // This is the outside in mode
    if (cfg.mode == 0) {
      // The step length logger for testing & end of world aborter
      using MaterialInteractor = Acts::MaterialInteractor;
      using SteppingLogger = Acts::detail::SteppingLogger;
      using EndOfWorld = Acts::EndOfWorldReached;

      // Action list and abort list
      using ActionList = Acts::ActionList<SteppingLogger, MaterialInteractor>;
      using AbortList = Acts::AbortList<EndOfWorld>;
      using PropagatorOptions =
          Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

      PropagatorOptions options(context.geoContext, context.magFieldContext);
      options.pathLimit = pathLength;

      // Activate loop protection at some pt value
      options.loopProtection =
          startParameters.transverseMomentum() < cfg.ptLoopers;

      // Switch the material interaction on/off & eventually into logging mode
      auto& mInteractor = options.actionList.get<MaterialInteractor>();
      mInteractor.multipleScattering = cfg.multipleScattering;
      mInteractor.energyLoss = cfg.energyLoss;
      mInteractor.recordInteractions = cfg.recordMaterialInteractions;

      // Switch the logger to sterile, e.g. for timing checks
      auto& sLogger = options.actionList.get<SteppingLogger>();
      sLogger.sterile = cfg.sterileLogger;
      // Set a maximum step size
      options.maxStepSize = cfg.maxStepSize;

      // Propagate using the propagator
      auto result = m_propagator.propagate(startParameters, options);
      if (result.ok()) {
        const auto& resultValue = result.value();
        auto steppingResults =
            resultValue.template get<SteppingLogger::result_type>();

        // Set the stepping result
        pOutput.first = std::move(steppingResults.steps);
        // Also set the material recording result - if configured
        if (cfg.recordMaterialInteractions) {
          auto materialResult =
              resultValue.template get<MaterialInteractor::result_type>();
          pOutput.second = std::move(materialResult);
        }
      }
    }
    return pOutput;
  }

 private:
  propagator_t m_propagator;
};

}  // namespace ActsExamples
