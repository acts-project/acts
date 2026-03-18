// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MaterialMapping/MaterialValidation.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

MaterialValidation::MaterialValidation(
    const MaterialValidation::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MaterialValidation", std::move(logger)), m_cfg(cfg) {
  // Prepare the I/O collections
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_outputMaterialTracks.initialize(m_cfg.outputMaterialTracks);
  // Check the configuration - material validator
  if (m_cfg.materialValidator == nullptr) {
    throw std::invalid_argument("Missing material validator.");
  }
}

ProcessCode MaterialValidation::execute(const AlgorithmContext& context) const {
  const auto& inputTrackParameters = m_inputTrackParameters(context);

  ACTS_DEBUG("Validation with " << inputTrackParameters.size()
                                << " input trackparameters");

  std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>
      recordedMaterialTracks;
  recordedMaterialTracks.reserve(inputTrackParameters.size());

  for (const auto& [it, parameters] : Acts::enumerate(inputTrackParameters)) {
    // Record the material
    auto rMaterial = m_cfg.materialValidator->recordMaterial(
        context.geoContext, context.magFieldContext,
        parameters.position(context.geoContext),
        parameters.momentum().normalized());

    recordedMaterialTracks.try_emplace(recordedMaterialTracks.end(), it,
                                       rMaterial);
  }

  // Write the mapped and unmapped material tracks to the output
  m_outputMaterialTracks(context, std::move(recordedMaterialTracks));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
