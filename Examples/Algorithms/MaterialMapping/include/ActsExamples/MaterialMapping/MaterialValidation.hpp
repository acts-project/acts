// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialValidator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <numbers>

namespace ActsExamples {

/// @class MaterialValidation
///
class MaterialValidation : public IAlgorithm {
 public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config {
    /// Input track parameters
    std::string inputTrackParameters = "InputTrackParameters";

    /// Output collection name
    std::string outputMaterialTracks = "ValidationMaterialTracks";

    // The validator
    std::shared_ptr<Acts::MaterialValidator> materialValidator = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  explicit MaterialValidation(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Destructor
  /// - it also writes out the file
  ~MaterialValidation() override = default;

  /// Framework execute method
  ///
  /// @param context The algorithm context for event consistency
  ProcessCode execute(const AlgorithmContext& context) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;  //!< internal config object

  ReadDataHandle<TrackParametersContainer> m_inputTrackParameters{
      this, "InputTrackParameters"};

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "ValidationMaterialTracks"};
};

}  // namespace ActsExamples
