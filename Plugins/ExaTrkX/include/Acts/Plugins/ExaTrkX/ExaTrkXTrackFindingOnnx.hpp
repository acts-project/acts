// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"

#include <memory>

namespace Ort {
class Env;
class Session;
class Value;
}  // namespace Ort

namespace Acts {

/// @brief Implementation of the Exa.TrkX track finding algorithm based on ONNX.
/// Uses cugraph as graph library.
///
class ExaTrkXTrackFindingOnnx final : public ExaTrkXTrackFindingBase {
 public:
  /// Configuration struct for the track finding.
  struct Config {
    // input model directory
    std::string modelDir;

    // hyperparameters in the pipeline.
    int64_t spacepointFeatures = 3;
    int embeddingDim = 8;
    float rVal = 1.6;
    int knnVal = 500;
    float filterCut = 0.21;
  };

  /// Constructor of the track finding module
  ///
  /// @param config is the config struct to configure the module
  ExaTrkXTrackFindingOnnx(const Config& config);

  /// Destructor
  ~ExaTrkXTrackFindingOnnx();

  /// Run the inference
  ///
  /// @param inputValues Spacepoint data as a flattened NxD array, where D is
  /// the dimensionality of a spacepoint (usually 3, but additional information
  /// like cell information can be provided).
  /// @param spacepointIDs The corresponding spacepoint IDs
  /// @param trackCandidates This vector is filled with the tracks as vectors
  /// of spacepoint IDs
  /// @param logger If provided, logging is enabled
  /// @param recordTiming If enabled, returns a @ref ExaTrkXTime object with
  /// measured timings
  /// @note The input values are not const, because the ONNX API
  /// takes only non-const pointers.
  std::optional<ExaTrkXTime> getTracks(
      std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
      std::vector<std::vector<int> >& trackCandidates,
      LoggerWrapper logger = getDummyLogger(),
      bool recordTiming = false) const override;

  /// Return the configuration object of the track finding module
  const Config& config() const { return m_cfg; }

 private:
  void runSessionWithIoBinding(Ort::Session& sess,
                               std::vector<const char*>& inputNames,
                               std::vector<Ort::Value>& inputData,
                               std::vector<const char*>& outputNames,
                               std::vector<Ort::Value>& outputData) const;

  void buildEdges(std::vector<float>& embedFeatures,
                  std::vector<int64_t>& edgeList, int64_t numSpacepoints) const;

 private:
  Config m_cfg;
  std::unique_ptr<Ort::Env> m_env;
  std::unique_ptr<Ort::Session> m_embeddingSession;
  std::unique_ptr<Ort::Session> m_filterSession;
  std::unique_ptr<Ort::Session> m_gnnSession;
};

}  // namespace Acts
