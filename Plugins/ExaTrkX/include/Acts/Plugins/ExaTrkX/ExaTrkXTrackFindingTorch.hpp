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

namespace torch::jit {
class Module;
}

namespace Acts {

/// @brief Class implementing the Exa.TrkX track finding algorithm based on
/// libtorch. Uses Boost.Graph for as graph library
///
class ExaTrkXTrackFindingTorch final : public ExaTrkXTrackFindingBase {
 public:
  /// Configuration struct for the track finding
  struct Config {
    std::string modelDir;

    // hyperparameters in the pipeline.
    int64_t spacepointFeatures = 3;
    int embeddingDim = 8;
    float rVal = 1.6;
    int knnVal = 500;
    float filterCut = 0.21;
    int n_chunks = 5;
    float edgeCut = 0.5;
  };

  /// Constructor of the track finding module
  ///
  /// @param config is the config struct to configure the module
  ExaTrkXTrackFindingTorch(const Config& config);

  /// Destructor
  ~ExaTrkXTrackFindingTorch();

  /// Run the inference
  ///
  /// @param inputValues Spacepoint data as a flattened NxD array, where D is
  /// the dimensionality of a spacepoint (usually 3, but additional information
  /// like cell information can be provided).
  /// @param spacepointIDs The corresponding spacepoint IDs
  /// @param trackCandidates This vector is filled with the tracks as vectors of spacepoint IDs
  /// @param logger If provided, logging is enabled
  /// @param recordTiming If enabled, returns a ExaTrkXTime object with measured timings
  /// @note The input values are not const, because the ONNX API
  /// takes only non-const pointers.
  std::optional<ExaTrkXTime> getTracks(
      std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
      std::vector<std::vector<int> >& trackCandidates,
      LoggerWrapper logger = getDummyLogger(),
      bool recordTiming = false) const override;

  /// Access the config struct
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::unique_ptr<torch::jit::Module> m_embeddingModel;
  std::unique_ptr<torch::jit::Module> m_filterModel;
  std::unique_ptr<torch::jit::Module> m_gnnModel;
};

}  // namespace Acts
