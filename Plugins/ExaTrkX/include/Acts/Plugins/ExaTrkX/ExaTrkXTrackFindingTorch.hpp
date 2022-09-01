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
#include <string>
#include <vector>

namespace torch::jit {
class Module;
}

namespace Acts {

class ExaTrkXTrackFindingTorch final : public ExaTrkXTrackFindingBase {
 public:
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

  ExaTrkXTrackFindingTorch(const Config& config);
  virtual ~ExaTrkXTrackFindingTorch();

  void getTracks(std::vector<float>& inputValues,
                 std::vector<int>& spacepointIDs,
                 std::vector<std::vector<int> >& trackCandidates,
                 ExaTrkXTime& timeInfo,
                 LoggerWrapper logger = getDummyLogger()) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::unique_ptr<torch::jit::Module> m_embeddingModel;
  std::unique_ptr<torch::jit::Module> m_filterModel;
  std::unique_ptr<torch::jit::Module> m_gnnModel;
};

}  // namespace Acts
