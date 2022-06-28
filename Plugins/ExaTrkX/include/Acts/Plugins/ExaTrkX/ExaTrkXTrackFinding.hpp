#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"

#include <memory>
#include <string>
#include <vector>

#include <torch/script.h>
#include <torch/torch.h>
using namespace torch::indexing;

namespace Acts {

class ExaTrkXTrackFinding : public ExaTrkXTrackFindingBase {
 public:
  struct Config {
    std::string modelDir;
    bool verbose = false;

    // hyperparameters in the pipeline.
    int64_t spacepointFeatures = 3;
    int embeddingDim = 8;
    float rVal = 1.6;
    int knnVal = 500;
    float filterCut = 0.21;
    int n_chunks = 5;
    float edgeCut = 0.5;
  };

  ExaTrkXTrackFinding(const Config& config);
  virtual ~ExaTrkXTrackFinding() {}

  void getTracks(std::vector<float>& inputValues,
                 std::vector<int>& spacepointIDs,
                 std::vector<std::vector<int> >& trackCandidates,
                 ExaTrkXTime& timeInfo) const final;

  const Config& config() const { return m_cfg; }

 private:
  void initTrainedModels();

 private:
  Config m_cfg;
  torch::jit::script::Module m_embeddingModel;
  torch::jit::script::Module m_filterModel;
  torch::jit::script::Module m_gnnModel;
};

}  // namespace Acts
