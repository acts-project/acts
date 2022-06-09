#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFindingBase.hpp"

#include <string>
#include <vector>
#include <memory>

#include <torch/torch.h>
#include <torch/script.h>
using namespace torch::indexing;

namespace Acts {

class ExaTrkXTrackFinding : public ExaTrkXTrackFindingBase
{
public:
    struct Config{
        std::string modelDir;
        bool verbose = false;

        // hyperparameters in the pipeline.
        int64_t spacepointFeatures = 3;
        int embeddingDim = 8;
        float rVal = 1.6;
        int knnVal = 500;
        float filterCut = 0.21;
        float edgeCut = 0.5;
    };


    ExaTrkXTrackFinding(const Config& config);
    virtual ~ExaTrkXTrackFinding() {}

    void getTracks(
        std::vector<float>& inputValues,
        std::vector<int>& spacepointIDs,
        std::vector<std::vector<int> >& trackCandidates,
        ExaTrkXTime& timeInfo) const final;

    const Config& config() const { return m_cfg; }

private:
    void initTrainedModels();

private:
    Config m_cfg;
    mutable torch::jit::script::Module e_model;
    mutable torch::jit::script::Module f_model;
    mutable torch::jit::script::Module g_model;
};

}
