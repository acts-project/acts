#pragma once
#include "ExaTrkXTriton.hpp"
#include "ExaTrkXTrackFindingBase.hpp"

#include <string>
#include <vector>
#include <memory>

class ExaTrkXTrackFindingTritonPython : public ExaTrkXTrackFindingBase
{
public:
    struct Config{
        std::string modelDir;
        std::string buildingModelName;
        std::string labelingModelName;
        std::string url;
        bool verbose = false;

        // hyperparameters in the pipeline.
        int64_t spacepointFeatures = 3;
        int embeddingDim = 8;
        float rVal = 1.6;
        int knnVal = 500;
        float filterCut = 0.21;
    };

    ExaTrkXTrackFindingTritonPython(const Config& config);
    virtual ~ExaTrkXTrackFindingTritonPython() {}

    void getTracks(
        std::vector<float>& inputValues,
        std::vector<int>& spacepointIDs,
        std::vector<std::vector<int> >& trackCandidates,
        ExaTrkXTime& timeInfo) const final;

    const Config& config() const { return m_cfg; }
    
private:
    Config m_cfg;
    std::unique_ptr<ExaTrkXTriton> b_client_; // building
    std::unique_ptr<ExaTrkXTriton> l_client_; // labeling
    mutable torch::jit::script::Module e_model;
    mutable torch::jit::script::Module f_model;
    mutable torch::jit::script::Module g_model;
};
