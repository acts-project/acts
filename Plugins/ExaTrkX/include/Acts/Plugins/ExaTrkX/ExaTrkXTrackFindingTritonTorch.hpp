#pragma once
#include "ExaTrkXTriton.hpp"
#include "ExaTrkXTrackFindingBase.hpp"

#include <string>
#include <vector>
#include <memory>

class ExaTrkXTrackFindingTritonTorch : public ExaTrkXTrackFindingBase
{
public:
    struct Config{
        std::string embedModelName;
        std::string filterModelName;
        std::string gnnModelName;
        std::string url;
        bool verbose = false;

        // hyperparameters in the pipeline.
        int64_t spacepointFeatures = 3;
        int embeddingDim = 8;
        float rVal = 1.6;
        int knnVal = 500;
        float filterCut = 0.21;
    };

    ExaTrkXTrackFindingTritonTorch(const Config& config);
    virtual ~ExaTrkXTrackFindingTritonTorch() {}

    void getTracks(
        std::vector<float>& inputValues,
        std::vector<int>& spacepointIDs,
        std::vector<std::vector<int> >& trackCandidates,
        ExaTrkXTime& timeInfo) const final;

    const Config& config() const { return m_cfg; }
    
private:
    Config m_cfg;
    std::unique_ptr<ExaTrkXTriton> e_client_;
    std::unique_ptr<ExaTrkXTriton> f_client_;
    std::unique_ptr<ExaTrkXTriton> g_client_;
};
