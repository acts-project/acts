#pragma once
#include "ExaTrkXTriton.hpp"
#include "ExaTrkXTrackFindingBase.hpp"


#include <string>
#include <vector>
#include <memory>

class ExaTrkXTrackFindingTritonOne : public ExaTrkXTrackFindingBase
{
public:
    struct Config{
        std::string embedModelName;
        std::string buildingModelName;
        std::string filterModelName;
        std::string gnnModelName;

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

    ExaTrkXTrackFindingTritonOne(const Config& config);
    virtual ~ExaTrkXTrackFindingTritonOne() {}

    void getTracks(
        std::vector<float>& inputValues,
        std::vector<int>& spacepointIDs,
        std::vector<std::vector<int> >& trackCandidates,
        ExaTrkXTime& timeInfo) const final;

    const Config& config() const { return m_cfg; }
    
private:
    Config m_cfg;
    std::unique_ptr<ExaTrkXTriton> e_client_; // embedding
    std::unique_ptr<ExaTrkXTriton> f_client_; // filtering
    std::unique_ptr<ExaTrkXTriton> g_client_; // gnn
    std::unique_ptr<ExaTrkXTriton> b_client_; // building
    std::unique_ptr<ExaTrkXTriton> l_client_; // labeling
};
