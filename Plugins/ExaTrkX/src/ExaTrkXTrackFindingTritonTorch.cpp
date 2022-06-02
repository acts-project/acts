#include "ExaTrkXTrackFindingTritonTorch.hpp"
#include "ExaTrkXUtils.hpp"

// #include "mmio_read.h"


#include "grpc_client.h"
#include "grpc_service.pb.h"
namespace tc = triton::client;


ExaTrkXTrackFindingTritonTorch::ExaTrkXTrackFindingTritonTorch(
    const ExaTrkXTrackFindingTritonTorch::Config& config):
    ExaTrkXTrackFindingBase("ExaTrkXTrackFindingTritonTorch", config.verbose), m_cfg(config)
{
    bool verbose = false;
    uint32_t client_timeout = 0;
    std::string model_version = "";
    e_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.embedModelName, m_cfg.url, model_version, client_timeout, verbose);
    f_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.filterModelName, m_cfg.url, model_version, client_timeout, verbose);
    g_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.gnnModelName, m_cfg.url, model_version, client_timeout, verbose);
}


// The main function that runs the Exa.TrkX pipeline
// Be care of sharpe corners.
void ExaTrkXTrackFindingTritonTorch::getTracks(
    std::vector<float>& inputValues,
    std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates,
    ExaTrkXTime& timeInfo) const{

    ExaTrkXTimer tot_timer;
    ExaTrkXTimer timer;
    tot_timer.start();
    // ************
    // Embedding
    // ************
    timer.start();
    int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
    std::vector<int64_t> embedInputShape{numSpacepoints, m_cfg.spacepointFeatures};
    // if (m_verbose) {
    //     std::cout << "My input shape is: " << inputValues.size() << std::endl;
    //     std::cout << "My embedding shape is: " << embedInputShape[0] << " " << embedInputShape[1] << std::endl;
    // }


    e_client_->ClearInput();
    e_client_->AddInput<float>("INPUT__0", embedInputShape, inputValues);
    std::vector<float> eOutputData;
    std::vector<int64_t> embedOutputShape{numSpacepoints, m_cfg.embeddingDim};
    e_client_->GetOutput<float>("OUTPUT__0", eOutputData, embedOutputShape);

    timeInfo.embedding = timer.stopAndGetElapsedTime();

    // std::cout <<"Embedding space of the first SP: ";
    // std::copy(eOutputData.begin(), eOutputData.begin() + m_cfg.embeddingDim,
    //           std::ostream_iterator<float>(std::cout, " "));

    // ************
    // Building Edges
    // ************
    timer.start();
    std::vector<int64_t> edgeList;
    buildEdges(
      eOutputData, edgeList, 
      numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);
    int64_t numEdges = edgeList.size() / 2;
    // std::cout << "Built " << numEdges<< " edges." << std::endl;
    // std::cout << edgeList.slice(1, 0, 5) << std::endl;

    timeInfo.building = timer.stopAndGetElapsedTime();


    // ************
    // Filtering
    // ************
    // std::cout << "Get scores for " << numEdges<< " edges." << std::endl;

    timer.start();
    f_client_->ClearInput();
    /// <TODO: reuse the embedding inputs?>
    f_client_->AddInput<float>("INPUT__0", embedInputShape, inputValues);
    std::vector<int64_t> fEdgeShape{2, numEdges};
    f_client_->AddInput<int64_t>("INPUT__1", fEdgeShape, edgeList);

    std::vector<float> fOutputData;
    std::vector<int64_t> fOutputShape{numEdges, 1};
    f_client_->GetOutput<float>("OUTPUT__0", fOutputData, fOutputShape);

    // However, I have to convert those numbers to a score by applying sigmoid!
    // Use torch::tensor
    torch::Tensor edgeListCTen = torch::tensor(edgeList, {torch::kInt64});
    edgeListCTen = edgeListCTen.reshape({2, numEdges});

    torch::Tensor fOutputCTen = torch::tensor(fOutputData, {torch::kFloat32});
    fOutputCTen = fOutputCTen.sigmoid();

    torch::Tensor filterMask = fOutputCTen > m_cfg.filterCut;
    torch::Tensor edgesAfterFCTen = edgeListCTen.index({Slice(), filterMask});

    std::vector<int64_t> edgesAfterFiltering;
    std::copy(
        edgesAfterFCTen.data_ptr<int64_t>(),
        edgesAfterFCTen.data_ptr<int64_t>() + edgesAfterFCTen.numel(),
        std::back_inserter(edgesAfterFiltering));

    int64_t numEdgesAfterF = edgesAfterFiltering.size() / 2;
    // std::cout << "After filtering: " << numEdgesAfterF << " edges." << std::endl;

    timeInfo.filtering = timer.stopAndGetElapsedTime();

    // ************
    // GNN
    // ************
    timer.start();

    g_client_->ClearInput();
    g_client_->AddInput<float>("INPUT__0", embedInputShape, inputValues);
    std::vector<int64_t> gEdgeShape{2, numEdgesAfterF};
    g_client_->AddInput<int64_t>("INPUT__1", gEdgeShape, edgesAfterFiltering);

    std::vector<float> gOutputData;
    std::vector<int64_t> gOutputShape{numEdgesAfterF, 1};
    g_client_->GetOutput<float>("OUTPUT__0", gOutputData, gOutputShape);

    torch::Tensor gOutputCTen = torch::tensor(gOutputData, {torch::kFloat32});
    gOutputCTen = gOutputCTen.sigmoid();

    timeInfo.gnn = timer.stopAndGetElapsedTime();
    // ************
    // Track Labeling with cugraph::connected_components
    // ************
    timer.start();

    std::vector<int32_t> rowIndices;
    std::vector<int32_t> colIndices;
    std::vector<float> edgeWeights;
    std::vector<int32_t> trackLabels(numSpacepoints);
    std::copy(
        edgesAfterFiltering.begin(),
        edgesAfterFiltering.begin()+numEdgesAfterF,
        std::back_insert_iterator(rowIndices));
    std::copy(
        edgesAfterFiltering.begin()+numEdgesAfterF,
        edgesAfterFiltering.end(),
        std::back_insert_iterator(colIndices));
    std::copy(
        gOutputCTen.data_ptr<float>(),
        gOutputCTen.data_ptr<float>() + numEdgesAfterF,
        std::back_insert_iterator(edgeWeights));

    // weakly_connected_components<int32_t,int32_t,float>(
    //     rowIndices, colIndices, edgeWeights, trackLabels);
    weaklyConnectedComponents<int32_t,int32_t,float>(
        numSpacepoints, 
        rowIndices, colIndices, edgeWeights, trackLabels);

    int idx = 0;
    // std::cout << "size of components: " << trackLabels.size() << std::endl;
    if (trackLabels.size() == 0)  return;


    trackCandidates.clear();

    int existTrkIdx = 0;
    // map labeling from MCC to customized track id.
    std::map<int32_t, int32_t> trackLableToIds;

    for(int32_t idx=0; idx < numSpacepoints; ++idx) {
        int32_t trackLabel = trackLabels[idx];
        int spacepointID = spacepointIDs[idx];

        int trkId;
        if(trackLableToIds.find(trackLabel) != trackLableToIds.end()) {
            trkId = trackLableToIds[trackLabel];
            trackCandidates[trkId].push_back(spacepointID);
        } else {
            // a new track, assign the track id
            // and create a vector
            trkId = existTrkIdx;
            trackCandidates.push_back(std::vector<int>{spacepointID});
            trackLableToIds[trackLabel] = trkId;
            existTrkIdx++;
        }
    }
    timeInfo.labeling = timer.stopAndGetElapsedTime();
    timeInfo.total = tot_timer.stopAndGetElapsedTime();
}
