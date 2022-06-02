#include "ExaTrkXTrackFindingTritonOne.hpp"
#include "ExaTrkXUtils.hpp"

#include "grpc_client.h"
#include "grpc_service.pb.h"
namespace tc = triton::client;


ExaTrkXTrackFindingTritonOne::ExaTrkXTrackFindingTritonOne(
    const ExaTrkXTrackFindingTritonOne::Config& config):
    ExaTrkXTrackFindingBase("ExaTrkXTrackFindingTritonOne", config.verbose), m_cfg(config)
{

    bool verbose = false;
    uint32_t client_timeout = 0;
    std::string model_version = "";

    e_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.embedModelName, m_cfg.url, model_version, client_timeout, verbose);
    b_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.buildingModelName, m_cfg.url, model_version, client_timeout, verbose);
    f_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.filterModelName, m_cfg.url, model_version, client_timeout, verbose);
    g_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.gnnModelName, m_cfg.url, model_version, client_timeout, verbose);
    l_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.labelingModelName, m_cfg.url, model_version, client_timeout, verbose);
}



// The main function that runs the Exa.TrkX ExaTrkXTrackFinding pipeline
void ExaTrkXTrackFindingTritonOne::getTracks(
    std::vector<float>& inputValues,
    std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates,
    ExaTrkXTime& timeInfo) const {

    ExaTrkXTimer tot_timer;
    ExaTrkXTimer timer;
    tot_timer.start();

    // ************
    // Embedding
    // ************
    timer.start();
    int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
    std::vector<int64_t> embedInputShape{numSpacepoints, m_cfg.spacepointFeatures};

    e_client_->ClearInput();
    e_client_->AddInput<float>("INPUT__0", embedInputShape, inputValues);
    std::vector<float> eOutput;
    std::vector<int64_t> embedOutputShape{numSpacepoints, m_cfg.embeddingDim};
    e_client_->GetOutput<float>("OUTPUT__0", eOutput, embedOutputShape);

    // std::cout <<"Embedding space of libtorch the first SP: \n";
    // std::cout << eOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/1) << std::endl;
    // std::cout << std::endl;

    timeInfo.embedding = timer.stopAndGetElapsedTime();
    
    /// ************
    /// Building Edges
    /// ************
    timer.start();
    b_client_->ClearInput();
    b_client_->AddInput<float>("INPUT0", embedOutputShape, eOutput);
    std::vector<int64_t> edgeListData;
    std::vector<int64_t> edgeListShape{2, -1};
    b_client_->GetOutput<int64_t>("OUTPUT0", edgeListData, edgeListShape);

    int64_t numEdges = edgeListData.size() / 2;
    edgeListShape[1] = numEdges;

    timeInfo.building = timer.stopAndGetElapsedTime();

    // ************
    // Filtering
    // ************
    // std::cout << "Get scores for " << numEdges<< " edges." << std::endl;
    
    timer.start();
    f_client_->ClearInput();
    /// <TODO: reuse the embedding inputs?>
    f_client_->AddInput<float>("INPUT__0", embedInputShape, inputValues);
    f_client_->AddInput<int64_t>("INPUT__1", edgeListShape, edgeListData);

    std::vector<float> fOutputData;
    std::vector<int64_t> fOutputShape{numEdges, 1};
    f_client_->GetOutput<float>("OUTPUT__0", fOutputData, fOutputShape);


    // std::cout << "After filtering: " << fOutput.size(0) << " " << fOutput.size(1) << std::endl;
    // std::cout << fOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/9) << std::endl;

    torch::Tensor edgeListCTen = torch::tensor(edgeListData, {torch::kInt64});
    edgeListCTen = edgeListCTen.reshape({2, numEdges});

    torch::Tensor fOutput = torch::tensor(fOutputData, {torch::kFloat32});
    fOutput.sigmoid_();

    torch::Tensor filterMask = fOutput > m_cfg.filterCut;
    at::Tensor edgesAfterF = edgeListCTen.index({Slice(), filterMask});
   
    // edgesAfterF = edgesAfterF.cpu();
    // std::cout << edgesAfterF.sizes().vec()[0]  << std::endl;
    std::vector<int64_t> edgesAfterFiltering;
    std::copy(
        edgesAfterF.data_ptr<int64_t>(),
        edgesAfterF.data_ptr<int64_t>() + edgesAfterF.numel(),
        std::back_inserter(edgesAfterFiltering));

    int64_t numEdgesAfterF = edgesAfterF.size(1);

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

    torch::Tensor gOutput = torch::tensor(gOutputData, {torch::kFloat32});
    gOutput.sigmoid_();

    timeInfo.gnn = timer.stopAndGetElapsedTime();

    // std::cout << "GNN scores for " << gOutput.size(0) << " edges." << std::endl;
    // std::cout << gOutput.slice(0, 0, 5) << std::endl;
    // ************
    // Track Labeling with cugraph::connected_components
    // ************
    timer.start();



    // l_client_->ClearInput();
    // l_client_->AddInput<int64_t>("NUMNODES", {1}, {numSpacepoints});
    // l_client_->AddInput<int64_t>("INPUT0", gEdgeShape, edgesAfterFiltering);
    // l_client_->AddInputFromTorch<float>("INPUT1", gOutput);
    // std::vector<int64_t> trackLabels;
    // std::vector<int64_t> trackLabelsShape{-1, 1};
    // l_client_->GetOutput<int64_t>("OUTPUT0", trackLabels, trackLabelsShape);

    ///*** replace the following block with python-backend.
    using vertex_t = int32_t;
    std::vector<vertex_t> rowIndices;
    std::vector<vertex_t> colIndices;
    std::vector<float> edgeWeights;
    std::vector<vertex_t> trackLabels(numSpacepoints);
    std::copy(
        edgesAfterF.data_ptr<int64_t>(),
        edgesAfterF.data_ptr<int64_t>()+numEdgesAfterF,
        std::back_insert_iterator(rowIndices));
    std::copy(
        edgesAfterF.data_ptr<int64_t>()+numEdgesAfterF,
        edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF+numEdgesAfterF,
        std::back_insert_iterator(colIndices));
    std::copy(
        gOutput.data_ptr<float>(),
        gOutput.data_ptr<float>() + numEdgesAfterF,
        std::back_insert_iterator(edgeWeights));

    weaklyConnectedComponents<int32_t,int32_t,float>(
        numSpacepoints, 
        rowIndices, colIndices, edgeWeights, trackLabels);

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