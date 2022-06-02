#include "ExaTrkXTrackFindingTritonPython.hpp"
#include "ExaTrkXUtils.hpp"

#include "grpc_client.h"
#include "grpc_service.pb.h"
namespace tc = triton::client;


ExaTrkXTrackFindingTritonPython::ExaTrkXTrackFindingTritonPython(
    const ExaTrkXTrackFindingTritonPython::Config& config):
    ExaTrkXTrackFindingBase("ExaTrkXTrackFindingTritonPython", config.verbose), m_cfg(config)
{

    std::string l_embedModelPath(m_cfg.modelDir + "/torchscript/embed.pt");
    std::string l_filterModelPath(m_cfg.modelDir + "/torchscript/filter.pt");
    std::string l_gnnModelPath(m_cfg.modelDir + "/torchscript/gnn.pt");
    c10::InferenceMode guard(true);
    try {   
        e_model = torch::jit::load(l_embedModelPath.c_str());
        e_model.eval();
        f_model = torch::jit::load(l_filterModelPath.c_str());
        f_model.eval();
        g_model = torch::jit::load(l_gnnModelPath.c_str());
        g_model.eval();
    } catch (const c10::Error& e) {
        throw std::invalid_argument("Failed to load models: " + e.msg()); 
    }

    bool verbose = false;
    uint32_t client_timeout = 0;
    std::string model_version = "";
    b_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.buildingModelName, m_cfg.url, model_version, client_timeout, verbose);
    l_client_ = std::make_unique<ExaTrkXTriton>(m_cfg.labelingModelName, m_cfg.url, model_version, client_timeout, verbose);
}

// The main function that runs the Exa.TrkX ExaTrkXTrackFindingence pipeline
// Be care of sharpe corners.
void ExaTrkXTrackFindingTritonPython::getTracks(
    std::vector<float>& inputValues,
    std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates,
    ExaTrkXTime& timeInfo) const {

    ExaTrkXTimer tot_timer;
    tot_timer.start();
    // hardcoded debugging information
    c10::InferenceMode guard(true);
    bool debug = true;
    torch::Device device(torch::kCUDA);

    /// printout the r,phi,z of the first spacepoint
    // std::cout <<"First spacepoint information: ";
    // std::copy(inputValues.begin(), inputValues.begin() + 3,
    //           std::ostream_iterator<float>(std::cout, " "));
    // std::cout << std::endl;

    ExaTrkXTimer timer;
    // ************
    // Embedding
    // ************

    timer.start();
    int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
    std::vector<torch::jit::IValue> eInputTensorJit;
    auto e_opts = torch::TensorOptions().dtype(torch::kFloat32);
    torch::Tensor eLibInputTensor = torch::from_blob(
        inputValues.data(),
        {numSpacepoints, m_cfg.spacepointFeatures},
        e_opts).to(torch::kFloat32);

    eInputTensorJit.push_back(eLibInputTensor.to(device));
    at::Tensor eOutput = e_model.forward(eInputTensorJit).toTensor();
    // std::cout <<"Embedding space of libtorch the first SP: \n";
    // std::cout << eOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/1) << std::endl;
    // std::cout << std::endl;

    timeInfo.embedding = timer.stopAndGetElapsedTime();
    
    /// ************
    /// Building Edges
    /// ************
    timer.start();
    b_client_->ClearInput();
    eOutput = eOutput.cpu();
    b_client_->AddInputFromTorch<float>("INPUT0", eOutput);
    std::vector<int64_t> edgeListData;
    std::vector<int64_t> edgeListShape{2, -1};
    b_client_->GetOutput<int64_t>("OUTPUT0", edgeListData, edgeListShape);

    int64_t numEdges = edgeListData.size() / 2;
    edgeListShape[1] = numEdges;
    auto edges_opts = torch::TensorOptions().dtype(torch::kInt64);
    auto edgeList = torch::from_blob(edgeListData.data(), edgeListShape, edges_opts);

    // torch::Tensor edgeList = buildEdges(
    //     eOutput, numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);
    // int64_t numEdges = edgeList.size(1);

    // std::cout << "Built " << edgeList.size(1) << " edges. " <<  edgeList.size(0) << std::endl;
    // std::cout << edgeList.slice(1, 0, 5) << std::endl;

    timeInfo.building = timer.stopAndGetElapsedTime();

    // ************
    // Filtering
    // ************
    // std::cout << "Get scores for " << numEdges<< " edges." << std::endl;
    
    timer.start();
    std::vector<torch::jit::IValue> fInputTensorJit;
    fInputTensorJit.push_back(eLibInputTensor.to(device));
    fInputTensorJit.push_back(edgeList.to(device));
    at::Tensor fOutput = f_model.forward(fInputTensorJit).toTensor();
    fOutput.squeeze_();
    fOutput.sigmoid_();

    // std::cout << "After filtering: " << fOutput.size(0) << " " << fOutput.size(1) << std::endl;
    // std::cout << fOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/9) << std::endl;

    torch::Tensor filterMask = fOutput > m_cfg.filterCut;
    torch::Tensor edgesAfterF = edgeList.index({Slice(), filterMask});
    edgesAfterF = edgesAfterF.to(torch::kInt64);
    int64_t numEdgesAfterF = edgesAfterF.size(1);

    // std::cout << "After filtering: " << numEdgesAfterF << " edges." << std::endl;

    timeInfo.filtering = timer.stopAndGetElapsedTime();

    // ************
    // GNN
    // ************
    timer.start();

    std::vector<torch::jit::IValue> gInputTensorJit;
    auto g_opts = torch::TensorOptions().dtype(torch::kInt64);
    gInputTensorJit.push_back(eLibInputTensor.to(device));
    gInputTensorJit.push_back(edgesAfterF.to(device));
    auto gOutput = g_model.forward(gInputTensorJit).toTensor();
    gOutput.sigmoid_();
    gOutput = gOutput.cpu();
    timeInfo.gnn = timer.stopAndGetElapsedTime();

    // std::cout << "GNN scores for " << gOutput.size(0) << " edges." << std::endl;
    // std::cout << gOutput.slice(0, 0, 5) << std::endl;
    // ************
    // Track Labeling with cugraph::connected_components
    // ************
    timer.start();

    ///*** replace the following block with python-backend.
    // using vertex_t = int32_t;
    // std::vector<vertex_t> rowIndices;
    // std::vector<vertex_t> colIndices;
    // std::vector<float> edgeWeights;
    // std::vector<vertex_t> trackLabels(numSpacepoints);
    // std::copy(
    //     edgesAfterF.data_ptr<int64_t>(),
    //     edgesAfterF.data_ptr<int64_t>()+numEdgesAfterF,
    //     std::back_insert_iterator(rowIndices));
    // std::copy(
    //     edgesAfterF.data_ptr<int64_t>()+numEdgesAfterF,
    //     edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF+numEdgesAfterF,
    //     std::back_insert_iterator(colIndices));
    // std::copy(
    //     gOutput.data_ptr<float>(),
    //     gOutput.data_ptr<float>() + numEdgesAfterF,
    //     std::back_insert_iterator(edgeWeights));

    // weaklyConnectedComponents<int32_t,int32_t,float>(
    //     numSpacepoints, 
    //     rowIndices, colIndices, edgeWeights, trackLabels);
    ///**********************

    // std::cout << "size of components: " << trackLabels.size() << std::endl;
    ///***********************************************************
    
    // The following two hard copy are needed in order to produce sensible results!
    edgesAfterF = edgesAfterF.cpu();
    std::vector<int64_t> gEdgeShape{2, numEdgesAfterF};
    std::vector<int64_t> gOutputShape{numEdgesAfterF};
    std::vector<int64_t> edgesAfterFiltering;
    std::copy(
        edgesAfterF.data_ptr<int64_t>(),
        edgesAfterF.data_ptr<int64_t>() + edgesAfterF.numel(),
        std::back_inserter(edgesAfterFiltering));

    std::vector<float> gOutputData;
    std::copy(
        gOutput.data_ptr<float>(),
        gOutput.data_ptr<float>() + gOutput.numel(),
        std::back_inserter(gOutputData));

    l_client_->ClearInput();
    std::vector<int64_t> numSpacepointsV{numSpacepoints};
    // l_client_->AddInput<int64_t>("NUMNODES", {1}, numSpacepointsV);
    l_client_->AddInput<int64_t>("INPUT0", gEdgeShape, edgesAfterFiltering);
    l_client_->AddInput<float>("INPUT1", gOutputShape, gOutputData);
    std::vector<int64_t> trackLabels;
    std::vector<int64_t> trackLabelsShape{-1, 1};
    l_client_->GetOutput<int64_t>("OUTPUT0", trackLabels, trackLabelsShape);
    ///***********************************************************
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