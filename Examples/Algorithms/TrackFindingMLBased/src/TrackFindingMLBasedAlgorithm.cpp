#include "ActsExamples/TrackFindingMLBased/TrackFindingMLBasedAlgorithm.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <torch/torch.h>
#include <torch/script.h>
using namespace torch::indexing;

#include <grid.h>
#include <insert_points.h>
#include <counting_sort.h>
#include <prefix_sum.h>
#include <find_nbrs.h>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "cuGraph/mmio_read.h"


ActsExamples::TrackFindingMLBasedAlgorithm::TrackFindingMLBasedAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }

  initTrainedModels();
}

ActsExamples::ProcessCode ActsExamples::TrackFindingMLBasedAlgorithm::execute(
  const ActsExamples::AlgorithmContext& ctx) const 
{
  // creating a thread state data structure,
  // storing thread state pointer
  // PyGILState_STATE gstate;
  // gstate = PyGILState_Ensure();

  // Read input data
  const auto& spacepoints =
    ctx.eventStore.get<SimSpacePointContainer>(m_cfg.inputSpacePoints);

  // Convert Input data to a list of size [num_measurements x measurement_features]
  size_t num_spacepoints = spacepoints.size();
  ACTS_INFO("Received " << num_spacepoints << " spacepoints");
  // input info are [idx, r, phi, z, 'cell_count', 'cell_val', 'leta', 'lphi', 'lx', 'ly', 'lz', 'geta', 'gphi']
  size_t num_features = 13; // <TODO> move this to the configuration
  // std::vector<float> inputMatrix(num_spacepoints * num_features, 0.0);
  // // <TODO> Make configurable which sp features to use.
  // // vectorization is applied automatically?
  // // <ERROR> we use <r, phi, z> as inputs!
  // for(size_t idx=0; idx < num_spacepoints; ++idx){
  //   auto sp = spacepoints[idx];
  //   inputMatrix[num_features*idx] = idx;
  //   // cluster position
  //   inputMatrix[num_features*idx+1] = sp.x();
  //   inputMatrix[num_features*idx+2] = sp.y();
  //   inputMatrix[num_features*idx+3] = sp.z();
  //   // cluster shape uses default values for now
  // }

  std::vector<float> inputMatrix{
    21215, 0.0321, 0.8898, -0.0026, 2.0000, 0.3177, 1.6235, 1.1526, 0.0500, 0.1125, 0.3000, 0.3619, 2.8577,
    21254, 0.0333, 0.8894, -0.0028, 5.0000, 0.2932, 1.0829, 0.5124, 0.2000, 0.1125, 0.3000, 0.3072, -2.6099,
    29500, 0.0716, 0.8774, -0.0079, 4.0000, 0.2621, 1.2490, 0.6435, 0.1500, 0.1125, 0.3000, 0.3294, -3.0288,
    29527, 0.0735, 0.8768, -0.0082, 4.0000, 0.3184, 1.2490, 0.6435, 0.1500, 0.1125, 0.3000, 0.3294, -2.8325,
    36548, 0.1164, 0.8631, -0.0142, 4.0000, 0.3258, 1.1633, 0.2742, 0.2000, 0.0563, 0.3000, 0.1554, -2.8365,
    42951, 0.1714, 0.8462, -0.0223, 5.0000, 0.3242, 1.0829, 0.5124, 0.2000, 0.1125, 0.3000, 0.3072, -2.9372,
    42964., 0.1734, 0.8456, -0.0226, 5.0000, 0.3406, 1.0829, 0.5124, 0.2000, 0.1125, 0.3000, 0.3072, -2.8566,
    75014., 0.2608, 0.8188, -0.0354, 3.0000, 3.0000, 0.3980, 1.3734, 0.2400, 1.2000, 0.5000, 1.5145, 2.8893,
    81389., 0.3625, 0.7872, -0.0508, 3.0000, 3.0000, 0.3980, 1.3734, 0.2400, 1.2000, 0.5000, 1.5145, 2.7098,
    81403., 0.3564, 0.7891, -0.0496, 3.0000, 3.0000, 0.3980, 1.3734, 0.2400, 1.2000, 0.5000, 1.5145, 2.8220,
    81905., 0.3653, 0.7863, -0.0510, 3.0000, 3.0000, 0.3980, 1.3734, 0.2400, 1.2000, 0.5000, 1.5145, 2.7098,
    81920., 0.3593, 0.7882, -0.0510, 4.0000, 4.0000, 0.3924, 1.3102, 0.3200, 1.2000, 0.5000, 1.4532, 2.9438,
    88050., 0.5018, 0.7431, -0.0736, 4.0000, 4.0000, 0.3924, 1.3102, 0.3200, 1.2000, 0.5000, 1.4532, 2.7151,
    88051., 0.4961, 0.7450, -0.0724, 5.0000, 5.0000, 0.3857, 1.2490, 0.4000, 1.2000, 0.5000, 1.3859, 2.9011,
    94093., 0.6620, 0.6866, -0.0988, 5.0000, 5.0000, 0.3857, 1.2490, 0.4000, 1.2000, 0.5000, 1.3859, 2.6499,
    94097., 0.6564, 0.6887, -0.0988, 5.0000, 5.0000, 0.3857, 1.2490, 0.4000, 1.2000, 0.5000, 1.3859, 2.7115,
    110072., 0.8173, 0.6219, -0.1300, 7.0000, 7.0000, 0.0646, 1.4932, 0.8400, 10.8000, 0.7000, 2.9859, 2.6896,
    114958., 1.0236, 0.4845, -0.1898, 14.0000, 14.0000, 0.0640, 1.4165, 1.6800, 10.8000, 0.7000, 2.4809, 2.5348,
    114959., 1.0195, 0.4892, -0.1898, 15.0000, 15.0000, 0.0639, 1.4056, 1.8000, 10.8000, 0.7000, 2.4224, 2.6000
  };

  // Convert C++ vector to Python list
  // <TODO> Check memory management
  PyObject *pArgs, *pValue;
  pArgs = PyTuple_New(1);
  pValue = PyList_New(inputMatrix.size());
  vector_to_pylist(inputMatrix, pValue);
  PyTuple_SetItem(pArgs, 0, pValue);

  // call ML-based track finding
  pValue = PyObject_CallObject(_pFunc, pArgs);

  // pValue contains a list of dynamic-size lists
  // Has to perform two loops to exatract the info.
  ProtoTrackContainer protoTracks;

  if (pValue && PyList_Check(pValue)) {
    Py_ssize_t num_trks = PyList_Size(pValue);
    protoTracks.reserve(num_trks);

    PyObject *pTrk, *pSP;
    for(Py_ssize_t i=0; i < num_trks; ++i){
      pTrk = PyList_GetItem(pValue, i);
      auto protoTrack = ProtoTrack();

      if (pTrk && PyList_Check(pTrk)){
        Py_ssize_t num_sps = PyList_Size(pTrk);
        for(Py_ssize_t j=0; j < num_sps; ++j) {
          pSP = PyList_GetItem(pTrk, j); // cannot fail
          if(!PyLong_Check(pSP)) continue; // Skip non-integers
          auto idx_sp = PyLong_AsLong(pSP);
          if (idx_sp == -1 && PyErr_Occurred()){
            // Integer too big to fit in a C long, bail out
            continue;
          }
          protoTrack.push_back(spacepoints[idx_sp].measurementIndex());
        }
      }
      protoTracks.push_back(protoTrack);
    }
  } else {
    ACTS_WARNING("ML-based Track finding failed");
    protoTracks.push_back(ProtoTrack());
  }
  Py_DECREF(pArgs);
  Py_DECREF(pValue);

  // PyGILState_Release(gstate);

  // Py_Finalize(); # there are no finalize function... where should I put this?
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}


void ActsExamples::TrackFindingMLBasedAlgorithm::initTrainedModels() {
    // Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ExaTrkX");
    m_env = new Ort::Env(ORT_LOGGING_LEVEL_WARNING, "ExaTrkX");
    std::string embedModelPath(m_cfg.modelDir + "/embedding.onnx");
    std::string filterModelPath(m_cfg.modelDir + "/filtering.onnx");
    std::string gnnModelPath(m_cfg.modelDir + "/gnn.onnx");
    // <TODO: improve the call to avoid calling copying construtors >

    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);
    OrtStatus* status = OrtSessionOptionsAppendExecutionProvider_CUDA(session_options, 0);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

    e_sess = new Ort::Session(*m_env, embedModelPath.c_str(), session_options);
    f_sess = new Ort::Session(*m_env, filterModelPath.c_str(), session_options);
    g_sess = new Ort::Session(*m_env, gnnModelPath.c_str(), session_options);
}


void ActsExamples::TrackFindingMLBasedAlgorithm::runSessionWithIoBinding(
    Ort::Session& sess,
    std::vector<const char*>& inputNames,
    std::vector<Ort::Value> & inputData,
    std::vector<const char*>& outputNames,
    std::vector<Ort::Value>&  outputData)
{
    // std::cout <<"In the runSessionWithIoBinding" << std::endl;
    if (inputNames.size() < 1) {
        throw std::runtime_error("Onnxruntime input data maping cannot be empty");
    }
    assert(inputNames.size() == inputData.size());

    Ort::IoBinding iobinding(sess);
    for(size_t idx = 0; idx < inputNames.size(); ++idx){
        iobinding.BindInput(inputNames[idx], inputData[idx]);
    }


    for(size_t idx = 0; idx < outputNames.size(); ++idx){
        iobinding.BindOutput(outputNames[idx], outputData[idx]);
    }

    // std::cout <<"Before running onnx" << std::endl;
    sess.Run(Ort::RunOptions{nullptr}, iobinding);

    // std::cout <<"Quitting the runSessionWithIoBinding" << std::endl;
    // return iobinding.GetOutputValues();

}


void ActsExamples::TrackFindingMLBasedAlgorithm::buildEdges(
  std::vector<float>& embedFeatures, std::vector<int64_t>& edgeList, int64_t numSpacepoints){
    torch::Device device(torch::kCUDA);
    auto options = torch::TensorOptions().dtype(torch::kFloat32).device(torch::kCUDA);

    int grid_params_size;
    int grid_delta_idx;
    int grid_total_idx;
    int grid_max_res;
    int grid_dim;
    int dim = m_cfg.embeddingDim;
    if (dim >= 3) {
        grid_params_size = 8;
        grid_delta_idx = 3;
        grid_total_idx = 7;
        grid_max_res = 128;
        grid_dim = 3;
    } else {
        throw std::runtime_error("DIM < 3 is not supported for now.\n");
    }


    float cell_size;
    float radius_cell_ratio = 2.0;
    int G = -1;
    int batch_size = 1;
    float rVal = m_cfg.rVal; // radius of nearest neighours
    int kVal = m_cfg.knnVal;  // maximum number of nearest neighbours.
    
    // Set up grid properties
    torch::Tensor grid_min;
    torch::Tensor grid_max;
    torch::Tensor grid_size;

    torch::Tensor embedTensor = torch::tensor(embedFeatures, options).reshape({1, numSpacepoints, m_cfg.embeddingDim});
    torch::Tensor gridParamsCuda = torch::zeros({batch_size, grid_params_size}, device).to(torch::kFloat32);
    torch::Tensor r_tensor = torch::full({batch_size}, rVal, device);
    torch::Tensor lengths = torch::full({batch_size}, numSpacepoints, device);
    
    
    // build the grid
    for(int i=0; i < batch_size; i++) {
        torch::Tensor allPoints = embedTensor.index({i, Slice(None, lengths.index({i}).item().to<long>()), Slice(None, grid_dim)});
        grid_min = std::get<0>(allPoints.min(0));
        grid_max = std::get<0>(allPoints.max(0));
        gridParamsCuda.index_put_({i, Slice(None, grid_delta_idx)}, grid_min);
        
        grid_size = grid_max - grid_min;
        
        cell_size = r_tensor.index({i}).item().to<float>() / radius_cell_ratio;
        
        if (cell_size < (grid_size.min().item().to<float>() / grid_max_res)) {
            cell_size = grid_size.min().item().to<float>() / grid_max_res;
        }
        
        gridParamsCuda.index_put_({i, grid_delta_idx}, 1 / cell_size);
        
        gridParamsCuda.index_put_({i, Slice(1 + grid_delta_idx, grid_total_idx)},
                                floor(grid_size / cell_size) + 1);
        
        gridParamsCuda.index_put_({i, grid_total_idx}, gridParamsCuda.index({i, Slice(1 + grid_delta_idx, grid_total_idx)}).prod());
        
        if (G < gridParamsCuda.index({i, grid_total_idx}).item().to<int>()) {
            G = gridParamsCuda.index({i, grid_total_idx}).item().to<int>();
        }
    }
    
    torch::Tensor pc_grid_cnt = torch::zeros({batch_size, G}, device).to(torch::kInt32);
    torch::Tensor pc_grid_cell = torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);
    torch::Tensor pc_grid_idx = torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);
    
    std::cout << "Inserting points" << std::endl;
    
    // put spacepoints into the grid
    InsertPointsCUDA(embedTensor, lengths.to(torch::kInt64), gridParamsCuda, 
                    pc_grid_cnt, pc_grid_cell, pc_grid_idx, G);
    
    torch::Tensor pc_grid_off = torch::full({batch_size, G}, 0, device).to(torch::kInt32);
    torch::Tensor grid_params = gridParamsCuda.to(torch::kCPU);
    
    std::cout << "Prefix Sum" << std::endl;
    
    for(int i = 0; i < batch_size ; i++) {
        PrefixSumCUDA(pc_grid_cnt.index({i}), grid_params.index({i, grid_total_idx}).item().to<int>(), pc_grid_off.index({i}));
    }
    
    torch::Tensor sorted_points = torch::zeros({batch_size, numSpacepoints, dim}, device).to(torch::kFloat32);
    torch::Tensor sorted_points_idxs = torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);
    
    CountingSortCUDA(embedTensor, lengths.to(torch::kInt64), pc_grid_cell,
                    pc_grid_idx, pc_grid_off,
                    sorted_points, sorted_points_idxs);
    
    std::cout << "Counting sorted" << std::endl;
    
    // torch::Tensor K_tensor = torch::full({batch_size}, kVal, device); 
    
    std::tuple<at::Tensor, at::Tensor> nbr_output = FindNbrsCUDA(sorted_points, sorted_points,
                                                        lengths.to(torch::kInt64), lengths.to(torch::kInt64),
                                                            pc_grid_off.to(torch::kInt32),
                                                        sorted_points_idxs, sorted_points_idxs,
                                                            gridParamsCuda.to(torch::kFloat32),
                                                        kVal, r_tensor, r_tensor*r_tensor);

    std::cout << "Neigbours to Edges" << std::endl;                             
    torch::Tensor positiveIndices = std::get<0>(nbr_output) >= 0;

    torch::Tensor repeatRange = torch::arange(positiveIndices.size(1), device).repeat({1, positiveIndices.size(2), 1}).transpose(1,2);
    
    torch::Tensor stackedEdges = torch::stack({repeatRange.index({positiveIndices}), std::get<0>(nbr_output).index({positiveIndices})});

    //  Remove self-loops:

    torch::Tensor selfLoopMask = stackedEdges.index({0}) != stackedEdges.index({1});
    stackedEdges = stackedEdges.index({Slice(), selfLoopMask});
    
    // Perform any other post-processing here. E.g. Can remove half of edge list with:
    
    torch::Tensor duplicate_mask = stackedEdges.index({0}) > stackedEdges.index({1});
    stackedEdges = stackedEdges.index({Slice(), duplicate_mask});
    
    // And randomly flip direction with:
    // torch::Tensor random_cut_keep = torch::randint(2, {stackedEdges.size(1)});
    // torch::Tensor random_cut_flip = 1-random_cut_keep;
    // torch::Tensor keep_edges = stackedEdges.index({Slice(), random_cut_keep.to(torch::kBool)});
    // torch::Tensor flip_edges = stackedEdges.index({Slice(), random_cut_flip.to(torch::kBool)}).flip({0});
    // stackedEdges = torch::cat({keep_edges, flip_edges}, 1);
    stackedEdges = stackedEdges.toType(torch::kInt64).to(torch::kCPU);

    std::cout << "copy edges to std::vector" << std::endl;
    std::copy(stackedEdges.data_ptr<int64_t>(), stackedEdges.data_ptr<int64_t>() + stackedEdges.numel(), std::back_inserter(edgeList)); 
}


void ActsExamples::TrackFindingMLBasedAlgorithm::getTracks(
  std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates)
{
    // hardcoded debugging information
    bool debug = true;
    const std::string embedding_outname = "debug_embedding_outputs.txt";
    const std::string edgelist_outname = "debug_edgelist_outputs.txt";
    const std::string filtering_outname = "debug_filtering_scores.txt";


    Ort::AllocatorWithDefaultOptions allocator;
    auto memoryInfo = Ort::MemoryInfo::CreateCpu(
        OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

    // printout the r,phi,z of the first spacepoint
    std::cout <<"First spacepoint information: ";
    std::copy(inputValues.begin(), inputValues.begin() + 3,
              std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;

    // ************
    // Embedding
    // ************

    int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
    std::vector<int64_t> eInputShape{numSpacepoints, m_cfg.spacepointFeatures};

    std::vector<const char*> eInputNames{"sp_features"};
    std::vector<Ort::Value> eInputTensor;
    eInputTensor.push_back(
        Ort::Value::CreateTensor<float>(
            memoryInfo, inputValues.data(), inputValues.size(),
            eInputShape.data(), eInputShape.size())
    );

    std::vector<float> eOutputData(numSpacepoints*m_cfg.embeddingDim);
    std::vector<const char*> eOutputNames{"embedding_output"};
    std::vector<int64_t> eOutputShape{numSpacepoints, m_cfg.embeddingDim};
    std::vector<Ort::Value> eOutputTensor;
    eOutputTensor.push_back(
        Ort::Value::CreateTensor<float>(
            memoryInfo, eOutputData.data(), eOutputData.size(),
            eOutputShape.data(), eOutputShape.size())
    );
    runSessionWithIoBinding(*e_sess, eInputNames, eInputTensor, eOutputNames, eOutputTensor);

    std::cout <<"Embedding space of the first SP: ";
    std::copy(eOutputData.begin(), eOutputData.begin() + m_cfg.embeddingDim,
              std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
    if (debug){
        std::fstream out(embedding_outname, out.out);
        if (!out.is_open()){
            std::cout << "failed to open " << embedding_outname << '\n';
        } else {
            std::copy(eOutputData.begin(), eOutputData.end(),
                    std::ostream_iterator<float>(out, " "));
        }
    }


    // ************
    // Building Edges
    // ************
    std::vector<int64_t> edgeList;
    buildEdges(eOutputData, edgeList, numSpacepoints);
    int64_t numEdges = edgeList.size() / 2;
    std::cout << "Built " << numEdges<< " edges." << std::endl;


    std::copy(edgeList.begin(), edgeList.begin() + 10,
              std::ostream_iterator<int64_t>(std::cout, " "));
    std::cout << std::endl;
    std::copy(edgeList.begin()+numEdges, edgeList.begin()+numEdges+10,
              std::ostream_iterator<int64_t>(std::cout, " "));
    std::cout << std::endl;

    if (debug){
        std::fstream out(edgelist_outname, out.out);
        if (!out.is_open()){
            std::cout << "failed to open " << edgelist_outname << '\n';
        } else {
            std::copy(edgeList.begin(), edgeList.end(),
                    std::ostream_iterator<int64_t>(out, " "));
        }
    }

    // ************
    // Filtering
    // ************
    std::vector<const char*> fInputNames{"f_nodes", "f_edges"};
    std::vector<Ort::Value> fInputTensor;
    fInputTensor.push_back(
        std::move(eInputTensor[0])
    );
    std::vector<int64_t> fEdgeShape{2, numEdges};
    fInputTensor.push_back(
        Ort::Value::CreateTensor<int64_t>(
            memoryInfo, edgeList.data(), edgeList.size(),
            fEdgeShape.data(), fEdgeShape.size())
    );

    // filtering outputs
    std::vector<const char*> fOutputNames{"f_edge_score"};
    std::vector<float> fOutputData(numEdges);
    std::vector<int64_t> fOutputShape{numEdges, 1};
    std::vector<Ort::Value> fOutputTensor;
    fOutputTensor.push_back(
        Ort::Value::CreateTensor<float>(
            memoryInfo, fOutputData.data(), fOutputData.size(), 
            fOutputShape.data(), fOutputShape.size())
    );
    runSessionWithIoBinding(*f_sess, fInputNames, fInputTensor, fOutputNames, fOutputTensor);

    std::cout << "Get scores for " << numEdges<< " edges." << std::endl;
    // However, I have to convert those numbers to a score by applying sigmoid!
    // Use torch::tensor
    torch::Tensor edgeListCTen = torch::tensor(edgeList, {torch::kInt64});
    edgeListCTen = edgeListCTen.reshape({2, numEdges});

    torch::Tensor fOutputCTen = torch::tensor(fOutputData, {torch::kFloat32});
    fOutputCTen = fOutputCTen.sigmoid();

    if (debug){
        std::fstream out(filtering_outname, out.out);
        if (!out.is_open()){
            std::cout << "failed to open " << filtering_outname << '\n';
        } else {
            std::copy(fOutputCTen.data_ptr<float>(), fOutputCTen.data_ptr<float>() + fOutputCTen.numel(),
                    std::ostream_iterator<float>(out, " "));
        }
    }
    // std::cout << fOutputCTen.slice(0, 0, 3) << std::endl;
    torch::Tensor filterMask = fOutputCTen > m_cfg.filterCut;
    torch::Tensor edgesAfterFCTen = edgeListCTen.index({Slice(), filterMask});


    std::vector<int64_t> edgesAfterFiltering;
    std::copy(
        edgesAfterFCTen.data_ptr<int64_t>(),
        edgesAfterFCTen.data_ptr<int64_t>() + edgesAfterFCTen.numel(),
        std::back_inserter(edgesAfterFiltering));

    int64_t numEdgesAfterF = edgesAfterFiltering.size() / 2;
    std::cout << "After filtering: " << numEdgesAfterF << " edges." << std::endl;

    // ************
    // GNN
    // ************
    std::vector<const char*> gInputNames{"g_nodes", "g_edges"};
    std::vector<Ort::Value> gInputTensor;
    gInputTensor.push_back(
        std::move(fInputTensor[0])
    );
    std::vector<int64_t> gEdgeShape{2, numEdgesAfterF};
    gInputTensor.push_back(
        Ort::Value::CreateTensor<int64_t>(
            memoryInfo, edgesAfterFiltering.data(), edgesAfterFiltering.size(),
            gEdgeShape.data(), gEdgeShape.size())
    );
    // gnn outputs
    std::vector<const char*> gOutputNames{"gnn_edge_score"};
    std::vector<float> gOutputData(numEdgesAfterF);
    std::vector<int64_t> gOutputShape{numEdgesAfterF};
    std::vector<Ort::Value> gOutputTensor;
    gOutputTensor.push_back(
        Ort::Value::CreateTensor<float>(
            memoryInfo, gOutputData.data(), gOutputData.size(), 
            gOutputShape.data(), gOutputShape.size())
    );
    runSessionWithIoBinding(*g_sess, gInputNames, gInputTensor, gOutputNames, gOutputTensor);

    torch::Tensor gOutputCTen = torch::tensor(gOutputData, {torch::kFloat32});
    gOutputCTen = gOutputCTen.sigmoid();
    // std::cout << gOutputCTen.slice(0, 0, 3) << std::endl;

    // ************
    // Track Labeling with cugraph::connected_components
    // ************
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

    weakly_connected_components<int32_t,int32_t,float>(
        rowIndices, colIndices, edgeWeights, trackLabels);

    int idx = 0;
    std::cout << "size of components: " << trackLabels.size() << std::endl;
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
}