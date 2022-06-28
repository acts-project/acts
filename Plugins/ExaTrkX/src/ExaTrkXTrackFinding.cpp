#include "Acts/Plugins/ExaTrkX/ExaTrkXTrackFinding.hpp"

#include "Acts/Plugins/ExaTrkX/ExaTrkXUtils.hpp"

#include <torch/script.h>
#include <torch/torch.h>
using namespace torch::indexing;

#include <grid/counting_sort.h>
#include <grid/find_nbrs.h>
#include <grid/grid.h>
#include <grid/insert_points.h>
#include <grid/prefix_sum.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
// #include "mmio_read.h"

namespace {
  void print_current_cuda_meminfo() {
    constexpr int kb = 1024;
    constexpr int mb = kb * kb;
    
    int device;
    std::size_t free, total;
    //c10::cuda::CUDACachingAllocator::emptyCache();
    cudaMemGetInfo( &free, &total );
    cudaGetDevice( &device );
    
    std::cout << "Current CUDA device - used / total [in MB]: " << (total - free) / mb << " / " << total / mb << std::endl;
  }
}

namespace Acts {

ExaTrkXTrackFinding::ExaTrkXTrackFinding(
    const ExaTrkXTrackFinding::Config& config)
    : ExaTrkXTrackFindingBase("ExaTrkXTrackFinding", config.verbose),
      m_cfg(config) {
  initTrainedModels();
}

void ExaTrkXTrackFinding::initTrainedModels() {
  std::string l_embedModelPath(m_cfg.modelDir + "/torchscript/embed.pt");
  std::string l_filterModelPath(m_cfg.modelDir + "/torchscript/filter.pt");
  std::string l_gnnModelPath(m_cfg.modelDir + "/torchscript/gnn.pt");
  c10::InferenceMode guard(true);
  try {
    m_embeddingModel = torch::jit::load(l_embedModelPath.c_str());
    m_embeddingModel.eval();
    m_filterModel = torch::jit::load(l_filterModelPath.c_str());
    m_filterModel.eval();
    m_gnnModel = torch::jit::load(l_gnnModelPath.c_str());
    m_gnnModel.eval();
  } catch (const c10::Error& e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

// The main function that runs the Exa.TrkX ExaTrkXTrackFindingence pipeline
// Be care of sharpe corners.
void ExaTrkXTrackFinding::getTracks(
    std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
    std::vector<std::vector<int> >& trackCandidates,
    ExaTrkXTime& timeInfo) const {
  ExaTrkXTimer tot_timer;
  tot_timer.start();
  // hardcoded debugging information
  c10::InferenceMode guard(true);
  torch::Device device(torch::kCUDA);
  
  // Clone models (solve memory leak? members can be const...)
  auto e_model = m_embeddingModel.clone();
  auto f_model = m_filterModel.clone();
  auto g_model = m_gnnModel.clone();

  // printout the r,phi,z of the first spacepoint
  if (m_cfg.verbose) {
    std::cout << "First spacepoint information: ";
    std::copy(inputValues.begin(), inputValues.begin() + 3,
              std::ostream_iterator<float>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Max and min spacepoint:"
              << *std::max_element(inputValues.begin(), inputValues.end())
              << " "
              << *std::min_element(inputValues.begin(), inputValues.end())
              << "\n";
    print_current_cuda_meminfo();
  }

  ExaTrkXTimer timer;
  
  // **********
  // Embedding
  // **********

  timer.start();
  int64_t numSpacepoints = inputValues.size() / m_cfg.spacepointFeatures;
  std::vector<torch::jit::IValue> eInputTensorJit;
  auto e_opts = torch::TensorOptions().dtype(torch::kFloat32);
  torch::Tensor eLibInputTensor =
      torch::from_blob(inputValues.data(),
                       {numSpacepoints, m_cfg.spacepointFeatures}, e_opts)
          .to(torch::kFloat32);

  eInputTensorJit.push_back(eLibInputTensor.to(device));
  std::optional<at::Tensor> eOutput = e_model.forward(eInputTensorJit).toTensor();
  eInputTensorJit.clear();

  if (m_cfg.verbose) {
    std::cout << "Embedding space of libtorch the first SP: \n";
    std::cout << eOutput->slice(/*dim=*/0, /*start=*/0, /*end=*/1) << std::endl;
    print_current_cuda_meminfo();
  }

  timeInfo.embedding = timer.stopAndGetElapsedTime();

  // ****************
  // Building Edges
  // ****************
  
  timer.start();
  
  std::optional<torch::Tensor> edgeList = buildEdges(
      *eOutput, numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal, m_cfg.knnVal);
  eOutput.reset();
  // torch::Tensor edgeList = buildEdgesBruteForce(
  //   eOutput, numSpacepoints, m_cfg.embeddingDim, m_cfg.rVal,
  //   m_cfg.knnVal);

  if (m_cfg.verbose) {
    std::cout << "Built " << edgeList->size(1) << " edges. " << edgeList->size(0)
              << std::endl;
    std::cout << edgeList->slice(1, 0, 5) << std::endl;
    print_current_cuda_meminfo();

    // std::ofstream file("reconstruction/edges_after_embedding.csv");
    // file << "e0,e1\n";
    // for(int i=0; i<edgeList.size(1); ++i) {
    //     file << edgeList[0][i].item<float>() << "," << edgeList[1][i].item<float>() << "\n";
    // }
  }

  timeInfo.building = timer.stopAndGetElapsedTime();

  // **********
  // Filtering
  // **********
  
  timer.start();
  
  const auto chunks = at::chunk(at::arange(edgeList->size(1)), m_cfg.n_chunks);
  std::vector<at::Tensor> results;
  
  for(const auto &chunk : chunks) {
    std::vector<torch::jit::IValue> fInputTensorJit;
    fInputTensorJit.push_back(eLibInputTensor.to(device));
    fInputTensorJit.push_back(edgeList->index({Slice(), chunk}).to(device));
  
    results.push_back(f_model.forward(fInputTensorJit).toTensor());
    results.back().squeeze_();
    results.back().sigmoid_();
  }
  
  auto fOutput = torch::cat(results);
  results.clear();

  if (m_cfg.verbose) {
    std::cout << "After filtering network: " << fOutput.size(0) << std::endl;
    std::cout << fOutput.slice(/*dim=*/0, /*start=*/0, /*end=*/9) << std::endl;
    print_current_cuda_meminfo();
  }

  torch::Tensor filterMask = fOutput > m_cfg.filterCut;
  torch::Tensor edgesAfterF = edgeList->index({Slice(), filterMask});
  edgeList.reset();
  edgesAfterF = edgesAfterF.to(torch::kInt64);
  int64_t numEdgesAfterF = edgesAfterF.size(1);

  if (m_cfg.verbose) {
    std::cout << "After filter cut: " << numEdgesAfterF << " edges."
              << std::endl;

    // std::ofstream file("reconstruction/edges_after_filter.csv");
    // file << "e0,e1\n";
    // for(int i=0; i<edgesAfterF.size(1); ++i) {
    //     file << edgesAfterF[0][i].item<float>() << "," << edgesAfterF[1][i].item<float>() << "\n";
    // }
    print_current_cuda_meminfo();
  }

  timeInfo.filtering = timer.stopAndGetElapsedTime();

  // ****
  // GNN
  // ****
  
  timer.start();

  auto bidirEdgesAfterF = torch::cat({edgesAfterF, edgesAfterF.flip(0)}, 1);

  if (m_cfg.verbose) {
    std::cout << "bidir edges shape " << bidirEdgesAfterF.size(0) << ", "
              << bidirEdgesAfterF.size(1) << "\n";
    print_current_cuda_meminfo();
  }

  std::vector<torch::jit::IValue> gInputTensorJit;
  // auto g_opts = torch::TensorOptions().dtype(torch::kInt64);
  gInputTensorJit.push_back(eLibInputTensor.to(device));
  gInputTensorJit.push_back(bidirEdgesAfterF.to(device));

  auto gOutputBidir = g_model.forward(gInputTensorJit).toTensor();
  gInputTensorJit.clear();
  gOutputBidir.sigmoid_();
  gOutputBidir = gOutputBidir.cpu();

  auto gOutput = gOutputBidir.index({Slice(None, gOutputBidir.size(0) / 2)});

  timeInfo.gnn = timer.stopAndGetElapsedTime();

  if (m_cfg.verbose) {
    std::cout << "GNN scores for " << gOutput.size(0) << " edges." << std::endl;
    std::cout << "(Bidir scores size: " << gOutputBidir.size(0) << std::endl;
    std::cout << gOutput.slice(0, 0, 5) << std::endl;
    // std::ofstream file("reconstruction/gnn_scores.csv");
    // file << "score\n";
    // for(int i=0; i<edgesAfterF.size(1); ++i) {
    //     file << gOutput[i].item<float>() << "\n";
    // }
    print_current_cuda_meminfo();
  }

  // ***************
  // Track Labeling
  // ***************
  
  timer.start();

  using vertex_t = int32_t;
  std::vector<vertex_t> rowIndices;
  std::vector<vertex_t> colIndices;
  std::vector<float> edgeWeights;
  std::vector<vertex_t> trackLabels(numSpacepoints);
  std::copy(edgesAfterF.data_ptr<int64_t>(),
            edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF,
            std::back_insert_iterator(rowIndices));
  std::copy(edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF,
            edgesAfterF.data_ptr<int64_t>() + numEdgesAfterF + numEdgesAfterF,
            std::back_insert_iterator(colIndices));
  std::copy(gOutput.data_ptr<float>(),
            gOutput.data_ptr<float>() + numEdgesAfterF,
            std::back_insert_iterator(edgeWeights));

  weaklyConnectedComponents<int32_t, int32_t, float>(
      numSpacepoints, rowIndices, colIndices, edgeWeights, trackLabels,
      m_cfg.edgeCut);

  // weakly_connected_components<int32_t,int32_t,float>(
  //     rowIndices, colIndices, edgeWeights, trackLabels);

  if (m_cfg.verbose) {
    std::cout << "size of components: " << trackLabels.size() << std::endl;
    std::vector<vertex_t> sorted(trackLabels);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::cout << "unique components: " << sorted.size() << std::endl;
    print_current_cuda_meminfo();
  }

  if (trackLabels.size() == 0)
    return;

  trackCandidates.clear();

  int existTrkIdx = 0;
  // map labeling from MCC to customized track id.
  std::map<int32_t, int32_t> trackLableToIds;

  for (int32_t idx = 0; idx < numSpacepoints; ++idx) {
    int32_t trackLabel = trackLabels[idx];
    int spacepointID = spacepointIDs[idx];

    int trkId;
    if (trackLableToIds.find(trackLabel) != trackLableToIds.end()) {
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
  c10::cuda::CUDACachingAllocator::emptyCache();
}

}  // namespace Acts
