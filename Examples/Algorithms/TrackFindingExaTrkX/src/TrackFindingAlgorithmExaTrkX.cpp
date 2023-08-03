// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <torch/torch.h>

using namespace ActsExamples;

namespace {

template <typename T>
auto cantor(const T& t) {
  return t[0] + (t[0] + t[1]) * (t[0] + t[1] + 1) / 2;
}

std::ostream &operator<<(std::ostream &os, const std::vector<int64_t> &track) {
   for(auto u : track) {
     os << u << " ";
   }
   return os;
}


}  // namespace

class TruthGraph {
  std::vector<int64_t> m_truthGraphCantor;
  const Acts::Logger &m_logger;

  auto &logger() const { return m_logger; }

 public:
  TruthGraph(const SimSpacePointContainer& spacepoints,
             const IndexMultimap<Index>& measHitMap,
             const SimHitContainer& truthHits, const Acts::Logger &l) : m_logger(l) {
    std::unordered_map<SimBarcode, std::pair<double, std::vector<std::size_t>>> tracks;

    for (auto i = 0ul; i < spacepoints.size(); ++i) {
      const auto measId = spacepoints[i]
                              .sourceLinks()[0]
                              .template get<IndexSourceLink>()
                              .index();

      auto [a, b] = measHitMap.equal_range(measId);
      for (auto it = a; it != b; ++it) {
        const auto &hit = *truthHits.nth(it->second);
        const auto pid = hit.particleId();
        tracks[pid].first = std::max(tracks[pid].first, hit.momentum4Before()[3]); 
        tracks[pid].second.push_back(i);
      }
    }

    // for(const auto &[pid, track] : tracks) {
    //     const auto &[m, t] = track;
    //     ACTS_VERBOSE(pid << ": " << t);
    // }
/*
    for(auto it=tracks.begin(); it != tracks.end();) {
        if( (it->second.first < 0.5) && (it->second.second.size() < 3) ) {
            it = tracks.erase(it);
        } else {
            ++it;
        }
    }
*/
    std::vector<int64_t> truthGraph;
    for (auto& [_, v] : tracks) {
      auto& [mom, track] = v;
      std::sort(track.begin(), track.end());
      for (auto i = 0ul; i < track.size() - 1; ++i) {
        truthGraph.push_back(track[i]);
        truthGraph.push_back(track[i + 1]);
      }
    }

    const auto n = std::min(10ul, truthGraph.size()/2);
    ACTS_INFO("Truth graph edges " << truthGraph.size() / 2);
    auto opts = torch::TensorOptions().dtype(torch::kInt64);
    ACTS_INFO("slice:\n" << torch::from_blob(truthGraph.data(),
                            {static_cast<long>(truthGraph.size() / 2), 2}, opts)
        .clone()
        .transpose(0, 1).slice(1, 0, n));

    m_truthGraphCantor.reserve(truthGraph.size() / 2);
    for (auto it = truthGraph.begin(); it != truthGraph.end(); it += 2) {
      m_truthGraphCantor.push_back(cantor(it));
    }
    std::sort(m_truthGraphCantor.begin(), m_truthGraphCantor.end());

    ACTS_INFO("Truth graph cantor:");
    for(auto i=0ul; i<n; ++i) { std::cout << m_truthGraphCantor[i] << " "; }
    std::cout << std::endl;
  }

  std::pair<float, float> effpur(const torch::Tensor& graph) const {
    const auto cantorTensor = cantor(std::get<0>(torch::sort(graph, 0)))
                                  .to(torch::kCPU)
                                  .to(torch::kInt64);

    std::vector<int64_t> predGraphCantor(
        cantorTensor.data_ptr<int64_t>(),
        cantorTensor.data_ptr<int64_t>() + cantorTensor.numel());
    std::sort(predGraphCantor.begin(), predGraphCantor.end());

    const auto n = std::min(10ul, predGraphCantor.size());
    ACTS_INFO("Slice of edgelist:\n" << std::get<0>(torch::sort(graph, 0)).slice(1, 0, n));

    ACTS_INFO("Pred graph cantor:");
    for(auto i=0ul; i<n; ++i) { std::cout << predGraphCantor[i] << " "; }
    std::cout << std::endl;

    ACTS_INFO("Pred graph edges " << predGraphCantor.size());

    std::vector<int64_t> intersection(std::max(predGraphCantor.size(), m_truthGraphCantor.size()));
    auto intersection_end =
        std::set_intersection(predGraphCantor.begin(), predGraphCantor.end(),
                              m_truthGraphCantor.begin(),
                              m_truthGraphCantor.end(), intersection.begin());

    float intersected = std::distance(intersection.begin(), intersection_end);

    ACTS_INFO("Intersection size " << intersected);
    float eff = intersected / m_truthGraphCantor.size();
    float pur = intersected / predGraphCantor.size();

    return {eff, pur};
  }
};

ActsExamples::TrackFindingAlgorithmExaTrkX::TrackFindingAlgorithmExaTrkX(
    Config config, Acts::Logging::Level level)
    : ActsExamples::IAlgorithm("TrackFindingMLBasedAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing spacepoint input collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing protoTrack output collection");
  }
  if (!m_cfg.graphConstructor) {
    throw std::invalid_argument("Missing graph construction module");
  }
  if (!m_cfg.trackBuilder) {
    throw std::invalid_argument("Missing track building module");
  }
  if (m_cfg.edgeClassifiers.empty() or
      not std::all_of(m_cfg.edgeClassifiers.begin(),
                      m_cfg.edgeClassifiers.end(),
                      [](const auto& a) { return static_cast<bool>(a); })) {
    throw std::invalid_argument("Missing graph construction module");
  }

  // Sanitizer run with dummy input to detect configuration issues
  // TODO This would be quite helpful I think, but currently it does not work in
  // general because the stages do not expose the number of node features.
  // However, this must be addressed anyways when we also want to allow to
  // configure this more flexible with e.g. cluster information as input. So for
  // now, we disable this.
#if 0
  if( m_cfg.sanitize ) {
  Eigen::VectorXf dummyInput = Eigen::VectorXf::Random(3 * 15);
  std::vector<float> dummyInputVec(dummyInput.data(),
                                   dummyInput.data() + dummyInput.size());
  std::vector<int> spacepointIDs;
  std::iota(spacepointIDs.begin(), spacepointIDs.end(), 0);
  
  runPipeline(dummyInputVec, spacepointIDs);
  }
#endif

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);

  m_inputSimHits.maybeInitialize(m_cfg.inputSimhits);
  m_inputMeasurementMap.maybeInitialize(m_cfg.inputMeasurementSimhitMap);
}

std::vector<std::vector<int>>
ActsExamples::TrackFindingAlgorithmExaTrkX::runPipeline(
    boost::multi_array<float, 2>& features, std::vector<int>& spacepointIDs,
    const TruthGraph* truthGraph) const {
  auto [nodes, edges] = (*m_cfg.graphConstructor)(features);
  std::any edge_weights;

  if (truthGraph) {
    auto [eff, pur] = truthGraph->effpur(std::any_cast<torch::Tensor>(edges));
    ACTS_INFO("After Graph construction: eff=" << eff << ", pur=" << pur);
  }

  for (auto edgeClassifier : m_cfg.edgeClassifiers) {
    auto [newNodes, newEdges, newWeights] =
        (*edgeClassifier)(std::move(nodes), std::move(edges));
    nodes = std::move(newNodes);
    edges = std::move(newEdges);
    edge_weights = std::move(newWeights);

    if (truthGraph) {
      auto [eff, pur] = truthGraph->effpur(std::any_cast<torch::Tensor>(edges));
      ACTS_INFO("After Graph construction: eff=" << eff << ", pur=" << pur);
    }
  }

  return (*m_cfg.trackBuilder)(std::move(nodes), std::move(edges),
                               std::move(edge_weights), spacepointIDs);
}

enum feat : std::size_t {
  eR,
  ePhi,
  eZ,
  eCellCount,
  eCellSum,
  eClusterX,
  eClusterY
};

ActsExamples::ProcessCode ActsExamples::TrackFindingAlgorithmExaTrkX::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  auto spacepoints = m_inputSpacePoints(ctx);
  std::sort(spacepoints.begin(), spacepoints.end(), [](const auto &a, const auto &b){ return std::hypot(a.x(), a.y()) < std::hypot(b.x(), b.y()); });

  std::optional<ClusterContainer> clusters;
  if (m_inputClusters.isInitialized()) {
    clusters = m_inputClusters(ctx);
  }

  std::optional<TruthGraph> truthGraph;
  TruthGraph* truthGraphPtr = nullptr;
  if (m_inputSimHits.isInitialized() && m_inputMeasurementMap.isInitialized()) {
    truthGraph.emplace(spacepoints, m_inputMeasurementMap(ctx),
                            m_inputSimHits(ctx), logger());
    truthGraphPtr = &*truthGraph;
  }

  // Convert Input data to a list of size [num_measurements x
  // measurement_features]
  const std::size_t num_spacepoints = spacepoints.size();
  const std::size_t numFeatures = clusters ? 7 : 3;
  ACTS_INFO("Received " << num_spacepoints << " spacepoints");

  boost::multi_array<float, 2> features(
      std::array<std::size_t, 2>{num_spacepoints, numFeatures});
  std::vector<int> spacepointIDs;

  spacepointIDs.reserve(spacepoints.size());

  double sumCells = 0.0;
  double sumActivation = 0.0;

  ACTS_INFO("First spacepoint x,y,z=" << spacepoints[0].x() << ", " << spacepoints[0].y() << ", " << spacepoints[0].z());

  for (auto i = 0ul; i < num_spacepoints; ++i) {
    const auto& sp = spacepoints[i];
    // For now just take the first index since does require one single index per
    // spacepoint
    const auto& sl = sp.sourceLinks()[0].template get<IndexSourceLink>();
    spacepointIDs.push_back(sl.index());

    features[i][eR] = std::hypot(sp.x(), sp.y()) / m_cfg.rScale;
    features[i][ePhi] = std::atan2(sp.y(), sp.x()) / m_cfg.phiScale;
    features[i][eZ] = sp.z() / m_cfg.zScale;

    if (clusters) {
      const auto& cluster = clusters->at(sl.index());
      const auto& chnls = cluster.channels;

      features[i][eCellCount] = cluster.channels.size() / m_cfg.cellCountScale;
      features[i][eCellSum] =
          std::accumulate(chnls.begin(), chnls.end(), 0.0,
                          [](double s, const Cluster::Cell& c) {
                            return s + c.activation;
                          }) /
          m_cfg.cellSumScale;
      features[i][eClusterX] = cluster.sizeLoc0 / m_cfg.clusterXScale;
      features[i][eClusterY] = cluster.sizeLoc1 / m_cfg.clusterYScale;

      sumCells += features[i][eCellCount];
      sumActivation += features[i][eCellSum];
    }
  }

  ACTS_DEBUG("Avg cell count: " << sumCells / spacepoints.size());
  ACTS_DEBUG("Avg activation: " << sumActivation / sumCells);

  // Run the pipeline
  const auto trackCandidates =
      runPipeline(features, spacepointIDs, truthGraphPtr);

  ACTS_DEBUG("Done with pipeline");

  // Make the prototracks
  std::vector<ProtoTrack> protoTracks;
  protoTracks.reserve(trackCandidates.size());

  int nShortTracks = 0;

  for (auto& x : trackCandidates) {
    ProtoTrack onetrack;
    std::copy(x.begin(), x.end(), std::back_inserter(onetrack));

    if (onetrack.size() < 3) {
      nShortTracks++;
      continue;
    }
    protoTracks.push_back(std::move(onetrack));
  }

  ACTS_INFO("Removed " << nShortTracks << " with less then 3 hits");
  ACTS_INFO("Created " << protoTracks.size() << " proto tracks");
  m_outputProtoTracks(ctx, std::move(protoTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
