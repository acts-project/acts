// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <Acts/Plugins/ExaTrkX/FullyConnectedGraphConstructor.hpp>
#include <Acts/Plugins/ExaTrkX/detail/TensorVectorConversion.hpp>

#ifdef ACTS_EXATRKX_WITH_CUDA
#include <cuda_runtime.h>
#endif

namespace Acts {

PipelineTensors FullyConnectedGraphConstructor::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> & /*moduleIds*/,
    const ExecutionContext &execContext) {
  if (inputValues.empty()) {
    throw NoEdgesError{};
  }

  // Not sure why I must reset this here...
#ifdef ACTS_EXATRKX_WITH_CUDA
  auto lastCudaError = cudaGetLastError();
  ACTS_DEBUG("Retrieved last CUDA error: " << lastCudaError);
#endif

  auto numAllFeatures = inputValues.size() / numNodes;

  // throw if the graph is too large
  std::size_t numEdges = numNodes * (numNodes - 1) / 2;
  if (numEdges > m_cfg.maxGraphSize) {
    ACTS_WARNING("Fully connected graph is larger than configured max edges: "
                 << numEdges << " > " << m_cfg.maxGraphSize);
    throw NoEdgesError();
  }

  // Build fully connected edges
  std::vector<std::tuple<std::int64_t, std::int64_t, float>> edgeData;
  edgeData.reserve(numEdges);
  std::size_t skipped = 0;
  for (auto i = 0ul; i < numNodes; ++i) {
    for (auto j = i + 1; j < numNodes; ++j) {
      auto ri = inputValues.at(i * numAllFeatures + m_cfg.rOffset);
      auto rj = inputValues.at(j * numAllFeatures + m_cfg.rOffset);
      auto zi = inputValues.at(i * numAllFeatures + m_cfg.zOffset);
      auto zj = inputValues.at(j * numAllFeatures + m_cfg.zOffset);

      auto dr = std::abs(ri - rj) * m_cfg.rScale;
      auto dz = std::abs(zi - zj) * m_cfg.zScale;
      if (dr > m_cfg.maxDeltaR || dz > m_cfg.maxDeltaZ) {
        skipped++;
        continue;
      }

      auto dist = std::sqrt(dr * dr + dz * dz);

      if (ri < rj) {
        edgeData.emplace_back(i, j, dist);
      } else {
        edgeData.emplace_back(j, i, dist);
      }
    }
  }

  if (m_cfg.maxOutEdges < std::numeric_limits<std::size_t>::max()) {
    ACTS_DEBUG("Filtering edges to only keep the "
               << m_cfg.maxOutEdges << " closest edges per node");

    // Sort the edges by src node and then by distance
    std::sort(edgeData.begin(), edgeData.end(),
              [](const auto &lhs, const auto &rhs) {
                if (std::get<0>(lhs) == std::get<0>(rhs)) {
                  return std::get<2>(lhs) < std::get<2>(rhs);
                }
                return std::get<0>(lhs) < std::get<0>(rhs);
              });

    // Only keep the N shortest edges per src node
    std::vector<std::size_t> nOutEdges(numNodes, 0);
    edgeData.erase(
        std::remove_if(edgeData.begin(), edgeData.end(),
                       [&nOutEdges, this](const auto &edge) {
                         auto src = std::get<0>(edge);
                         nOutEdges.at(src)++;
                         if (nOutEdges.at(src) <= m_cfg.maxOutEdges) {
                           return false;
                         }
                         return true;
                       }),
        edgeData.end());
  }

  // Check if we have any edges
  if (edgeData.empty()) {
    ACTS_WARNING("No edges created, skipping graph construction");
    throw NoEdgesError{};
  }

  numEdges = edgeData.size();
  ACTS_DEBUG("Built " << numEdges << " edges, skipped " << skipped);

  // Bring the edge list into the right format
  auto edgeList =
      Acts::Tensor<std::int64_t>::Create({2, numEdges}, {Device::Cpu(), {}});
  std::transform(edgeData.begin(), edgeData.end(), edgeList.data(),
                 [](const auto &edge) { return std::get<0>(edge); });
  std::transform(edgeData.begin(), edgeData.end(), edgeList.data() + numEdges,
                 [](const auto &edge) { return std::get<1>(edge); });

  // TODO I think this is already somewhere in the codebase
  const float pi = std::numbers::pi_v<float>;
  auto resetAngle = [pi](float angle) {
    if (angle > pi) {
      return angle - 2.f * pi;
    }
    if (angle < -pi) {
      return angle + 2.f * pi;
    }
    return angle;
  };

  // Shorthand for the feature offsets
  const auto r = m_cfg.rOffset;
  const auto phi = m_cfg.phiOffset;
  const auto z = m_cfg.zOffset;
  const auto eta = m_cfg.etaOffset;

  // TODO Unify edge feature building, this is only to get it in fast
  constexpr static std::size_t numEdgeFeatures = 6;
  auto edgeFeatures = Acts::Tensor<float>::Create({numEdges, numEdgeFeatures},
                                                  {Device::Cpu(), {}});
  for (auto i = 0ul; i < numEdges; ++i) {
    auto src = *(edgeList.data() + i);
    auto dst = *(edgeList.data() + numEdges + i);

    // Edge features
    // See
    // https://gitlab.cern.ch/gnn4itkteam/acorn/-/blob/dev/acorn/utils/loading_utils.py?ref_type=heads#L288
    const float *srcFeatures = inputValues.data() + src * numAllFeatures;
    const float *dstFeatures = inputValues.data() + dst * numAllFeatures;

    const float deltaR = dstFeatures[r] - srcFeatures[r];
    const float deltaPhi =
        resetAngle((dstFeatures[phi] - srcFeatures[phi]) * m_cfg.phiScale) /
        m_cfg.phiScale;
    const float deltaZ = dstFeatures[z] - srcFeatures[z];
    const float deltaEta = dstFeatures[eta] - srcFeatures[eta];
    float phislope = 0.0;
    float rphislope = 0.0;

    if (deltaR != 0.0) {
      phislope = std::clamp(deltaPhi / deltaR, -100.f, 100.f);
      float avgR = 0.5f * (dstFeatures[r] + srcFeatures[r]);
      rphislope = avgR * phislope;
    }

    *(edgeFeatures.data() + i * numEdgeFeatures + 0) = deltaR;
    *(edgeFeatures.data() + i * numEdgeFeatures + 1) = deltaPhi;
    *(edgeFeatures.data() + i * numEdgeFeatures + 2) = deltaZ;
    *(edgeFeatures.data() + i * numEdgeFeatures + 3) = deltaEta;
    *(edgeFeatures.data() + i * numEdgeFeatures + 4) = phislope;
    *(edgeFeatures.data() + i * numEdgeFeatures + 5) = rphislope;
  }

  auto nodeFeatures = Acts::Tensor<float>::Create({numNodes, numAllFeatures},
                                                  {Device::Cpu(), {}});
  std::copy(inputValues.begin(), inputValues.end(), nodeFeatures.data());

  ACTS_DEBUG("Move data to " << execContext.device);

  auto edgeListDevice = edgeList.clone(execContext);
  auto edgeFeaturesDevice = edgeFeatures.clone(execContext);
  auto nodeFeaturesDevice = nodeFeatures.clone(execContext);

  return {std::move(nodeFeaturesDevice),
          std::move(edgeListDevice),
          std::move(edgeFeaturesDevice),
          {}};
}

}  // namespace Acts
