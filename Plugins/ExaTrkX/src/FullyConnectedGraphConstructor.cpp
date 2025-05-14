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

#ifndef ACTS_EXATRKX_CPUONLY
#include <c10/cuda/CUDAGuard.h>
#endif

namespace Acts {

std::tuple<std::any, std::any, std::any>
FullyConnectedGraphConstructor::operator()(
    std::vector<float> &inputValues, std::size_t numNodes,
    const std::vector<std::uint64_t> & /*moduleIds*/,
    const ExecutionContext &execContext) {
  if (inputValues.empty()) {
    throw NoEdgesError{};
  }

  const auto device =
      execContext.device.type == Acts::Device::Type::eCUDA
          ? torch::Device(torch::kCUDA, execContext.device.index)
          : torch::kCPU;

  // Not sure why I must reset this here...
  auto lastCudaError = cudaGetLastError();
  ACTS_DEBUG("Retrieved last CUDA error: " << lastCudaError);
  c10::InferenceMode guard(true);

  // add a protection to avoid calling for kCPU
#ifdef ACTS_EXATRKX_CPUONLY
  assert(device == torch::Device(torch::kCPU));
#else
  std::optional<c10::cuda::CUDAGuard> device_guard;
  // At least under torch 2.3 and below stream guard causes a memory leak I
  // think We instead just synchronize the stream and use the default torch
  // stream
  // std::optional<c10::cuda::CUDAStreamGuard> streamGuard;
  if (execContext.device.type == Acts::Device::Type::eCUDA) {
    device_guard.emplace(execContext.device.index);
    // streamGuard.emplace(execContext.stream.value());
    execContext.stream->synchronize();
  }
#endif

  torch::NoGradGuard noGradGuard;

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
  std::vector<std::int64_t> edgeListVector;
  edgeListVector.reserve(numEdges * 2);
  std::transform(edgeData.begin(), edgeData.end(),
                 std::back_inserter(edgeListVector),
                 [](const auto &edge) { return std::get<0>(edge); });
  std::transform(edgeData.begin(), edgeData.end(),
                 std::back_inserter(edgeListVector),
                 [](const auto &edge) { return std::get<1>(edge); });

  // Convert the edge list to a tensor
  auto edgeList =
      detail::vectorToTensor2D(edgeListVector, numEdges).contiguous();
  assert(edgeList.size(0) == 2);
  assert(edgeList.size(1) == static_cast<std::int64_t>(numEdges));

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
  std::vector<float> edgeFeatureVector;
  edgeFeatureVector.reserve(numEdgeFeatures * edgeList.size(1));
  for (auto i = 0; i < edgeList.size(1); ++i) {
    auto src = edgeList.index({0, i}).item<int>();
    auto dst = edgeList.index({1, i}).item<int>();

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

    for (auto f : {deltaR, deltaPhi, deltaZ, deltaEta, phislope, rphislope}) {
      edgeFeatureVector.push_back(f);
    }
  }

  auto edgeFeatures =
      detail::vectorToTensor2D(edgeFeatureVector, numEdgeFeatures);
  assert(edgeFeatures.size(0) == static_cast<std::int64_t>(numEdges));
  assert(edgeFeatures.size(1) == numEdgeFeatures);

  auto inputTensor = detail::vectorToTensor2D(inputValues, numAllFeatures);
  assert(inputTensor.size(0) == static_cast<std::int64_t>(numNodes));
  assert(inputTensor.size(1) == static_cast<std::int64_t>(numAllFeatures));

  ACTS_DEBUG("Move data to " << execContext.device);

  auto inputTensorCuda = inputTensor.to(device);
  auto edgeListCuda = edgeList.to(device);
  auto edgeFeaturesCuda = edgeFeatures.to(device);

  ACTS_VERBOSE("inputTensor: " << inputTensorCuda);
  ACTS_VERBOSE("edgeList: " << edgeListCuda);
  ACTS_VERBOSE("edgeFeatures: " << edgeFeaturesCuda);

  return {inputTensorCuda, edgeListCuda, edgeFeaturesCuda};
}

}  // namespace Acts
