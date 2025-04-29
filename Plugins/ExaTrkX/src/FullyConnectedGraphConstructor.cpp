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

std::tuple<std::any, std::any, std::any> FullyConnectedGraphConstructor::operator()(std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<std::uint64_t> &/*moduleIds*/,
      const ExecutionContext &execContext) {
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
  if (execContext.device.is_cuda()) {
    device_guard.emplace(execContext.device.index());
    // streamGuard.emplace(execContext.stream.value());
    execContext.stream->synchronize();
  }
#endif

  torch::NoGradGuard noGradGuard;


  auto numAllFeatures = inputValues.size() / numNodes;

  // throw if the graph is too large
  std::size_t numEdges = numNodes * (numNodes - 1) / 2;
  if (numEdges > m_cfg.maxGraphSize) {
    ACTS_WARNING("Fully connected graph is larger than configured max edges: " << numEdges
                  << " > " << m_cfg.maxGraphSize);
    throw NoEdgesError();
  }

  // Build fully connected edges
  std::vector<int64_t> edgeListVector;
  edgeListVector.reserve(numEdges * 2);
  std::size_t skipped = 0;
  for(auto i = 0ul; i < numNodes; ++i) {
    for(auto j = i + 1; j < numNodes; ++j) {
      auto ri = inputValues.at(i * numAllFeatures + m_cfg.rOffset);
      auto rj = inputValues.at(j * numAllFeatures + m_cfg.rOffset);
      if( std::abs(ri - rj) * m_cfg.rScale > m_cfg.maxDeltaR ) {
        skipped++;
        continue;
      }

      if( ri < rj) {
        edgeListVector.push_back(i);
        edgeListVector.push_back(j);
      } else {
        edgeListVector.push_back(j);
        edgeListVector.push_back(i);
      }
    } 
  }

  numEdges = edgeListVector.size() / 2;
  ACTS_DEBUG("Built " << numEdges << " edges, skipped " << skipped);

  auto edgeList = detail::vectorToTensor2D(edgeListVector, 2).transpose(0, 1).contiguous();
  assert(edgeList.size(0) == 2);
  assert(edgeList.size(1) == static_cast<int64_t>(numEdges));

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

    const float deltaR = dstFeatures[0] - srcFeatures[0];
    const float deltaPhi =
        resetAngle((dstFeatures[1] - srcFeatures[1]) * m_cfg.phiScale) /
        m_cfg.phiScale;
    const float deltaZ = dstFeatures[2] - srcFeatures[2];
    const float deltaEta = dstFeatures[3] - srcFeatures[3];
    float phislope = 0.0;
    float rphislope = 0.0;

  if (deltaR != 0.0) {
    phislope = std::clamp(deltaPhi / deltaR, -100.f, 100.f);
    float avgR = 0.5f * (dstFeatures[0] + srcFeatures[0]);
    rphislope = avgR * phislope;
  }

    for (auto f : {deltaR, deltaPhi, deltaZ, deltaEta, phislope, rphislope}) {
      edgeFeatureVector.push_back(f);
    }
  }

  auto edgeFeatures =
      detail::vectorToTensor2D(edgeFeatureVector, numEdgeFeatures);
  assert(edgeFeatures.size(0) == static_cast<int64_t>(numEdges));
  assert(edgeFeatures.size(1) == numEdgeFeatures);

  auto inputTensor = detail::vectorToTensor2D(inputValues, numAllFeatures);
  assert(inputTensor.size(0) == static_cast<int64_t>(numNodes));
  assert(inputTensor.size(1) == static_cast<int64_t>(numAllFeatures));

  ACTS_DEBUG("Move data to " << execContext.device);

  auto inputTensorCuda = inputTensor.to(execContext.device);
  auto edgeListCuda = edgeList.to(execContext.device);
  auto edgeFeaturesCuda = edgeFeatures.to(execContext.device);

  ACTS_DEBUG("inputTensor: " << inputTensor);
  ACTS_DEBUG("edgeList: " << edgeList);
  ACTS_DEBUG("edgeFeatures: " << edgeFeatures);

  //std::cout << "clone edgeList" << std::endl; auto edgeList2 = edgeList.clone();
  //std::cout << "clone edgeFeatures" << std::endl;
  //auto edgeFeatures2 = edgeFeatures.clone();
  //std::cout << "clone input" << std::endl;
  //auto inputTensor2 = inputTensor.clone();
  return {inputTensorCuda, edgeListCuda, edgeFeaturesCuda};
}


}
