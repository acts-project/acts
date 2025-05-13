// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/Tensor.hpp>

template <typename T>
void fillTensor(Acts::Tensor<T>& tensor, const std::vector<T>& data,
                Acts::ExecutionContext execContext) {
  BOOST_CHECK(tensor.size() == data.size());

  if (execContext.device.type == Acts::Device::Type::eCPU) {
    std::copy(data.begin(), data.end(), tensor.data());
  } else {
    ACTS_CUDA_CHECK(cudaMemcpyAsync(tensor.data(), data.data(), tensor.nbytes(),
                                    cudaMemcpyHostToDevice,
                                    *execContext.stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(*execContext.stream));
  }
}

template <typename T>
std::vector<T> copyToHost(const Acts::Tensor<T>& tensor,
                          Acts::ExecutionContext execContext) {
  std::vector<T> data(tensor.size());
  if (execContext.device.type == Acts::Device::Type::eCPU) {
    std::copy(tensor.data(), tensor.data() + tensor.size(), data.begin());
  } else {
    ACTS_CUDA_CHECK(cudaMemcpyAsync(data.data(), tensor.data(), tensor.nbytes(),
                                    cudaMemcpyDeviceToHost,
                                    *execContext.stream));
    ACTS_CUDA_CHECK(cudaStreamSynchronize(*execContext.stream));
  }
  return data;
}

void testSigmoid(std::vector<float> input, Acts::ExecutionContext execContext) {
  auto tensor = Acts::Tensor<float>::Create({input.size(), 1}, execContext);
  fillTensor(tensor, input, execContext);

  Acts::sigmoid(tensor, execContext);

  std::vector<float> expected(input.size());
  std::transform(input.begin(), input.end(), expected.begin(),
                 [](float x) { return 1.f / (1.f + std::exp(-x)); });

  auto result = copyToHost(tensor, execContext);
  BOOST_CHECK(result.size() == expected.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    BOOST_CHECK_CLOSE(result[i], expected[i], 1e-4);
  }
}

void testEdgeSelection(const std::vector<float>& scores,
                       const std::vector<std::int64_t>& edgeIndex,
                       const std::vector<std::int64_t>& edgeIndexExpected,
                       Acts::ExecutionContext execContext) {
  auto scoreTensor =
      Acts::Tensor<float>::Create({scores.size(), 1}, execContext);
  auto edgeTensor = Acts::Tensor<std::int64_t>::Create(
      {2, edgeIndex.size() / 2}, execContext);

  fillTensor(scoreTensor, scores, execContext);
  fillTensor(edgeTensor, edgeIndex, execContext);

  auto [selectedScores, selectedEdges] =
      Acts::applyScoreCut(scoreTensor, edgeTensor, 0.5f, execContext);

  auto selectedScoresHost = copyToHost(selectedScores, execContext);
  auto selectedEdgesHost = copyToHost(selectedEdges, execContext);

  BOOST_CHECK(selectedScoresHost.size() == 2);

  BOOST_CHECK(selectedEdgesHost.size() == edgeIndexExpected.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      selectedEdgesHost.begin(), selectedEdgesHost.end(),
      edgeIndexExpected.begin(), edgeIndexExpected.end());
}

const Acts::ExecutionContext execContextCpu{Acts::Device::Cpu(), {}};

BOOST_AUTO_TEST_CASE(tensor_create_move_cpu) {
  auto tensor = Acts::Tensor<float>::Create({10, 1}, execContextCpu);

  BOOST_CHECK(tensor.shape()[1] == 1);
  BOOST_CHECK(tensor.shape()[0] == 10);

  auto tensor2 = std::move(tensor);
  BOOST_CHECK(tensor2.shape()[1] == 1);
  BOOST_CHECK(tensor2.shape()[0] == 10);
  BOOST_CHECK(tensor2.data() != nullptr);
  BOOST_CHECK(tensor.data() == nullptr);
}

BOOST_AUTO_TEST_CASE(tensor_sigmoid_cpu) {
  testSigmoid({-2.f, -1.f, 0.f, 1.f, 2.f}, execContextCpu);
}

const std::vector<float> scores = {0.1f, 0.4f, 0.6f, 0.9f};
const std::vector<std::int64_t> edgeIndex = {0, 1, 2, 3, 4, 5, 6, 7};
const std::vector<std::int64_t> edgeIndexExpected = {2, 3, 6, 7};

BOOST_AUTO_TEST_CASE(tensor_edge_selection_cpu) {
  testEdgeSelection(scores, edgeIndex, edgeIndexExpected, execContextCpu);
}

#ifdef ACTS_EXATRKX_WITH_CUDA

#include <cuda_runtime_api.h>

const Acts::ExecutionContext execContextCuda{Acts::Device::Cuda(0),
                                             cudaStreamLegacy};

BOOST_AUTO_TEST_CASE(tensor_create_move_cuda) {
  auto tensor = Acts::Tensor<float>::Create({10, 1}, execContextCuda);

  BOOST_CHECK(tensor.shape()[1] == 1);
  BOOST_CHECK(tensor.shape()[0] == 10);

  auto tensor2 = std::move(tensor);
  BOOST_CHECK(tensor2.shape()[1] == 1);
  BOOST_CHECK(tensor2.shape()[0] == 10);
  BOOST_CHECK(tensor2.data() != nullptr);
  BOOST_CHECK(tensor.data() == nullptr);
}

BOOST_AUTO_TEST_CASE(tensor_sigmoid_cuda) {
  testSigmoid({-2.f, -1.f, 0.f, 1.f, 2.f}, execContextCuda);
}

BOOST_AUTO_TEST_CASE(tensor_edge_selection_cuda) {
  testEdgeSelection(scores, edgeIndex, edgeIndexExpected, execContextCuda);
}

#endif
