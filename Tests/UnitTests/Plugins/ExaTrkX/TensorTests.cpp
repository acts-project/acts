// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <Acts/Plugins/ExaTrkX/Tensor.hpp>

#ifdef ACTS_EXATRKX_WITH_CUDA
#include <cuda_runtime_api.h>
#endif

const Acts::ExecutionContext execContextCpu{Acts::Device::Cpu(), {}};

template <typename T>
Acts::Tensor<T> createCpuTensor(const std::vector<T>& data,
                                std::array<std::size_t, 2> shape) {
  auto tensor = Acts::Tensor<T>::Create(shape, execContextCpu);
  std::copy(data.begin(), data.end(), tensor.data());
  return tensor;
}

void testSigmoid(std::vector<float> input, Acts::ExecutionContext execContext) {
  auto tensor = createCpuTensor(input, {input.size(), 1ul});

  auto tensorTarget = tensor.clone(execContext);
  Acts::sigmoid(tensorTarget, execContext.stream);
  auto result = tensorTarget.clone({Acts::Device::Cpu(), execContext.stream});

  std::vector<float> expected(input.size());
  std::transform(input.begin(), input.end(), expected.begin(),
                 [](float x) { return 1.f / (1.f + std::exp(-x)); });

  BOOST_CHECK(result.size() == expected.size());
  for (std::size_t i = 0; i < result.size(); ++i) {
    BOOST_CHECK_CLOSE(result.data()[i], expected[i], 1e-4);
  }
}

void testEdgeSelection(const std::vector<float>& scores,
                       const std::vector<std::int64_t>& edgeIndex,
                       const std::vector<std::int64_t>& edgeIndexExpected,
                       Acts::ExecutionContext execContext) {
  auto scoreTensor = createCpuTensor<float>(scores, {scores.size(), 1ul});
  auto edgeTensor = createCpuTensor(edgeIndex, {2, edgeIndex.size() / 2});

  auto scoreTensorTarget = scoreTensor.clone(execContext);
  auto edgeTensorTarget = edgeTensor.clone(execContext);

  auto [selectedScores, selectedEdges] = Acts::applyScoreCut(
      scoreTensorTarget, edgeTensorTarget, 0.5f, execContext.stream);

  auto selectedScoresHost =
      selectedScores.clone({Acts::Device::Cpu(), execContext.stream});
  auto selectedEdgesHost =
      selectedEdges.clone({Acts::Device::Cpu(), execContext.stream});

  BOOST_CHECK(selectedScoresHost.size() == 2);

  BOOST_CHECK(selectedEdgesHost.size() == edgeIndexExpected.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      selectedEdgesHost.data(),
      selectedEdgesHost.data() + selectedEdgesHost.size(),
      edgeIndexExpected.begin(), edgeIndexExpected.end());
}

void testConstructionAndMove(Acts::ExecutionContext execContext) {
  auto tensor = Acts::Tensor<float>::Create({10, 1}, execContext);

  BOOST_CHECK(tensor.shape()[1] == 1);
  BOOST_CHECK(tensor.shape()[0] == 10);

  auto tensor2 = std::move(tensor);
  BOOST_CHECK(tensor2.shape()[1] == 1);
  BOOST_CHECK(tensor2.shape()[0] == 10);
  BOOST_CHECK(tensor2.data() != nullptr);
  BOOST_CHECK(tensor.data() == nullptr);
}

void testEdgeLimit(Acts::ExecutionContext execContext) {
  // Create a sample edge tensor with 10 edges on CPU
  auto edgeTensor = Acts::Tensor<std::int64_t>::Create({2, 10}, execContextCpu);
  for (std::size_t i = 0; i < 10; ++i) {
    edgeTensor.data()[i] = i;
    edgeTensor.data()[i + 10] = 2 * i;
  }

  // Create a sample edge feature tensor with 3 };
  auto edgeFeatureTensor = Acts::Tensor<float>::Create({10, 3}, execContextCpu);
  for (std::size_t i = 0; i < 10; ++i) {
    edgeFeatureTensor.data()[i * 3] = static_cast<float>(i);  // Feature 1
    edgeFeatureTensor.data()[i * 3 + 1] =
        static_cast<float>(i + 1);  // Feature 2
    edgeFeatureTensor.data()[i * 3 + 2] =
        static_cast<float>(i + 2);  // Feature 3
  }

  // Clone to execContext
  auto edgeTensorTarget = edgeTensor.clone(execContext);
  std::optional<Acts::Tensor<float>> edgeFeatureTensorTarget =
      edgeFeatureTensor.clone(execContext);

  // Apply edge limit
  auto [limitedEdges, limitedEdgeFeatures] = Acts::applyEdgeLimit(
      edgeTensorTarget, edgeFeatureTensorTarget, 5, execContext.stream);

  // Clone results back to CPU
  auto limitedEdgesHost =
      limitedEdges.clone({Acts::Device::Cpu(), execContext.stream});
  auto limitedEdgeFeaturesHost =
      limitedEdgeFeatures->clone({Acts::Device::Cpu(), execContext.stream});

  // Check the size of the limited edges
  BOOST_CHECK(limitedEdgesHost.shape()[1] == 5);
  BOOST_CHECK(limitedEdgesHost.shape()[0] == 2);
  BOOST_CHECK(limitedEdgeFeaturesHost.shape()[0] == 5);
  BOOST_CHECK(limitedEdgeFeaturesHost.shape()[1] == 3);

  // Check the data of the limited edges
  for (std::size_t i = 0; i < 5; ++i) {
    BOOST_CHECK(limitedEdgesHost.data()[i] == edgeTensor.data()[i]);
    BOOST_CHECK(limitedEdgesHost.data()[i + 5] == edgeTensor.data()[i + 10]);

    BOOST_CHECK(limitedEdgeFeaturesHost.data()[i * 3] ==
                edgeFeatureTensor.data()[i * 3]);
    BOOST_CHECK(limitedEdgeFeaturesHost.data()[i * 3 + 1] ==
                edgeFeatureTensor.data()[i * 3 + 1]);
    BOOST_CHECK(limitedEdgeFeaturesHost.data()[i * 3 + 2] ==
                edgeFeatureTensor.data()[i * 3 + 2]);
  }
}

BOOST_AUTO_TEST_CASE(tensor_create_move_cpu) {
  testConstructionAndMove(execContextCpu);
}

BOOST_AUTO_TEST_CASE(test_clone_cpu) {
  std::vector<float> data = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f};
  auto tensor = createCpuTensor(data, {3, 2});
  auto tensorClone = tensor.clone(execContextCpu);

  BOOST_CHECK(tensorClone.shape()[0] == 3);
  BOOST_CHECK(tensorClone.shape()[1] == 2);
  BOOST_CHECK(tensorClone.data() != nullptr);
  BOOST_CHECK(tensorClone.data() != tensor.data());
  BOOST_CHECK(tensorClone.size() == tensor.size());
  BOOST_CHECK(tensorClone.nbytes() == tensor.nbytes());

  BOOST_CHECK_EQUAL_COLLECTIONS(tensorClone.data(),
                                tensorClone.data() + tensorClone.size(),
                                data.begin(), data.end());
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

BOOST_AUTO_TEST_CASE(tensor_edge_limit_cpu) {
  testEdgeLimit(execContextCpu);
}

#ifdef ACTS_EXATRKX_WITH_CUDA

const Acts::ExecutionContext execContextCuda{Acts::Device::Cuda(0),
                                             cudaStreamLegacy};

BOOST_AUTO_TEST_CASE(tensor_create_move_cuda) {
  testConstructionAndMove(execContextCuda);
}

BOOST_AUTO_TEST_CASE(tensor_clone_roundtrip) {
  std::vector<float> data = {1.f, 2.f, 3.f, 4.f, 5.f, 6.f};
  auto tensorOrigHost = createCpuTensor(data, {3, 2});
  BOOST_CHECK(tensorOrigHost.device().isCpu());

  auto tensorClone = tensorOrigHost.clone(execContextCuda);
  BOOST_CHECK(tensorClone.device().isCuda());

  auto tensorCloneCuda = tensorClone.clone(execContextCuda);
  BOOST_CHECK(tensorCloneCuda.device().isCuda());

  auto tensorCloneHost =
      tensorCloneCuda.clone({Acts::Device::Cpu(), cudaStreamLegacy});
  BOOST_CHECK(tensorCloneHost.device().isCpu());

  BOOST_CHECK(tensorCloneHost.shape()[0] == 3);
  BOOST_CHECK(tensorCloneHost.shape()[1] == 2);
  BOOST_CHECK(tensorCloneHost.data() != nullptr);
  BOOST_CHECK(tensorCloneHost.data() != tensorCloneCuda.data());
  BOOST_CHECK(tensorCloneHost.size() == tensorCloneCuda.size());
  BOOST_CHECK(tensorCloneHost.nbytes() == tensorCloneCuda.nbytes());
  BOOST_CHECK_EQUAL_COLLECTIONS(tensorCloneHost.data(),
                                tensorCloneHost.data() + tensorCloneHost.size(),
                                data.begin(), data.end());
}

BOOST_AUTO_TEST_CASE(tensor_sigmoid_cuda) {
  testSigmoid({-2.f, -1.f, 0.f, 1.f, 2.f}, execContextCuda);
}

BOOST_AUTO_TEST_CASE(tensor_edge_selection_cuda) {
  testEdgeSelection(scores, edgeIndex, edgeIndexExpected, execContextCuda);
}

BOOST_AUTO_TEST_CASE(tensor_edge_limit_cuda) {
  testEdgeLimit(execContextCuda);
}

#endif
