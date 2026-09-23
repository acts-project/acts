// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Gnn/OnnxEdgeClassifier.hpp"
#include "ActsTests/CommonHelpers/DataDirectory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <vector>

#ifdef ACTS_GNN_WITH_CUDA
#include <cuda_runtime_api.h>
#endif

using namespace ActsPlugins;

namespace ActsTests {

namespace {

const ExecutionContext execContextCpu{Device::Cpu(), {}};
#ifdef ACTS_GNN_WITH_CUDA
const ExecutionContext execContextCuda{Device::Cuda(0), cudaStreamLegacy};
#endif

// The classifier runs echo_edge_classifier.onnx (see
// make_echo_edge_classifier.py), which returns the node features it is given,
// flattened row by row, as its edge scores. Undoing the sigmoid the classifier
// applies to them gives back the feature matrix that reached the model.

// Two nodes with four features each, e.g. (r, phi, z, t). All values differ,
// so any swapped column or misapplied scale shows up.
constexpr std::size_t nNodes = 2;
constexpr std::size_t nFeatures = 4;
const std::vector<float> nodes = {
    0.5f, -1.5f, 2.5f,  0.25f,  // node 0
    1.0f, 0.75f, -3.0f, 1.75f   // node 1
};

/// Run the echo classifier with @p selectedFeatures and @p featureScales and
/// return the (row-major) node feature matrix the model got to see
std::vector<float> modelInputs(const std::vector<int> &selectedFeatures,
                               const std::vector<float> &featureScales,
                               const ExecutionContext &execContext) {
  OnnxEdgeClassifier::Config cfg;
  cfg.modelPath = getDataPath("echo_edge_classifier.onnx");
  cfg.selectedFeatures = selectedFeatures;
  cfg.featureScales = featureScales;
  // Keep every edge, the sigmoid is always above
  cfg.cut = 0.f;
  cfg.device = execContext.device;
  OnnxEdgeClassifier classifier(
      cfg,
      Acts::getDefaultLogger("OnnxEdgeClassifier", Acts::Logging::WARNING));

  auto nodeFeatures =
      Tensor<float>::Create({nNodes, nFeatures}, execContextCpu);
  std::ranges::copy(nodes, nodeFeatures.data());

  // One edge per value the model reads. The edge index itself only sets the
  // number of scores the model returns.
  const std::size_t nModelFeatures =
      selectedFeatures.empty() ? nFeatures : selectedFeatures.size();
  const std::size_t nEdges = nNodes * nModelFeatures;
  auto edgeIndex = Tensor<std::int64_t>::Create({2, nEdges}, execContextCpu);
  std::fill_n(edgeIndex.data(), edgeIndex.size(), 0);

  PipelineTensors tensors{
      nodeFeatures.clone(execContext), edgeIndex.clone(execContext), {}, {}};
  auto out = classifier(std::move(tensors), execContext);

  // The stage passes the full, unscaled node features on to the next one
  auto outNodes = out.nodeFeatures.clone(execContextCpu);
  BOOST_CHECK_EQUAL_COLLECTIONS(outNodes.data(),
                                outNodes.data() + outNodes.size(),
                                nodes.begin(), nodes.end());

  BOOST_REQUIRE(out.edgeScores.has_value());
  auto scores = out.edgeScores->clone(execContextCpu);
  BOOST_REQUIRE_EQUAL(scores.size(), nEdges);
  std::vector<float> inputs(nEdges);
  std::transform(scores.data(), scores.data() + nEdges, inputs.begin(),
                 [](float s) { return std::log(s / (1.f - s)); });
  return inputs;
}

void checkInputs(const std::vector<float> &actual,
                 const std::vector<float> &expected) {
  BOOST_REQUIRE_EQUAL(actual.size(), expected.size());
  for (std::size_t i = 0; i < actual.size(); ++i) {
    BOOST_TEST_CONTEXT("model input " << i) {
      BOOST_CHECK_SMALL(actual[i] - expected[i], 1e-4f);
    }
  }
}

void testFeatureSelection(const ExecutionContext &ctx) {
  // Without a selection the model gets all features
  checkInputs(modelInputs({}, {}, ctx), nodes);

  // The features reach the model in the configured order: (z, r), not (r, z)
  // as they are ordered in the node features
  checkInputs(modelInputs({2, 0}, {}, ctx), {2.5f, 0.5f, -3.f, 1.f});

  // Each scale divides its own feature: z / 10 and r / 2
  checkInputs(modelInputs({2, 0}, {10.f, 2.f}, ctx),
              {0.25f, 0.25f, -0.3f, 0.5f});

  // Without a selection the scales apply to the features in order
  checkInputs(modelInputs({}, {1.f, 2.f, 4.f, 0.5f}, ctx),
              {0.5f, -0.75f, 0.625f, 0.5f, 1.f, 0.375f, -0.75f, 3.5f});

  // A feature can be selected more than once: phi twice, with a different
  // scale each time, around t
  checkInputs(modelInputs({1, 3, 1}, {1.f, 1.f, 3.f}, ctx),
              {-1.5f, 0.25f, -0.5f, 0.75f, 1.75f, 0.25f});

  // A selected feature has to exist
  BOOST_CHECK_THROW(modelInputs({0, 4}, {}, ctx), std::runtime_error);
  BOOST_CHECK_THROW(modelInputs({-1}, {}, ctx), std::runtime_error);

  // There has to be one scale per model input feature
  BOOST_CHECK_THROW(modelInputs({2, 0}, {1.f}, ctx), std::runtime_error);
  BOOST_CHECK_THROW(modelInputs({}, {1.f}, ctx), std::runtime_error);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GnnOnnxEdgeClassifierSuite)

BOOST_AUTO_TEST_CASE(test_feature_selection_cpu) {
  testFeatureSelection(execContextCpu);
}

#ifdef ACTS_GNN_WITH_CUDA
BOOST_AUTO_TEST_CASE(test_feature_selection_cuda) {
  testFeatureSelection(execContextCuda);
}
#endif

BOOST_AUTO_TEST_CASE(test_missing_model_is_named) {
  OnnxEdgeClassifier::Config cfg;
  cfg.modelPath = "this-model-does-not-exist.onnx";
  cfg.device = Device::Cpu();
  BOOST_CHECK_EXCEPTION(
      OnnxEdgeClassifier(cfg, Acts::getDefaultLogger("OnnxEdgeClassifier",
                                                     Acts::Logging::WARNING)),
      std::runtime_error, [&cfg](const std::runtime_error &e) {
        return std::string_view{e.what()}.find(cfg.modelPath) !=
               std::string_view::npos;
      });
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
