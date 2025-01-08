// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

// #include "Acts/Plugins/ExaTrkX/detail/GraphCreatorWrapper.hpp"
#include "Acts/Plugins/ExaTrkX/ModuleMapCuda.hpp"

#include <filesystem>

#include <torch/torch.h>

using namespace Acts;
using namespace std::string_literals;

const static std::vector<std::uint64_t> moduleIds = {
    1152288185909248,   1154487209164800,   1156686232420352,
    19465753858146304,  19467952881401856,  19470151904657408,
    19472350927912960,  38025510135005184,  38027709158260736,
    38029908181516288,  56750193156030464,  75334138688700416,
    145381825471053824, 145384024494309376, 670113304466685952};

const static std::vector<float> etas = {
    1.1035867, 1.4966115, 1.4966115,   0.8026956,   0.45476568,
    1.1420625, 1.1420625, 1.4920081,   1.1011527,   1.4920081,
    1.4920081, 1.1011527, -0.09873319, -0.13823439, 0.01685799};

const static std::vector<float> x = {
    9.97385031,  10.80896891, 11.76712838, 20.34355448, 22.49253009,
    24.70090843, 27.28929414, 36.765937,   39.81633959, 43.23767641,
    52.2216424,  69.01141738, 8.70596993,  8.67200961,  88.41495007};

const static std::vector<float> y = {
    -41.44774342,  -45.03218228,  -49.02648202,  -85.18263063,  -94.22355498,
    -103.77090933, -114.81141496, -155.49001862, -168.77518535, -183.66889409,
    -222.83900877, -297.01520898, -36.1706782,   -36.1622108,   -383.94170063};

const static std::vector<float> z = {
    -265.,   -293.,   -324.,   -606.,   -677.,         -751.,  -837.,    -1154.,
    -1257.5, -1373.5, -1678.5, -2255.5, -223.92666667, -223.1, -2860.715};

auto makeFeatureVector() {
  std::vector<float> f;
  for (auto i = 0ul; i < z.size(); ++i) {
    f.push_back(std::hypot(x[i], y[i]) / 1000.0);
    f.push_back(std::atan2(y[i], x[i]) / 3.1415926);
    f.push_back(z[i] / 1000.);
    f.push_back(etas[i]);
  }
  return f;
}

auto getModuleMapPath() {
  auto &mts = boost::unit_test::framework::master_test_suite();
  BOOST_REQUIRE(mts.argc == 2);
  auto moduleMapPath = mts.argv[1];
  BOOST_REQUIRE(std::filesystem::exists(mts.argv[1] + ".doublets.root"s));
  BOOST_REQUIRE(std::filesystem::exists(mts.argv[1] + ".triplets.root"s));
  return moduleMapPath;
}

#if 0
BOOST_AUTO_TEST_CASE(test_graph_creator_wrapper) {
  Acts::detail::GraphCreatorWrapperCuda wrapper(moduleMapPath, 0, 512);

  auto logger = Acts::getDefaultLogger("test", Acts::Logging::INFO);
  auto features = makeFeatureVector();
  auto [a, b] = wrapper.build(features, moduleIds, *logger);

  std::cout << "Edges:\n" << a << std::endl;
  std::cout << "Edge features:\n" << b << std::endl;
}
#endif

BOOST_AUTO_TEST_CASE(test_cuda_module_map) {
  Acts::ModuleMapCuda::Config cfg;
  cfg.rScale = 1000.f;
  cfg.phiScale = 3.1415926f;
  cfg.zScale = 1000.f;
  cfg.moduleMapPath = getModuleMapPath();

  auto logger = Acts::getDefaultLogger("test", Acts::Logging::VERBOSE);

  Acts::ModuleMapCuda moduleMap(cfg, std::move(logger));

  cudaStream_t stream;
  cudaStreamCreate(&stream);

  Acts::ExecutionContext ctx;
  ctx.device = torch::Device(torch::kCUDA);
  ctx.stream = c10::cuda::getStreamFromExternal(stream, 0);

  auto features = makeFeatureVector();
  auto [nodes, edgeIndex, edgeFeatures] =
      moduleMap(features, moduleIds.size(), moduleIds, ctx);

  ctx.stream->synchronize();
  // cudaStreamDestroy(stream);
  BOOST_CHECK(cudaGetLastError() == cudaSuccess);

  auto nodeTensor = std::any_cast<torch::Tensor>(nodes);
  auto edgeIndexTensor = std::any_cast<torch::Tensor>(edgeIndex);
  auto edgeFeatureTensor = std::any_cast<torch::Tensor>(edgeFeatures);
}
