// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include <boost/program_options.hpp>

#include "StepperBenchmarkCommons.hpp"

namespace po = boost::program_options;
using namespace Acts;
using namespace ActsTests;

constexpr std::size_t nComponents = 72;  // 12 max components GSF, 6 material

template <typename Fun>
void test(const std::vector<std::vector<GsfComponent>> &inputData,
          const Acts::Surface &surface, Fun &&fun) {
  auto num_runs = 1000;
  auto res = microBenchmark(
      [&](const std::vector<GsfComponent> &input) {
        // Avoid reallocation cost every iteration
        static std::vector<GsfComponent> inputCopy(nComponents);
        inputCopy = input;
        fun(inputCopy, 12, surface);
        assumeRead(inputCopy);
      },
      inputData, num_runs);

  std::cout << res << std::endl;
}

int main(int argc, char *argv[]) {
  bool weight = true;
  bool kl_optimized = true;
  bool kl_naive = true;

  try {
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")(
        "weight", po::value<bool>(&weight)->default_value(true),
        "run weight-based reduction (reduceMixtureLargestWeights)")(
        "kl-optimized", po::value<bool>(&kl_optimized)->default_value(true),
        "run optimized KL distance reduction (reduceMixtureWithKLDistance)")(
        "kl-naive", po::value<bool>(&kl_naive)->default_value(true),
        "run naive KL distance reduction (reduceMixtureWithKLDistanceNaive)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.contains("help")) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  std::vector<std::vector<GsfComponent>> data(100);

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Transform3::Identity());

  for (auto &sample : data) {
    sample.reserve(nComponents);
    for (auto i = 0ul; i < nComponents; ++i) {
      double weight_val = 1.0 / nComponents;
      auto cov = Acts::BoundMatrix::Identity();
      auto pars = Acts::BoundVector::Random();
      sample.push_back({weight_val, pars, cov});
    }
  }

  if (weight) {
    std::cout << "reduceMixtureLargestWeights" << std::endl;
    test(data, *surface, reduceMixtureLargestWeights);
    std::cout << std::endl;
  }

  if (kl_optimized) {
    std::cout << "reduceMixtureWithKLDistance (optimized)" << std::endl;
    test(data, *surface, reduceMixtureWithKLDistance);
    std::cout << std::endl;
  }

  if (kl_naive) {
    std::cout << "reduceMixtureWithKLDistanceNaive (baseline)" << std::endl;
    test(data, *surface, reduceMixtureWithKLDistanceNaive);
    std::cout << std::endl;
  }
}
