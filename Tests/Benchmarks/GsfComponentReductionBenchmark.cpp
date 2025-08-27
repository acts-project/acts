// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GsfMixtureReduction.hpp"

#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"

using namespace Acts;
using namespace Acts::Test;

constexpr std::size_t nComponents = 72; // 12 max components GSF, 6 material

template<typename Fun>
void test(const std::vector<std::vector<GsfComponent>> &inputData, const Acts::Surface &surface, Fun &&fun) {
  auto num_runs = 1000;
  auto res = microBenchmark([&](const std::vector<GsfComponent> &input) {
    // Avoid reallocation cost every iteration
    static std::vector<GsfComponent> inputCopy(nComponents);
    inputCopy = input;
    fun(inputCopy, 12, surface);
    assumeRead(inputCopy);
  }, inputData, num_runs);

  std::cout << res << std::endl;

}



int main() {
  std::vector<std::vector<GsfComponent>> data(100);

  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Transform3::Identity());

  for(auto &sample : data) {
    sample.reserve(nComponents);
    for(auto i=0ul; i<nComponents; ++i) {
      double weight = 1.0/nComponents;
      auto cov = Acts::BoundMatrix::Identity();
      auto pars = Acts::BoundVector::Random();
      sample.push_back({weight, pars, cov});
    }
  }

  std::cout << "reduceMixtureLargestWeights" << std::endl;
  test(data, *surface, reduceMixtureLargestWeights);
  std::cout << std::endl;
  
  std::cout << "reduceMixtureWithKLDistance" << std::endl;
  test(data, *surface, reduceMixtureWithKLDistance);
  std::cout << std::endl;
}


