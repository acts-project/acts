// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchTruthGraphMetricsHook.hpp"

#include <torch/torch.h>

namespace {

template <typename T>
auto cantor(const T& t) {
  return t[0] + (t[0] + t[1]) * (t[0] + t[1] + 1) / 2;
}

}  // namespace

Acts::TorchTruthGraphMetricsHook::TorchTruthGraphMetricsHook(
    const std::vector<int64_t>& truthGraph,
    std::unique_ptr<const Acts::Logger> l)
    : m_logger(std::move(l)) {
  // Sort the edges, so we get predictable cantor pairs
  auto truthGraphSorted = truthGraph;
  for (auto it = truthGraphSorted.begin(); it != truthGraphSorted.end();
       it += 2) {
    std::sort(it, it + 2);
  }

  // Use cantor pairing to store truth graph, so we can easily use set
  // operations to compute efficiency and purity
  m_truthGraphCantor.reserve(truthGraphSorted.size() / 2);
  for (auto it = truthGraphSorted.begin(); it != truthGraphSorted.end();
       it += 2) {
    m_truthGraphCantor.push_back(cantor(it));
  }

  std::sort(m_truthGraphCantor.begin(), m_truthGraphCantor.end());

  // Check if unique (should be!)
  auto new_end =
      std::unique(m_truthGraphCantor.begin(), m_truthGraphCantor.end());
  if (new_end != m_truthGraphCantor.end()) {
    ACTS_WARNING("Truth graph not unique ("
                 << std::distance(new_end, m_truthGraphCantor.end())
                 << " duplicates)");
    m_truthGraphCantor.erase(new_end, m_truthGraphCantor.end());
  }
}

void Acts::TorchTruthGraphMetricsHook::operator()(const std::any&,
                                                  const std::any& edges) const {
  const auto edgeTensor = std::any_cast<torch::Tensor>(edges);

  if (edgeTensor.size(0) != 2) {
    throw std::invalid_argument("input tensor must have shape (2,N)");
  }

  const auto cantorTensor = cantor(std::get<0>(torch::sort(edgeTensor, 0)))
                                .to(torch::kCPU)
                                .to(torch::kInt64);

  assert(cantorTensor.size(0) == edgeTensor.size(1));

  std::vector<int64_t> predGraphCantor(
      cantorTensor.data_ptr<int64_t>(),
      cantorTensor.data_ptr<int64_t>() + cantorTensor.numel());

  // Sort
  std::sort(predGraphCantor.begin(), predGraphCantor.end());

  // Check if unique (should be!)
  auto new_end = std::unique(predGraphCantor.begin(), predGraphCantor.end());
  if (new_end != predGraphCantor.end()) {
    ACTS_WARNING("Input edges not unique ("
                 << std::distance(new_end, predGraphCantor.end())
                 << " duplicates)");
    predGraphCantor.erase(new_end, predGraphCantor.end());
  }

  // Calculate intersection
  std::vector<int64_t> intersection(
      std::max(predGraphCantor.size(), m_truthGraphCantor.size()));

  auto intersection_end =
      std::set_intersection(predGraphCantor.begin(), predGraphCantor.end(),
                            m_truthGraphCantor.begin(),
                            m_truthGraphCantor.end(), intersection.begin());

  float intersected = std::distance(intersection.begin(), intersection_end);

  ACTS_DEBUG("Intersection size " << intersected);
  float eff = intersected / m_truthGraphCantor.size();
  float pur = intersected / predGraphCantor.size();

  ACTS_INFO("Efficiency=" << eff << ", purity=" << pur);
}
