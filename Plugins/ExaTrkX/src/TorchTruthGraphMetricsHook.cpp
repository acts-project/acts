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
auto cantor(T a, T b) {
  return a + (a + b) * (a + b + 1) / 2;
}

auto cantorize(std::vector<int64_t> edgeIndex) {
  // Sort the edges, so we get predictable cantor pairs
  for (auto it = edgeIndex.begin(); it != edgeIndex.end(); it += 2) {
    std::sort(it, it + 2);
  }

  // Use cantor pairing to store truth graph, so we can easily use set
  // operations to compute efficiency and purity
  std::vector<int64_t> cantorEdgeIndex;
  cantorEdgeIndex.reserve(edgeIndex.size() / 2);
  for (auto it = edgeIndex.begin(); it != edgeIndex.end(); it += 2) {
    cantorEdgeIndex.push_back(cantor(*it, *std::next(it)));
  }

  std::sort(cantorEdgeIndex.begin(), cantorEdgeIndex.end());

  return cantorEdgeIndex;
}

}  // namespace

Acts::TorchTruthGraphMetricsHook::TorchTruthGraphMetricsHook(
    const std::vector<int64_t>& truthGraph,
    std::unique_ptr<const Acts::Logger> l)
    : m_logger(std::move(l)) {
  // Compute truth cantor graph
  m_truthGraphCantor = cantorize(truthGraph);

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
  auto edgeTensor = std::any_cast<torch::Tensor>(edges);

  if (edgeTensor.size(0) != 2) {
    throw std::invalid_argument("input tensor must have shape (2,N)");
  }

  // bring to shape so we can put it to vector
  edgeTensor = edgeTensor.t().flatten();

  std::vector<int64_t> edgeIndex(
      edgeTensor.data_ptr<int64_t>(),
      edgeTensor.data_ptr<int64_t>() + edgeTensor.numel());

  auto predGraphCantor = cantorize(edgeIndex);

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
