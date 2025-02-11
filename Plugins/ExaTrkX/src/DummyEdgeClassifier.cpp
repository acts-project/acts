// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/DummyEdgeClassifier.hpp"

#include <torch/torch.h>

namespace Acts {

DummyEdgeClassifier::DummyEdgeClassifier(const Config& cfg,
                                         std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {}

DummyEdgeClassifier::~DummyEdgeClassifier() {}

std::tuple<std::any, std::any, std::any, std::any>
DummyEdgeClassifier::operator()(std::any inNodeFeatures, std::any inEdgeIndex,
                                std::any inEdgeFeatures,
                                const ExecutionContext& /*execContext*/) {
  if (m_cfg.keepAll) {
    auto scores = torch::ones(std::any_cast<at::Tensor>(inEdgeIndex).size(1))
                      .to(torch::kFloat32);

    return {std::move(inNodeFeatures), std::move(inEdgeFeatures),
            std::move(inEdgeFeatures), std::move(scores)};
  } else {
    auto noEdges =
        torch::empty({2, 0}, at::TensorOptions{}.dtype(torch::kLong));
    at::Tensor scores =
        torch::empty({0}, at::TensorOptions{}.dtype(torch::kFloat32));

    return {std::move(inNodeFeatures), std::move(noEdges),
            std::move(inEdgeFeatures), std::move(scores)};
  }
}

}  // namespace Acts
