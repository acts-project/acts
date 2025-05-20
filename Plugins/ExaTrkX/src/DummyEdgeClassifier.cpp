// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/DummyEdgeClassifier.hpp"

namespace Acts {

DummyEdgeClassifier::DummyEdgeClassifier(const Config& cfg,
                                         std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {}

DummyEdgeClassifier::~DummyEdgeClassifier() {}

PipelineTensors DummyEdgeClassifier::operator()(
    PipelineTensors tensors, const ExecutionContext& execContext) {
  if (m_cfg.keepAll) {
    auto scores = Acts::Tensor<float>::Create(
        {1, tensors.edgeIndex.shape().at(1)}, {Device::Cpu(), {}});
    std::fill(scores.data(), scores.data() + scores.size(), 1.f);

    tensors.edgeScores = scores.clone(execContext);
  } else {
    tensors.edgeIndex = Acts::Tensor<std::int64_t>::Create({2, 0}, execContext);
    tensors.edgeScores = Acts::Tensor<float>::Create({0, 1}, execContext);
  }
  return tensors;
}

}  // namespace Acts
