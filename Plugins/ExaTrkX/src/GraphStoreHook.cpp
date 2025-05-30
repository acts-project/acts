// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/GraphStoreHook.hpp"

Acts::GraphStoreHook::GraphStoreHook() {
  m_storedGraph = std::make_unique<Graph>();
}

void Acts::GraphStoreHook::operator()(const PipelineTensors &tensors,
                                      const ExecutionContext &execCtx) const {
  auto edgeIndexTensor =
      tensors.edgeIndex.clone({Device::Cpu(), execCtx.stream});

  // We need to transpose the edges here for the right memory layout
  m_storedGraph->first.reserve(edgeIndexTensor.size());
  for (auto i = 0ul; i < edgeIndexTensor.shape().at(1); ++i) {
    m_storedGraph->first.push_back(*(edgeIndexTensor.data() + i));
    m_storedGraph->first.push_back(
        *(edgeIndexTensor.data() + edgeIndexTensor.shape().at(1) + i));
  }

  if (!tensors.edgeScores.has_value()) {
    return;
  }

  auto scoreTensor = tensors.edgeScores->clone({Device::Cpu(), execCtx.stream});

  m_storedGraph->second =
      std::vector(scoreTensor.data(), scoreTensor.data() + scoreTensor.size());
}
