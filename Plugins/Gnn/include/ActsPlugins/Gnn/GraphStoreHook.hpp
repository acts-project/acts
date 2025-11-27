// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/GnnPipeline.hpp"

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

class GraphStoreHook : public GnnHook {
 public:
  using Graph = std::pair<std::vector<std::int64_t>, std::vector<float>>;

 private:
  std::unique_ptr<Graph> m_storedGraph;

 public:
  GraphStoreHook();

  void operator()(const PipelineTensors &tensors,
                  const ExecutionContext &execCtx) const override;

  const Graph &storedGraph() const { return *m_storedGraph; }
};

/// @}
}  // namespace ActsPlugins
