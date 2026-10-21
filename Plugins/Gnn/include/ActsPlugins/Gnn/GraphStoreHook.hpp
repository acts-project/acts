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

/// Hook for storing graph data during GNN pipeline execution
class GraphStoreHook : public GnnHook {
 public:
  /// Graph data structure (edges, edge features)
  using Graph = std::pair<std::vector<std::int64_t>, std::vector<float>>;

 private:
  std::unique_ptr<Graph> m_storedGraph;

 public:
  GraphStoreHook();

  /// @brief Store graph data during pipeline execution
  /// @param tensors Pipeline tensor data
  /// @param execCtx Execution context
  void operator()(const PipelineTensors &tensors,
                  const ExecutionContext &execCtx) const override;

  /// @brief Get the stored graph data
  /// @return Reference to the stored graph
  const Graph &storedGraph() const { return *m_storedGraph; }
};

/// @}
}  // namespace ActsPlugins
