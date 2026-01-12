// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Gnn/GnnPipeline.hpp"
#include "ActsPlugins/Gnn/detail/CantorEdge.hpp"

namespace ActsPlugins {
/// @addtogroup gnn_plugin
/// @{

class TruthGraphMetricsHook : public GnnHook {
  std::unique_ptr<const Acts::Logger> m_logger;
  std::vector<detail::CantorEdge<std::int64_t>> m_truthGraphCantor;

  const Acts::Logger &logger() const { return *m_logger; }

 public:
  TruthGraphMetricsHook(const std::vector<std::int64_t> &truthGraph,
                        std::unique_ptr<const Acts::Logger> l);

  void operator()(const PipelineTensors &tensors,
                  const ExecutionContext &execCtx) const override;
};

/// @}
}  // namespace ActsPlugins
