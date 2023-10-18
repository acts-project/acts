// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXPipeline.hpp"
#include "Acts/Plugins/ExaTrkX/detail/CantorEdge.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class TorchTruthGraphMetricsHook : public ExaTrkXHook {
  std::unique_ptr<const Logger> m_logger;
  std::vector<detail::CantorEdge<int64_t>> m_truthGraphCantor;

  const Logger &logger() const { return *m_logger; }

 public:
  TorchTruthGraphMetricsHook(const std::vector<int64_t> &truthGraph,
                             std::unique_ptr<const Acts::Logger> l);
  ~TorchTruthGraphMetricsHook() override {}

  void operator()(const std::any &, const std::any &edges) const override;
};

}  // namespace Acts
