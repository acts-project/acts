// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>

namespace ActsExamples {

/// A simple algorithm that just prints hello world.
class HelloLoggerAlgorithm : public IAlgorithm {
 public:
  explicit HelloLoggerAlgorithm(
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  // Log a few messages.
  ProcessCode execute(const AlgorithmContext& ctx) const override;
};

}  // namespace ActsExamples
