// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>

namespace ActsExamples {

/// A simple algorithm that just prints hello world.
class HelloLoggerAlgorithm : public ActsExamples::IAlgorithm {
 public:
  HelloLoggerAlgorithm(Acts::Logging::Level level);

  // Log a few messages.
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const override;
};

}  // namespace ActsExamples
