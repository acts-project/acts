// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

namespace podio {
class Frame;
}

namespace ActsExamples {

class EDM4hepWriter;
struct AlgorithmContext;

class EDM4hepOutputConverter {
 public:
  EDM4hepOutputConverter(std::unique_ptr<const Acts::Logger> logger)
      : m_logger(std::move(logger)) {}

  virtual void convert(const AlgorithmContext& ctx,
                       podio::Frame& output) const = 0;

  virtual void initialize(EDM4hepWriter& writer) = 0;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples