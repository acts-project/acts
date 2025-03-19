// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

namespace podio {
class Frame;
}

namespace ActsExamples {

class EDM4hepInputConverter : public IAlgorithm {
 public:
  EDM4hepInputConverter(const std::string& name, Acts::Logging::Level level,
                        const std::string& inputFrame);

  ProcessCode execute(const ActsExamples::AlgorithmContext& ctx) const final;

  virtual ProcessCode convert(const AlgorithmContext& ctx,
                              const podio::Frame& frame) const = 0;

 private:
  ReadDataHandle<podio::Frame> m_inputFrame{this, "InputFrame"};
};

}  // namespace ActsExamples
