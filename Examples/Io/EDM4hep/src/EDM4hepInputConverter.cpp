// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepInputConverter.hpp"

#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepInputConverter::EDM4hepInputConverter(const std::string& name,
                                             Acts::Logging::Level level,
                                             const std::string& inputFrame)
    : IAlgorithm(name, level), m_inputFrame(this, "InputFrame") {
  if (inputFrame.empty()) {
    throw std::invalid_argument("Missing input frame");
  }
  m_inputFrame.initialize(inputFrame);
}

ProcessCode EDM4hepInputConverter::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const podio::Frame& frame = m_inputFrame(ctx);
  return convert(ctx, frame);
}
}  // namespace ActsExamples
