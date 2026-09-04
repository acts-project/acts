// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <podio/Frame.h>

namespace ActsExamples {

class PodioInputConverter::Impl {
 public:
  Impl(PodioInputConverter& parent, const std::string& inputFrame)
      : m_inputFrame(&parent, "InputFrame") {
    m_inputFrame.initialize(inputFrame);
  }

  ReadDataHandle<podio::Frame> m_inputFrame;
};

PodioInputConverter::PodioInputConverter(
    const std::string& name, const std::string& inputFrame,
    std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm(name, std::move(logger)),
      m_impl(std::make_unique<Impl>(*this, inputFrame)) {}

ProcessCode PodioInputConverter::execute(const AlgorithmContext& ctx) const {
  const podio::Frame& frame = m_impl->m_inputFrame(ctx);
  return convert(ctx, frame);
}

PodioInputConverter::~PodioInputConverter() = default;

}  // namespace ActsExamples
