// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioMeasurementInputConverter.hpp"

#include "ActsPodioEdm/MeasurementCollection.h"

#include <podio/Frame.h>

namespace ActsExamples {

PodioMeasurementInputConverter::PodioMeasurementInputConverter(
    const Config& config, Acts::Logging::Level level)
    : PodioInputConverter{"PodioMeasurementInputConverter", level,
                          config.inputFrame},
      m_cfg(config) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument(
        "PodioMeasurementInputConverter: "
        "inputMeasurements is not set");
  }
}

ProcessCode PodioMeasurementInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  const auto& measurements =
      frame.get<ActsPodioEdm::MeasurementCollection>(m_cfg.inputMeasurements);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
