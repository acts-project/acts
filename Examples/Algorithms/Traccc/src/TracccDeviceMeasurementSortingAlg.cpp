// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccDeviceMeasurementSortingAlg.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

TracccDeviceMeasurementSortingAlg::TracccDeviceMeasurementSortingAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TracccDeviceMeasurementSortingAlg", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument(
        "TracccDeviceMeasurementSortingAlg: detector is not set");
  }
  m_sorting = m_cfg.detector->backend().makeMeasurementSorting(
      this->logger().clone("TracccMeasurementSorting"));

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
}

ProcessCode TracccDeviceMeasurementSortingAlg::execute(
    const AlgorithmContext& ctx) const {
  const TracccDeviceMeasurements& unsorted = m_inputMeasurements(ctx);

  // The sorting takes the number of measurements from the size of its input,
  // so the input buffer has to be resizable (the default of the copy tool)
  TracccDeviceMeasurements sorted = (**m_sorting)(unsorted);
  m_cfg.detector->backend().synchronize();

  ACTS_DEBUG("Sorted " << m_sorting->copy().get_size(sorted)
                       << " measurements on the device");

  m_outputMeasurements(ctx, std::move(sorted));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
