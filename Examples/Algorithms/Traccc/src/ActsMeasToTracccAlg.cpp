// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/ActsMeasToTracccAlg.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

#include <detray/geometry/identifier.hpp>

namespace ActsExamples {

ActsMeasToTracccAlg::ActsMeasToTracccAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ActsMeasToTracccAlg", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.detector == nullptr) {
    throw std::invalid_argument("ActsMeasToTracccAlg: detector is not set");
  }
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
}

ProcessCode ActsMeasToTracccAlg::execute(const AlgorithmContext& ctx) const {
  const MeasurementContainer& measurements = m_inputMeasurements(ctx);
  const TracccDetector& detector = *m_cfg.detector;

  TracccHostMeasurements tracccMeasurements{m_mr};

  std::size_t nUnmapped = 0;
  std::size_t nNoLocal = 0;
  std::size_t n2D = 0;

  for (std::size_t i = 0; i < measurements.size(); ++i) {
    const auto meas = measurements.getMeasurement(i);

    const auto detrayId = detector.surfaceIds().toDetray(meas.geometryId());
    if (!detrayId.has_value()) {
      ACTS_VERBOSE("No detray surface for " << meas.geometryId());
      ++nUnmapped;
      continue;
    }

    // Collect the local position components, keep the time if present
    std::array<float, 2> position{0.f, 0.f};
    std::array<float, 2> variance{0.f, 0.f};
    std::array<std::uint8_t, 2> subspace{0u, 1u};
    unsigned int dimensions = 0u;
    float time = 0.f;

    const auto indices = meas.subspaceIndexVector();
    for (std::size_t j = 0; j < indices.size(); ++j) {
      const auto index = static_cast<std::uint8_t>(indices[j]);
      if (index == Acts::eBoundLoc0 || index == Acts::eBoundLoc1) {
        position[index] = static_cast<float>(meas.parameters()[j]);
        variance[index] = static_cast<float>(meas.covariance()(j, j));
        subspace[dimensions] = index;
        ++dimensions;
      } else if (index == Acts::eBoundTime) {
        time = static_cast<float>(meas.parameters()[j]);
      }
    }
    if (dimensions == 0u) {
      ++nNoLocal;
      continue;
    }
    if (dimensions == 2u) {
      ++n2D;
    }

    tracccMeasurements.push_back({});
    auto tm = tracccMeasurements.at(tracccMeasurements.size() - 1);
    tm.local_position() = position;
    tm.local_variance() = variance;
    tm.dimensions() = dimensions;
    tm.set_subspace(subspace);
    tm.time() = time;
    tm.diameter() = 0.f;
    tm.identifier() = static_cast<unsigned int>(i);
    tm.surface_link() = detray::geometry::identifier{*detrayId};
    tm.cluster_index() = 0u;
  }

  if (nUnmapped > 0 || nNoLocal > 0) {
    ACTS_WARNING("Skipped "
                 << nUnmapped << " measurements without detray surface and "
                 << nNoLocal << " measurements without local coordinates");
  }
  ACTS_DEBUG("Converted " << tracccMeasurements.size() << " / "
                          << measurements.size() << " measurements, " << n2D
                          << " of them 2D");

  m_outputMeasurements(ctx, std::move(tracccMeasurements));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
