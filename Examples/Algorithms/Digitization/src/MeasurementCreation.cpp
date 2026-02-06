// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MeasurementCreation.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <stdexcept>
#include <string>

ActsExamples::VariableBoundMeasurementProxy ActsExamples::createMeasurement(
    MeasurementContainer& container, Acts::GeometryIdentifier geometryId,
    const DigitizedParameters& dParams) {
  if (dParams.indices.size() > 4u) {
    std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                           std::to_string(dParams.indices.size());
    throw std::runtime_error(errorMsg.c_str());
  }

  return Acts::visit_measurement(
      dParams.indices.size(), [&](auto dim) -> VariableBoundMeasurementProxy {
        auto [indices, par, cov] = measurementConstituents<dim>(dParams);
        return VariableBoundMeasurementProxy{
            container.emplaceMeasurement<dim>(geometryId, indices, par, cov)};
      });
}

Acts::Vector3 ActsExamples::measurementGlobalPosition(
    const DigitizedParameters& digitizedParameters,
    const Acts::Surface& surface, const Acts::GeometryContext& gctx) {
  // we need a regular surface to perform the local to global transformation
  // without direction input. the direction could be obtained from truth
  // information but we can leave it out here.
  const Acts::RegularSurface* regularSurface =
      dynamic_cast<const Acts::RegularSurface*>(&surface);
  if (regularSurface == nullptr) {
    throw std::invalid_argument(
        "Cannot make global measurement position for a non-regular surface");
  }

  Acts::Vector2 locPos = regularSurface->bounds().center();
  for (auto i = 0ul; i < digitizedParameters.indices.size(); ++i) {
    auto idx = digitizedParameters.indices.at(i);
    if (idx == Acts::eBoundLoc0 || idx == Acts::eBoundLoc1) {
      locPos[idx] = digitizedParameters.values.at(i);
    }
  }

  return regularSurface->localToGlobal(gctx, locPos);
}
