// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MeasurementCreation.hpp"

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <stdexcept>
#include <string>
#include <utility>

ActsExamples::Measurement ActsExamples::createMeasurement(
    const DigitizedParameters& dParams, const IndexSourceLink& isl) {
  Acts::SourceLink sl{isl};
  switch (dParams.indices.size()) {
    case 1u: {
      auto [indices, par, cov] = measurementConstituents<1>(dParams);
      return Acts::Measurement<Acts::BoundIndices, 1>(std::move(sl), indices,
                                                      par, cov);
    }
    case 2u: {
      auto [indices, par, cov] = measurementConstituents<2>(dParams);
      return Acts::Measurement<Acts::BoundIndices, 2>(std::move(sl), indices,
                                                      par, cov);
    };
    case 3u: {
      auto [indices, par, cov] = measurementConstituents<3>(dParams);
      return Acts::Measurement<Acts::BoundIndices, 3>(std::move(sl), indices,
                                                      par, cov);
    };
    case 4u: {
      auto [indices, par, cov] = measurementConstituents<4>(dParams);
      return Acts::Measurement<Acts::BoundIndices, 4>(std::move(sl), indices,
                                                      par, cov);
    };
    default:
      std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                             std::to_string(dParams.indices.size());
      throw std::runtime_error(errorMsg.c_str());
  }
}
