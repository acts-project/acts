// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <utility>
#include <variant>

#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

std::pair<SimSpacePointContainer, std::map<std::size_t, std::size_t>>
convertTracccToActsSpacePoints(
    std::vector<traccc::spacepoint,
                std::pmr::polymorphic_allocator<traccc::spacepoint>>&
        spacePoints,
    const std::map<std::size_t, std::size_t>& tracccToActsMeasurementIndexMap,
    const MeasurementContainer& actsMeasurements) {
  SimSpacePointContainer convertedSpacePoints;
  std::map<std::size_t, std::size_t> outMap;

  for (std::size_t i = 0; i < spacePoints.size(); ++i) {
    const traccc::spacepoint& spacePoint = spacePoints.at(i);

    const Acts::Vector3 globalPos(spacePoint.x(), spacePoint.y(),
                                  spacePoint.z());
    const std::optional<Acts::ActsScalar> t = std::nullopt;
    const Acts::ActsScalar varRho = 0;
    const Acts::ActsScalar varZ = 0;
    const std::optional<Acts::ActsScalar> varT = std::nullopt;

    if (auto measurementIndexIt = tracccToActsMeasurementIndexMap.find(
            spacePoint.meas.measurement_id);
        measurementIndexIt != tracccToActsMeasurementIndexMap.cend()) {
      const Acts::SourceLink sourceLink =
          actsMeasurements.at(measurementIndexIt->second).sourceLink();

      boost::container::static_vector<Acts::SourceLink, 2> sourceLinks = {
          sourceLink};

      outMap[i] = convertedSpacePoints.size();
      convertedSpacePoints.emplace_back(globalPos, t, varRho, varZ, varT,
                                        std::move(sourceLinks));
    }
  }

  return std::make_pair(std::move(convertedSpacePoints), std::move(outMap));
}

std::pair<std::vector<traccc::spacepoint,
                      std::pmr::polymorphic_allocator<traccc::spacepoint>>,
          std::map<std::size_t, std::size_t>>
convertActsToTracccSpacePoints(
    const SimSpacePointContainer& actsSpacePoints,
    const std::map<std::size_t, std::size_t>& actsToTracccMeasurementIndexMap,
    const std::vector<traccc::measurement,
                      std::pmr::polymorphic_allocator<traccc::measurement>>&
        tracccMeasurements) {
  std::vector<traccc::spacepoint,
              std::pmr::polymorphic_allocator<traccc::spacepoint>>
      outputContainer;

  std::map<std::size_t, std::size_t> outputMap;

  for (std::size_t i = 0; i < actsSpacePoints.size(); ++i) {
    const SimSpacePoint& spacePoint = actsSpacePoints.at(i);

    auto idx = spacePoint.sourceLinks()[0]
                   .get<ActsExamples::IndexSourceLink>()
                   .index();

    traccc::point3 global{spacePoint.x(), spacePoint.y(), spacePoint.z()};

    outputMap[i] = outputContainer.size();
    outputContainer.emplace_back(
        global, tracccMeasurements.at(actsToTracccMeasurementIndexMap.at(idx)));
  }

  return std::make_pair(std::move(outputContainer), std::move(outputMap));
}

}  // namespace ActsExamples::Traccc::Common::Conversion
