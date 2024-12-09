// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <vector>

#include "traccc/edm/cell.hpp"
#include "traccc/edm/measurement.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

namespace {
template <typename C1, typename C2, typename M1GeoIdProj, typename M2GeoIdProj,
          typename M1LocalProj, typename M2LocalProj>
std::map<std::size_t, std::size_t> makeGenericMeasurementIndexMap(
    const C1& fromMeasurements, const C2& toMeasurements,
    M1GeoIdProj&& fromMeasurementGeoIdProjector,
    M2GeoIdProj&& toMeasurementGeoIdProjector,
    M1LocalProj&& fromMeasurementLocalProjector,
    M2LocalProj&& toMeasurementLocalProjector, Acts::ActsScalar tol = 0.001) {
  std::map<std::size_t, std::size_t> outputMap;

  std::map<std::size_t, std::vector<std::size_t>> geometryMap;

  for (auto i = toMeasurements.cbegin(); i != toMeasurements.cend(); ++i) {
    std::size_t key = toMeasurementGeoIdProjector(*i);

    geometryMap[key].push_back(std::distance(toMeasurements.cbegin(), i));
  }

  auto toProjectX = [&toMeasurementLocalProjector,
                     &toMeasurements](const std::size_t& i) {
    return toMeasurementLocalProjector(toMeasurements.at(i)).first;
  };

  for (auto& i : geometryMap) {
    std::ranges::sort(i.second, std::less{}, toProjectX);
  }

  for (auto i = fromMeasurements.cbegin(); i != fromMeasurements.cend(); ++i) {
    std::size_t key = fromMeasurementGeoIdProjector(*i);

    const std::vector<std::size_t>& v = geometryMap[key];

    auto targetX = fromMeasurementLocalProjector(*i).first;

    auto lo =
        std::ranges::lower_bound(v, targetX - tol, std::less{}, toProjectX);
    auto hi =
        std::ranges::upper_bound(v, targetX + tol, std::less{}, toProjectX);

    bool found = false;

    for (auto j = lo; j != hi && !found; ++j) {
      decltype(auto) c = toMeasurements.at(*j);

      if (std::abs(toMeasurementLocalProjector(c).first -
                   fromMeasurementLocalProjector(*i).first) <= tol &&
          std::abs(toMeasurementLocalProjector(c).second -
                   fromMeasurementLocalProjector(*i).second) <= tol) {
        outputMap[std::distance(fromMeasurements.cbegin(), i)] = *j;
        found = true;
      }
    }
  }

  return outputMap;
}

template <typename detector_t>
std::map<std::size_t, std::size_t> makeTracccToActsMeasurementIndexMap(
    const std::pmr::vector<traccc::measurement>& tracccMeasurements,
    const ActsExamples::MeasurementContainer& actsMeasurements,
    const detector_t& detector) {
  return makeGenericMeasurementIndexMap(
      tracccMeasurements, actsMeasurements,
      [&detector](const traccc::measurement& i) {
        return detector.surface(i.surface_link).source;
      },
      [](const ActsExamples::MeasurementContainer::ConstVariableProxy& i) {
        return i.sourceLink()
            .template get<ActsExamples::IndexSourceLink>()
            .geometryId()
            .value();
      },
      [](const traccc::measurement& i) {
        return std::make_pair(i.local[0], i.local[1]);
      },
      [](const ActsExamples::MeasurementContainer::ConstVariableProxy& i) {
        return std::make_pair(i.fullParameters()[0], i.fullParameters()[1]);
      });
}

template <typename detector_t>
std::map<std::size_t, std::size_t> makeActsToTracccMeasurementIndexMap(
    const ActsExamples::MeasurementContainer& actsMeasurements,
    const std::pmr::vector<traccc::measurement>& tracccMeasurements,
    const detector_t& detector) {
  return makeGenericMeasurementIndexMap(
      actsMeasurements, tracccMeasurements,
      [](const ActsExamples::MeasurementContainer::ConstVariableProxy& i) {
        return i.sourceLink()
            .template get<ActsExamples::IndexSourceLink>()
            .geometryId()
            .value();
      },
      [&detector](const traccc::measurement& i) {
        return detector.surface(i.surface_link).source;
      },
      [](const ActsExamples::MeasurementContainer::ConstVariableProxy& i) {
        return std::make_pair(i.fullParameters()[0], i.fullParameters()[1]);
      },
      [](const traccc::measurement& i) {
        return std::make_pair(i.local[0], i.local[1]);
      });
}
}  // namespace

}  // namespace ActsExamples::Traccc::Common::Conversion
