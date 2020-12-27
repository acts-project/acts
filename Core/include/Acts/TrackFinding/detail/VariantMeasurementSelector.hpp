// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <algorithm>
#include <variant>
#include <vector>

namespace Acts {
namespace detail {

/// Sort variant measurements
///
/// @tparam source_link_t Source link type to connect to the detector readout
/// @tparam indices_t Parameter index type, determines the full parameter space
///
/// @param measurements The variant measurements to be sorted
/// @param kIndex Enum index value to identify a parameter
template <typename source_link_t, typename indices_t>
void sortVariantMeasurements(
    std::vector<VariantMeasurement<source_link_t, indices_t>>& measurements,
    indices_t kIndex) {
  std::sort(measurements.begin(), measurements.end(),
            [&](const auto& lhs, const auto& rhs) {
              bool less = false;
              std::visit(
                  [&](const auto& x, const auto& y) {
                    if (x.contains(kIndex) and y.contains(kIndex)) {
                      less = x.parameters()[kIndex] < y.parameters()[kIndex];
                    }
                  },
                  lhs, rhs);
              return less;
            });
}

/// Find the variant measurements close to a provided track parameters
///
/// @tparam source_link_t Source link type to connect to the detector readout
/// @tparam indices_t Parameter index type, determines the full parameter space
///
/// @param measurements The variant measurements to be sorted
/// @param params The bound or free track parameters vector
///
/// @return the measurements which are neighbors of the track parameters
template <typename source_link_t, typename indices_t>
auto findVariantMeasurements(
    std::vector<VariantMeasurement<source_link_t, indices_t>>& measurements,
    const ActsVector<kParametersSize<indices_t>>& params)
    -> std::vector<VariantMeasurement<source_link_t, indices_t>> {
  std::vector<VariantMeasurement<source_link_t, indices_t>> newMeasurements =
      measurements;
  std::vector<VariantMeasurement<source_link_t, indices_t>>
      selectedMeasurements;
  selectedMeasurements.reserve(10);
  for (size_t ipara = 0; ipara < kParametersSize<indices_t>; ++ipara) {
    indices_t kIndex = static_cast<indices_t>(ipara);
    // Sort measurements using their parameters with kIndex
    sortVariantMeasurements<source_link_t, indices_t>(newMeasurements, kIndex);

    // Find the measurement not less than the provided parameter
    auto it = std::find_if(
        newMeasurements.begin(), newMeasurements.end(), [&](const auto& meas) {
          bool found = false;
          std::visit(
              [&](const auto& m) {
                found = m.parameters()[kIndex] >= params[kIndex];
              },
              meas);
          return found;
        });

    // Make a range with the found measurement
    auto itLow =
        it - newMeasurements.begin() < 5 ? newMeasurements.begin() : it - 5;
    auto itUp = newMeasurements.end() - it < 5 ? newMeasurements.end() : it + 5;
    selectedMeasurements.clear();
    std::copy(itLow, itUp, std::back_inserter(selectedMeasurements));
    newMeasurements = selectedMeasurements;
  }

  return newMeasurements;
}

}  // namespace detail
}  // namespace Acts
