// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLinkConcept.hpp"

#include <variant>

namespace Acts {

class Surface;

namespace MeasurementHelpers {

/// @brief Extract surface from a type erased measurement object
/// @tparam T The FittableMeasurement type
/// @return const pointer to the extracted surface
template <typename T>
const Surface* getSurface(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return &meas.referenceSurface(); },
                    fittable_measurement);
}

template <typename T>
size_t getSize(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return meas.size(); },
                    fittable_measurement);
}
}  // namespace MeasurementHelpers

struct MinimalSourceLink {
  const FittableMeasurement<MinimalSourceLink>* meas{nullptr};

  bool operator==(const MinimalSourceLink& rhs) const {
    return meas == rhs.meas;
  }

  const Surface& referenceSurface() const {
    return *MeasurementHelpers::getSurface(*meas);
  }

  const FittableMeasurement<MinimalSourceLink>& operator*() const {
    return *meas;
  }
};

inline std::ostream& operator<<(std::ostream& os, const MinimalSourceLink& sl) {
  os << "SourceLink(" << sl.meas << ")";
  return os;
}

static_assert(SourceLinkConcept<MinimalSourceLink>,
              "MinimalSourceLink does not fulfill SourceLinkConcept");
}  // namespace Acts
