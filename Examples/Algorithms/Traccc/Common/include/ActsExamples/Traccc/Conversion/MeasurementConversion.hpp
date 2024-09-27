// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Traccc/Detail/AlgebraConversion.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <tuple>
#include <variant>
#include <vector>

#include "detray/core/detector.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_state.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Converts a traccc bound index to an Acts bound index.
/// @param tracccBoundIndex the traccc bound index.
/// @returns an Acts bound index.
Acts::BoundIndices boundIndex(const traccc::bound_indices tracccBoundIndex);

/// @brief Converts traccc measurements to acts measurements.
/// @param detector The detray detector,
/// @param measurements The traccc measurements,
/// @return A vector of Acts bound variant measurements.
/// @note The type IndexSourceLink is used for the measurements' source links.
template <typename detector_t, std::forward_iterator iterator_t>
inline MeasurementContainer convertTracccToActsMeasurements(
    const detector_t& detector, iterator_t measurements_begin,
    iterator_t measurements_end) {
  MeasurementContainer measurementContainer;

  for (iterator_t i = measurements_begin; i != measurements_end; ++i) {
    Acts::GeometryIdentifier moduleGeoId(
        detector.surface((*i).surface_link).source);
    Index measurementIdx = static_cast<Index>(measurementContainer.size());
    IndexSourceLink idxSourceLink{moduleGeoId, measurementIdx};

    Eigen::Matrix<MeasurementContainer::VariableProxy::SubspaceIndex,
                  Eigen::Dynamic, 1>
        indices(2);

    for (unsigned int j = 0; j < 2; j++) {
      indices[j] =
          boundIndex(traccc::bound_indices((*i).subs.get_indices()[j]));
    }

    measurementContainer.emplaceMeasurement<2>(
        Acts::SourceLink{idxSourceLink}, indices,
        Acts::TracccPlugin::detail::toActsVector<2>((*i).local),
        Eigen::DiagonalMatrix<Acts::ActsScalar, 2>(
            Acts::TracccPlugin::detail::toActsVector<2>((*i).variance))
            .toDenseMatrix());
  }

  return measurementContainer;
}

}  // namespace ActsExamples::Traccc::Common::Conversion
