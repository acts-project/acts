// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/Detail/AlgebraConversion.hpp"

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

// Acts Examples include(s)
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Traccc/Util/MapUtil.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

// Traccc include(s)
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_state.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <variant>
#include <vector>

namespace ActsExamples::Traccc::Common::Conversion {

// Custom hash and equals functions
// as some are not defined by std::hash and std::equal_to

struct TracccMeasurementHash {
  std::size_t operator()(const traccc::measurement& s) const noexcept {
    return s.measurement_id;
  }
};

struct ActsMeasurementHash {
  std::size_t operator()(
      const ActsExamples::MeasurementContainer::ConstVariableProxy& s)
      const noexcept {
    return static_cast<std::size_t>(
        s.sourceLink().get<ActsExamples::IndexSourceLink>().index());
  }
};

struct ActsMeasurementEquals {
  bool operator()(
      const ActsExamples::MeasurementContainer::ConstVariableProxy& m1,
      const ActsExamples::MeasurementContainer::ConstVariableProxy& m2) const {
    auto lhsIdx = m1.sourceLink().get<ActsExamples::IndexSourceLink>().index();
    auto rhsIdx = m2.sourceLink().get<ActsExamples::IndexSourceLink>().index();
    return lhsIdx == rhsIdx;
  }
};

/// @brief Converts a traccc bound index to an Acts bound index.
/// @param tracccBoundIndex the traccc bound index.
/// @returns an Acts bound index.
Acts::BoundIndices boundIndex(const traccc::bound_indices tracccBoundIndex);

/// @brief Gets the the local position of the measurement.
/// @param measurement the Acts measurement.
/// @returns A two-dimensional vector containing the local position.
/// The first item in the vector is local position on axis 0 and
/// I.e., [local position (axis 0), local position (axis 1)].
/// @note if the dimension is less than 2 then the remaining values are set to 0.
template <std::size_t dim>
inline Acts::ActsVector<2> getLocal(
    const ActsExamples::MeasurementContainer::ConstFixedProxy<dim>&
        measurement) {
  traccc::scalar loc0 = 0;
  traccc::scalar loc1 = 0;
  if constexpr (dim > Acts::BoundIndices::eBoundLoc0) {
    loc0 = measurement.parameters()(Acts::BoundIndices::eBoundLoc0);
  }
  if constexpr (dim > Acts::BoundIndices::eBoundLoc1) {
    loc1 = measurement.parameters()(Acts::BoundIndices::eBoundLoc1);
  }
  return Acts::ActsVector<2>(loc0, loc1);
}

/// @brief Get the the local position of the measurement.
/// @param measurement the Acts bound variant measurement.
/// @return A two-dimensional vector containing the local position.
/// I.e., [local position (axis 0), local position (axis 1)].
/// @note if the dimension is less than 2 then the remaining values are set to 0.
inline Acts::ActsVector<2> getLocal(
    const ActsExamples::MeasurementContainer::ConstVariableProxy& measurement) {
  traccc::scalar loc0 = 0;
  traccc::scalar loc1 = 0;
  if (measurement.size() > Acts::BoundIndices::eBoundLoc0) {
    loc0 = measurement.parameters()(Acts::BoundIndices::eBoundLoc0);
  }
  if (measurement.size() > Acts::BoundIndices::eBoundLoc1) {
    loc1 = measurement.parameters()(Acts::BoundIndices::eBoundLoc1);
  }
  return Acts::ActsVector<2>(loc0, loc1);
}

/// @brief Get the the variance of the measurement.
/// @param measurement the Acts measurement.
/// @return A two-dimensional vector containing the variance.
/// I.e., [variance (axis 0), variance (axis 1)].
/// @note if the dimension is less than 2 then the remaining values are set to 0.
template <std::size_t dim>
inline Acts::ActsVector<2> getVariance(
    const ActsExamples::MeasurementContainer::ConstFixedProxy<dim>&
        measurement) {
  traccc::scalar var0 = 0;
  traccc::scalar var1 = 0;
  if constexpr (dim >= Acts::BoundIndices::eBoundLoc0) {
    var0 = measurement.covariance()(Acts::BoundIndices::eBoundLoc0,
                                    Acts::BoundIndices::eBoundLoc0);
  }
  if constexpr (dim > Acts::BoundIndices::eBoundLoc1) {
    var1 = measurement.covariance()(Acts::BoundIndices::eBoundLoc1,
                                    Acts::BoundIndices::eBoundLoc1);
  }
  return Acts::ActsVector<2>(var0, var1);
}

/// @brief Get the the variance of the measurement.
/// @param measurement the Acts bound variant measurement.
/// @return A two-dimensional vector containing the variance.
/// I.e., [variance (axis 0), variance (axis 1)].
/// @note if the dimension is less than 2 then the remaining values are set to 0.
inline Acts::ActsVector<2> getVariance(
    const ActsExamples::MeasurementContainer::ConstVariableProxy& measurement) {
  traccc::scalar var0 = 0;
  traccc::scalar var1 = 0;
  if (measurement.size() >= Acts::BoundIndices::eBoundLoc0) {
    var0 = measurement.covariance()(Acts::BoundIndices::eBoundLoc0,
                                    Acts::BoundIndices::eBoundLoc0);
  }
  if (measurement.size() > Acts::BoundIndices::eBoundLoc1) {
    var1 = measurement.covariance()(Acts::BoundIndices::eBoundLoc1,
                                    Acts::BoundIndices::eBoundLoc1);
  }
  return Acts::ActsVector<2>(var0, var1);
}

/// @brief Get the geometry ID from the measurement through its source link.
/// @note The sourcelink is assumed to be of type IndexSourceLink.
inline Acts::GeometryIdentifier getGeometryID(
    const ActsExamples::MeasurementContainer::ConstVariableProxy& measurement) {
  return measurement.sourceLink()
      .template get<ActsExamples::IndexSourceLink>()
      .geometryId();
}

/// @brief Converts traccc measurements to acts measurements.
/// @param detector The detray detector,
/// @param measurements The traccc measurements,
/// @return A vector of Acts bound variant measurements.
/// @note The type IndexSourceLink is used for the measurements' source links.
template <typename detector_t, std::forward_iterator iterator_t>
inline auto convertMeasurements(const detector_t& detector,
                                iterator_t measurements_begin,
                                iterator_t measurements_end,
                                MeasurementContainer& measurementContainer) {
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

  return Util::makeConversionOneToOne(
      measurements_begin, measurements_end, measurementContainer.cbegin(),
      measurementContainer.cend(), TracccMeasurementHash{},
      std::equal_to<traccc::measurement>{});
}

}  // namespace ActsExamples::Traccc::Common::Conversion
