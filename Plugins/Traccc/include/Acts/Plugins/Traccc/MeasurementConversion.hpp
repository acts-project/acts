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
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

// Traccc include(s)
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/track_state.hpp"

// System include(s)
#include <memory>
#include <variant>

namespace Acts::TracccPlugin {

/// @brief Converts a traccc bound index to an Acts bound index.
/// @param tracccBoundIndex the traccc bound index.
/// @returns an Acts bound index.
inline Acts::BoundIndices boundIndex(
    const traccc::bound_indices tracccBoundIndex) {
  switch (tracccBoundIndex) {
    case traccc::bound_indices::e_bound_loc0:
      return Acts::BoundIndices::eBoundLoc0;
    case traccc::bound_indices::e_bound_loc1:
      return Acts::BoundIndices::eBoundLoc1;
    case traccc::bound_indices::e_bound_phi:
      return Acts::BoundIndices::eBoundPhi;
    case traccc::bound_indices::e_bound_theta:
      return Acts::BoundIndices::eBoundTheta;
    case traccc::bound_indices::e_bound_qoverp:
      return Acts::BoundIndices::eBoundQOverP;
    case traccc::bound_indices::e_bound_time:
      return Acts::BoundIndices::eBoundTime;
    case traccc::bound_indices::e_bound_size:
      return Acts::BoundIndices::eBoundSize;
    default:
      throw std::runtime_error("Could not convert traccc bound index");
  }
}

/// @brief Creates an Acts measurement from a traccc measurement.
/// @tparam the dimension of the Acts measurement (subspace size).
/// @param m the traccc measurement.
/// @param sl the Acts source link to use for the Acts measurement.
/// @returns an Acts measurement with data copied from the traccc measurement
/// and with its source link set to the one provided to the function.
template <std::size_t dim>
inline Acts::Measurement<Acts::BoundIndices, dim> measurement(
    const traccc::measurement& m, const Acts::SourceLink sl) {
  auto params = Detail::newVector<dim>(m.local);
  std::array<Acts::BoundIndices, dim> indices;
  for (unsigned int i = 0; i < dim; i++) {
    indices[i] = boundIndex(traccc::bound_indices(m.subs.get_indices()[i]));
  }
  auto cov = Eigen::DiagonalMatrix<Acts::ActsScalar, static_cast<int>(dim)>(
                 Detail::newVector<dim>(m.variance))
                 .toDenseMatrix();
  return Acts::Measurement<Acts::BoundIndices, dim>(std::move(sl), indices,
                                                    params, cov);
}

/// @brief Creates an Acts bound variant measurement from a traccc measurement.
/// Using recursion, the functions determines the dimension of the traccc
/// measurement which is used for the Acts measurement that the bound variant
/// measurement holds. The dimension must lie between [0; max_dim].
/// @tparam max_dim the largest possible dimension of any measurement type in the variant (default = 4).
/// @param m the traccc measurement.
/// @param sl the Acts source link to use for the Acts measurement.
/// @returns an Acts bound variant measurement with data copied from the traccc measurement
/// and with its source link set to the one provided to the function.
template <std::size_t max_dim = 4UL>
inline Acts::BoundVariantMeasurement boundVariantMeasurement(
    const traccc::measurement& m, const Acts::SourceLink sl) {
  if constexpr (max_dim == 0UL) {
    std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                           std::to_string(m.meas_dim);
    throw std::runtime_error(errorMsg.c_str());
  } else {
    if (m.meas_dim == max_dim) {
      return measurement<max_dim>(m, sl);
    }
    return boundVariantMeasurement<max_dim - 1>(m, sl);
  }
}

/// @brief Gets the the local position of the measurement.
/// @param measurement the Acts measurement.
/// @returns A two-dimensional vector containing the local position.
/// The first item in the vector is local position on axis 0 and
/// I.e., [local position (axis 0), local position (axis 1)].
template <std::size_t dim>
inline Acts::ActsVector<2> getLocal(
    const Acts::Measurement<Acts::BoundIndices, dim>& measurement) {
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
inline Acts::ActsVector<2> getLocal(
    const Acts::BoundVariantMeasurement& measurement) {
  return std::visit([](auto& m) { return getLocal(m); }, measurement);
}

/// @brief Get the the variance of the measurement.
/// @param measurement the Acts measurement.
/// @return A two-dimensional vector containing the variance.
/// I.e., [variance (axis 0), variance (axis 1)].
template <std::size_t dim>
inline Acts::ActsVector<2> getVariance(
    const Acts::Measurement<Acts::BoundIndices, dim>& measurement) {
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
inline Acts::ActsVector<2> getVariance(
    const Acts::BoundVariantMeasurement& measurement) {
  return std::visit([](auto& m) { return getVariance(m); }, measurement);
}

}  // namespace Acts::TracccPlugin
