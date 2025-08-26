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
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <array>
#include <cstddef>
#include <tuple>
#include <vector>

namespace ActsExamples {

/// Struct to identify digitized parameters
///
/// Is public as it is also used by the I/O system
struct DigitizedParameters {
  std::vector<Acts::BoundIndices> indices = {};
  std::vector<double> values = {};
  std::vector<double> variances = {};

  Cluster cluster;
};

/// Helper method for created a measurement from digitized parameters
///
/// @param container The measurement container to insert into
/// @param geometryId The geometry ID of the measurement surface
/// @param dParams The digitized parameters of variable size
///
/// To be used also by the e I/O system
///
/// @return the measurement proxy
ActsExamples::VariableBoundMeasurementProxy createMeasurement(
    MeasurementContainer& container, Acts::GeometryIdentifier geometryId,
    const DigitizedParameters& dParams) noexcept(false);

/// Construct the constituents of a measurement.
///
/// @tparam kMeasDIM the full dimension of the measurement
///
/// @param dParams the struct of arrays of parameters to be created
///
/// @return a tuple of constituents for a measurement
template <std::size_t kMeasDIM>
std::tuple<std::array<Acts::BoundIndices, kMeasDIM>, Acts::ActsVector<kMeasDIM>,
           Acts::ActsSquareMatrix<kMeasDIM>>
measurementConstituents(const DigitizedParameters& dParams) {
  std::array<Acts::BoundIndices, kMeasDIM> indices{};
  Acts::ActsVector<kMeasDIM> par;
  Acts::ActsSquareMatrix<kMeasDIM> cov =
      Acts::ActsSquareMatrix<kMeasDIM>::Identity();
  for (Eigen::Index ei = 0; ei < static_cast<Eigen::Index>(kMeasDIM); ++ei) {
    indices[ei] = dParams.indices[ei];
    par[ei] = dParams.values[ei];
    cov(ei, ei) = dParams.variances[ei];
  }
  return {indices, par, cov};
}

/// Function that computes a global position for a measurement.
/// For 1D measurements, the center of the module is used for the missing
/// dimension.
inline Acts::Vector3 measurementGlobalPosition(
    const DigitizedParameters& digitizedParameters,
    const Acts::Surface& surface, const Acts::GeometryContext& gctx) {
  // we need a regular surface to perform the local to global transformation
  // without direction input. the direction could be obtained from truth
  // information but we can leave it out here.
  const Acts::RegularSurface* regularSurface =
      dynamic_cast<const Acts::RegularSurface*>(&surface);
  if (regularSurface == nullptr) {
    throw std::invalid_argument("Expected a regular surface");
  }

  Acts::Vector2 locPos = Acts::Vector2::Zero();
  for (auto i = 0ul; i < digitizedParameters.indices.size(); ++i) {
    auto idx = digitizedParameters.indices.at(i);
    if (idx == Acts::eBoundLoc0 || idx == Acts::eBoundLoc1) {
      locPos[idx] = digitizedParameters.values.at(i);
    } else {
      locPos[idx] = regularSurface->bounds().center()[idx];
    }
  }

  return regularSurface->localToGlobal(gctx, locPos);
}

}  // namespace ActsExamples
