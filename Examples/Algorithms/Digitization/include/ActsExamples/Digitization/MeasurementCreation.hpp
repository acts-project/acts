// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

namespace ActsExamples {

/// Struct to identify digitized parameters
///
/// Is public as it is also used by the I/O system
struct DigitizedParameters {
  std::vector<Acts::BoundIndices> indices = {};
  std::vector<Acts::ActsScalar> values = {};
  std::vector<Acts::ActsScalar> variances = {};

  Cluster cluster;
};

/// Helper method for created a measurement from digitized parameters
///
/// @param dParams The digitized parameters of variable size
/// @param isl The indexed source link for the measurement
///
/// To be used also by the e I/O system
///
/// @return a variant measurement
Measurement createMeasurement(const DigitizedParameters& dParams,
                              const IndexSourceLink& isl) noexcept(false);

/// Contruct the constituents of a measurement.
///
/// @tparam kMeasDIM the full dimension of the measurement
///
/// @param dParams the struct of arrays of parameters to be created
///
/// @return a tuple of constituents for a measurement
template <size_t kMeasDIM>
std::tuple<std::array<Acts::BoundIndices, kMeasDIM>, Acts::ActsVector<kMeasDIM>,
           Acts::ActsSymMatrix<kMeasDIM>>
measurementConstituents(const DigitizedParameters& dParams) {
  std::array<Acts::BoundIndices, kMeasDIM> indices{};
  Acts::ActsVector<kMeasDIM> par;
  Acts::ActsSymMatrix<kMeasDIM> cov = Acts::ActsSymMatrix<kMeasDIM>::Identity();
  for (Eigen::Index ei = 0; ei < static_cast<Eigen::Index>(kMeasDIM); ++ei) {
    indices[ei] = dParams.indices[ei];
    par[ei] = dParams.values[ei];
    cov(ei, ei) = dParams.variances[ei];
  }
  return {indices, par, cov};
}

}  // namespace ActsExamples
