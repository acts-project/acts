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
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <ranges>
#include <span>

namespace Acts {

/// A calibrated measurement in the bound parameter space, stored inline.
///
/// The measurement dimension is a runtime value in `[1, eBoundSize]`, but the
/// storage is allocated for the maximum size, so producing one never touches
/// the heap. This makes it suitable as the return value of a per-measurement
/// calibration which is called in the inner loop of track finding.
///
/// Only the leading `size()` entries of the parameter storage and the leading
/// `size() x size()` block of the (column major, stride `size()`) covariance
/// storage are meaningful.
class CalibratedBoundMeasurement {
 public:
  /// Construct an invalid, zero dimensional measurement.
  CalibratedBoundMeasurement() = default;

  /// Construct from fixed size parameters and covariance.
  ///
  /// The measurement dimension is deduced from @p parameters, so it has to be
  /// known at compile time. Use @ref visit_measurement to get there from a
  /// runtime dimension.
  ///
  /// @tparam parameters_t the parameter vector type
  /// @tparam covariance_t the covariance matrix type
  /// @tparam index_range_t the subspace index range type
  ///
  /// @param parameters the measured parameters
  /// @param covariance the covariance of the measured parameters
  /// @param subspaceIndices the bound indices the measurement constrains
  template <typename parameters_t, typename covariance_t,
            std::ranges::sized_range index_range_t>
  CalibratedBoundMeasurement(const Eigen::DenseBase<parameters_t>& parameters,
                             const Eigen::DenseBase<covariance_t>& covariance,
                             const index_range_t& subspaceIndices) {
    constexpr int kSize = parameters_t::RowsAtCompileTime;
    static_assert(kSize != Eigen::Dynamic,
                  "Measurement dimension must be known at compile time");
    static_assert(kSize >= 1 && kSize <= static_cast<int>(eBoundSize),
                  "Invalid measurement dimension");
    assert(static_cast<int>(subspaceIndices.size()) == kSize &&
           "Subspace size mismatch");

    m_size = kSize;
    Eigen::Map<Vector<kSize>>{m_parameters.data()} = parameters;
    Eigen::Map<SquareMatrix<kSize>>{m_covariance.data()} = covariance;
    std::ranges::transform(
        subspaceIndices, m_subspaceIndices.begin(),
        [](auto index) { return static_cast<SubspaceIndex>(index); });
  }

  /// The measurement dimension.
  /// @return the number of constrained bound parameters
  std::uint8_t size() const { return m_size; }

  /// Pointer to the contiguous parameter storage.
  /// @return pointer to `size()` contiguous values
  const double* parametersData() const { return m_parameters.data(); }

  /// Pointer to the contiguous covariance storage.
  /// @return pointer to the column major `size() x size()` block
  const double* covarianceData() const { return m_covariance.data(); }

  /// Fixed size view of the measured parameters.
  /// @tparam kSize the measurement dimension, must equal @ref size
  /// @return a map onto the parameter storage
  template <std::size_t kSize>
  Eigen::Map<const Vector<kSize>> parameters() const {
    assert(kSize == m_size && "Size mismatch");
    return Eigen::Map<const Vector<kSize>>{m_parameters.data()};
  }

  /// Fixed size view of the measurement covariance.
  /// @tparam kSize the measurement dimension, must equal @ref size
  /// @return a map onto the covariance storage
  template <std::size_t kSize>
  Eigen::Map<const SquareMatrix<kSize>> covariance() const {
    assert(kSize == m_size && "Size mismatch");
    return Eigen::Map<const SquareMatrix<kSize>>{m_covariance.data()};
  }

  /// The bound indices the measurement constrains, exactly @ref size of them.
  /// @return a view of the subspace indices
  std::span<const SubspaceIndex> subspaceIndices() const {
    return {m_subspaceIndices.data(), m_size};
  }

  /// The subspace indices zero padded to `eBoundSize`.
  ///
  /// This is the form expected by @ref calculatePredictedChi2, which only
  /// looks at the leading @ref size entries.
  ///
  /// @return the padded subspace indices
  const BoundSubspaceIndices& paddedSubspaceIndices() const {
    return m_subspaceIndices;
  }

  /// Write the calibrated data into a track state.
  ///
  /// @note This passes exactly @ref size subspace indices, matching what the
  ///       calibrators do. The projector is serialized including its unused
  ///       slots, so passing a padded range here would not be equivalent.
  ///
  /// @tparam track_state_proxy_t the track state proxy type
  /// @param trackState the track state to fill
  template <typename track_state_proxy_t>
  void fill(track_state_proxy_t trackState) const {
    visit_measurement(m_size, [&]<std::size_t kSize>(
                                  std::integral_constant<std::size_t, kSize>) {
      trackState.allocateCalibrated(parameters<kSize>().eval(),
                                    covariance<kSize>().eval());
    });
    trackState.setProjectorSubspaceIndices(subspaceIndices());
  }

 private:
  std::uint8_t m_size{0};
  /// Zero padded beyond `m_size`, matching what the calibrators store.
  BoundSubspaceIndices m_subspaceIndices{};
  std::array<double, eBoundSize> m_parameters{};
  /// Column major with stride `m_size`.
  std::array<double, eBoundSize * eBoundSize> m_covariance{};
};

}  // namespace Acts
