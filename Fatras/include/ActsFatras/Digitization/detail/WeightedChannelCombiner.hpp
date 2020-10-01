// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/DigitizationData.hpp"

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

#include <climits>
#include <utility>

namespace ActsFatras {

namespace detail {

/// Combine the channels
///
/// @tparam kParameters Parameter indices pack that defines the used parameter
template <typename signal_t, Acts::BoundIndices... kParameters>
struct WeightedChannelCombiner {
  using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;

  /// Run the channel combination over the parameter pack
  ///
  /// @tparam values_t The type of the values vector
  ///
  /// @param values [in,out] Combined parameters, initialized input expected
  /// @param cSize [in,out] The cluster size estimation, done on the fly
  /// @param channels The channel list to combine
  ///
  /// This function runs over the channels and creates a weighted mean to
  /// estiamte the cluster position. Digital clustering can be achieved by
  /// providing a channel list with equal weight for each channel.
  ///
  /// At the same time the cluster size is estimated by finding the channel
  /// extremas min/max and fill the difference
  template <typename values_t>
  static void run(
      Eigen::MatrixBase<values_t>& values,
      std::array<unsigned int, sizeof...(kParameters)>& cSize,
      const std::vector<Channel<signal_t, kParameters...>>& channels) {
    runImpl(values, cSize, channels, StorageSequence{});
  }

  template <typename values_t, std::size_t... kStorage>
  static void runImpl(
      Eigen::MatrixBase<values_t>& values,
      std::array<unsigned int, sizeof...(kParameters)>& cSize,
      const std::vector<Channel<signal_t, kParameters...>>& channels,
      std::index_sequence<kStorage...>) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(values_t, sizeof...(kParameters));
    static_assert(sizeof...(kParameters) == sizeof...(kStorage),
                  "Parameters and storage index packs must have the same size");

    double fullWeight = 0.;
    std::array<std::pair<int, int>, sizeof...(kParameters)> cExtremas;
    ((cExtremas[kStorage] = {std::numeric_limits<int>::max(),
                             std::numeric_limits<int>::lowest()}),
     ...);

    for (auto& ch : channels) {
      double weight = ch.value;
      ((values[kStorage] += weight * ch.cellId[kStorage].second), ...);
      ((cExtremas[kStorage].first =
            std::min(cExtremas[kStorage].first,
                     static_cast<int>(ch.cellId[kStorage].first))),
       ...);
      ((cExtremas[kStorage].second =
            std::max(cExtremas[kStorage].second,
                     static_cast<int>(ch.cellId[kStorage].first))),
       ...);
      fullWeight += weight;
    }

    // Normalize and & estimate size
    ((values[kStorage] = values[kStorage] / fullWeight), ...);
    ((cSize[kStorage] = 1 + (unsigned int)(cExtremas[kStorage].second -
                                           cExtremas[kStorage].first)),
     ...);
  }
};

}  // namespace detail
}  // namespace ActsFatras
