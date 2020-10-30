// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include <Acts/EventData/ParameterSet.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/SurfaceError.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

#include <array>
#include <utility>

namespace ActsFatras {

/// Channel combination to a resulting parameter set.
///
struct ChannelMerger {
  /// Generic implementation of a channel merger, currently only additive
  /// channel merging
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @param channels The channels from one cluster
  ///
  /// @return A cluster containing the parameter set and cluster size
  template <typename signal_t, Acts::BoundIndices... kParameters>
  const std::vector<Channel<signal_t, kParameters...>> merge(
      const std::vector<Channel<signal_t, kParameters...>>& channels) const {
    using StorageSequence = std::make_index_sequence<sizeof...(kParameters)>;
    using Channel = Channel<signal_t, kParameters...>;
    using ChannelKey = std::array<unsigned int, sizeof...(kParameters)>;

    // Fill a channel map
    std::map<ChannelKey, Channel> channelMap;
    for (const auto& ch : channels) {
      ChannelKey key;
      extractKeys(ch, key, StorageSequence{});
      auto chItr = channelMap.find(key);
      if (chItr != channelMap.end()) {
        chItr->second.value += ch.value;
        chItr->second.links.insert(ch.links.begin(), ch.links.end());
      } else {
        channelMap.insert(std::pair<ChannelKey, Channel>(key, ch));
      }
    }
    // Unroll the channels after merging
    std::vector<Channel> mergedChannels;
    mergedChannels.reserve(channelMap.size());
    for (auto& [key, value] : channelMap) {
      mergedChannels.push_back(value);
    }
    return mergedChannels;
  }

  /// Helper method to extract the map keys from the channels
  ///
  /// @tparam kStorage is the access sequence to the channel
  template <typename channel_t, typename key_t, std::size_t... kStorage>
  static void extractKeys(const channel_t& channel, key_t& key,
                          std::index_sequence<kStorage...>) {
    ((key[kStorage] = channel.cellId[kStorage].first), ...);
  }
};

}  // namespace ActsFatras
