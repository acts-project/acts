// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/detail/WeightedChannelCombiner.hpp"
#include <Acts/EventData/ParameterSet.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/SurfaceError.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

namespace ActsFatras {

/// Channel combination to a resulting parameter set.
///
struct ChannelCombiner {
  /// Generic implementation of a channel combiner, the channels are
  /// expected to be clustered and pre-processed, only a simple channel
  /// combination to a measurement is done.
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @param channels The channels from one cluster
  ///
  /// @return A cluster containing the parameter set and cluster size
  template <Acts::BoundIndices... kParameters>
  Cluster<kParameters...> combine(
      const std::vector<const Channel<kParameters...>>& channels) const {
    std::array<unsigned int, sizeof...(kParameters)> cSize;
    using ParSet = Acts::ParameterSet<Acts::BoundIndices, kParameters...>;

    typename ParSet::ParametersVector cVec;
    cVec.setZero();
    typename ParSet::CovarianceMatrix cCov;

    detail::WeightedChannelCombiner<kParameters...> combiner;
    combiner.run(cVec, cSize, channels);

    return Cluster<kParameters...>(ParSet(std::move(cCov), std::move(cVec)),
                                   std::move(cSize), channels);
  }
};

}  // namespace ActsFatras
