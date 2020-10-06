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
struct ChannelCreator {
  /// Smearing functions definition:
  ///
  /// @tparam generator_t The type of the random generator.
  ///
  /// - it takes the unsmeared parameter
  /// - it returns the smeared parameter and a covariance
  template <typename generator_t>
  using SmearFunction = std::function<Acts::Result<std::pair<double, double>>(
      double, generator_t&)>;

  /// Generic implementation of a channel creator.
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @return A vector of channels that are geometrically hit
  template <typename signal_t, typename generator_t,
            Acts::BoundIndices... kParameters>
  const std::vector<Channel<signal_t, kParameters...>> channels(
      const DigitizationInput& sInput, generator_t& sRandom,
      SmearFunction<generator_t>& sFunction) const {
    std::vector<Channel<signal_t, kParameters...>> channels;

    return channels;
  }
};

}  // namespace ActsFatras
