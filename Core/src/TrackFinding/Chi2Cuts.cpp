// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/Chi2Cuts.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

#include <boost/math/distributions/chi_squared.hpp>

namespace Acts {

Chi2CutPerDimension chi2CutsWithSameProbability(double cut1D) {
  Chi2CutPerDimension result{};
  result.fill(cut1D);
  if (!std::isfinite(cut1D) || cut1D == std::numeric_limits<double>::max() ||
      cut1D <= 0) {
    return result;
  }

  const double pValue = boost::math::cdf(boost::math::complement(
      boost::math::chi_squared_distribution<>(1), cut1D));
  if (!(pValue > std::numeric_limits<double>::min())) {
    return result;
  }
  for (std::size_t dim = 2; dim < result.size(); ++dim) {
    result[dim] = boost::math::quantile(boost::math::complement(
        boost::math::chi_squared_distribution<>(static_cast<double>(dim)),
        pValue));
  }
  return result;
}

}  // namespace Acts
