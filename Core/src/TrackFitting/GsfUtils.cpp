// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "GsfUtils.hpp"

namespace Acts {
namespace detail {


auto extractMultiComponentState(const MultiTrajectory &traj,
                                const std::vector<size_t> &tips,
                                const std::map<size_t, ActsScalar> &weights,
                                StatesType type) -> MultiComponentBoundTrackParameters<SinglyCharged> {
  throw_assert(
      !tips.empty(),
      "need at least one component to extract trajectory of type " << type);

  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> cmps;
  std::shared_ptr<const Surface> surface;

  for (auto &tip : tips) {
    const auto proxy = traj.getTrackState(tip);

    throw_assert(weights.find(tip) != weights.end(),
                 "Could not find weight for idx " << tip);

    switch (type) {
      case StatesType::ePredicted:
        cmps.push_back(
            {weights.at(tip), proxy.predicted(), proxy.predictedCovariance()});
        break;
      case StatesType::eFiltered:
        cmps.push_back(
            {weights.at(tip), proxy.filtered(), proxy.filteredCovariance()});
        break;
      case StatesType::eSmoothed:
        cmps.push_back(
            {weights.at(tip), proxy.smoothed(), proxy.smoothedCovariance()});
    }

    if (!surface) {
      surface = proxy.referenceSurface().getSharedPtr();
    } else {
      throw_assert(
          surface->geometryId() == proxy.referenceSurface().geometryId(),
          "surface mismatch");
    }
  }

  return MultiComponentBoundTrackParameters<SinglyCharged>(surface, cmps);
}


}
}
