// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

namespace Acts {
namespace detail_lt {
template <typename D, size_t M, bool ReadOnly>
inline TrackStateProxy<D, M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory<D>, ReadOnly>& trajectory, IndexType istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <typename D, size_t M, bool ReadOnly>
TrackStatePropMask TrackStateProxy<D, M, ReadOnly>::getMask() const {
  using PM = TrackStatePropMask;

  PM mask = PM::None;
  if (hasPredicted()) {
    mask |= PM::Predicted;
  }
  if (hasFiltered()) {
    mask |= PM::Filtered;
  }
  if (hasSmoothed()) {
    mask |= PM::Smoothed;
  }
  if (hasJacobian()) {
    mask |= PM::Jacobian;
  }
  if (hasCalibrated()) {
    mask |= PM::Calibrated;
  }
  return mask;
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::parameters() const -> Parameters {
  if (hasSmoothed()) {
    return smoothed();
  } else if (hasFiltered()) {
    return filtered();
  } else {
    return predicted();
  }
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::covariance() const -> Covariance {
  if (hasSmoothed()) {
    return smoothedCovariance();
  } else if (hasFiltered()) {
    return filteredCovariance();
  } else {
    return predictedCovariance();
  }
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::projector() const -> Projector {
  assert(has<hashString("projector")>());
  return bitsetToMatrix<Projector>(
      component<ProjectorBitset, hashString("projector")>());
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::uncalibratedSourceLink() const
    -> const SourceLink& {
  assert(has<hashString("uncalibratedSourceLink")>());
  return component<std::optional<SourceLink>,
                   hashString("uncalibratedSourceLink")>()
      .value();
}

template <typename D, size_t M, bool ReadOnly>
inline auto TrackStateProxy<D, M, ReadOnly>::calibratedSourceLink() const
    -> const SourceLink& {
  assert(has<hashString("calibratedSourceLink")>());
  return component<std::optional<SourceLink>,
                   hashString("calibratedSourceLink")>()
      .value();
}

}  // namespace detail_lt

template <typename D>
template <typename F>
void MultiTrajectory<D>::visitBackwards(IndexType iendpoint,
                                        F&& callable) const {
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  if (iendpoint == MultiTrajectoryTraits::kInvalid) {
    throw std::runtime_error(
        "Cannot visit backwards with kInvalid as endpoint");
  }

  while (true) {
    auto ts = getTrackState(iendpoint);
    if constexpr (std::is_same_v<std::invoke_result_t<F, ConstTrackStateProxy>,
                                 bool>) {
      bool proceed = callable(ts);
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (!proceed || !ts.hasPrevious()) {
        break;
      }
    } else {
      callable(ts);
      // this point has no parent and ends the trajectory
      if (!ts.hasPrevious()) {
        break;
      }
    }
    iendpoint = ts.previous();
  }
}

}  // namespace Acts
