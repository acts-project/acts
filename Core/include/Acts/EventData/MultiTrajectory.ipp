// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
namespace detail_lt {
template <typename SL, size_t N, size_t M, bool ReadOnly>
inline TrackStateProxy<SL, N, M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory<SL>, ReadOnly>& trajectory, size_t istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <typename SL, size_t N, size_t M, bool ReadOnly>
TrackStatePropMask TrackStateProxy<SL, N, M, ReadOnly>::getMask() const {
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
  if (hasUncalibrated()) {
    mask |= PM::Uncalibrated;
  }
  if (hasCalibrated()) {
    mask |= PM::Calibrated;
  }

  return mask;
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::parameters() const
    -> Parameters {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
  } else if (hasFiltered()) {
    idx = data().ifiltered;
  } else {
    idx = data().ipredicted;
  }

  return Parameters(m_traj->m_params.data.col(idx).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::covariance() const
    -> Covariance {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
  } else if (hasFiltered()) {
    idx = data().ifiltered;
  } else {
    idx = data().ipredicted;
  }
  return Covariance(m_traj->m_cov.data.col(idx).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::predicted() const
    -> Parameters {
  assert(data().ipredicted != IndexData::kInvalid);
  return Parameters(m_traj->m_params.col(data().ipredicted).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::predictedCovariance() const
    -> Covariance {
  assert(data().ipredicted != IndexData::kInvalid);
  return Covariance(m_traj->m_cov.col(data().ipredicted).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::filtered() const
    -> Parameters {
  assert(data().ifiltered != IndexData::kInvalid);
  return Parameters(m_traj->m_params.col(data().ifiltered).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::filteredCovariance() const
    -> Covariance {
  assert(data().ifiltered != IndexData::kInvalid);
  return Covariance(m_traj->m_cov.col(data().ifiltered).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::smoothed() const
    -> Parameters {
  assert(data().ismoothed != IndexData::kInvalid);
  return Parameters(m_traj->m_params.col(data().ismoothed).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::smoothedCovariance() const
    -> Covariance {
  assert(data().ismoothed != IndexData::kInvalid);
  return Covariance(m_traj->m_cov.col(data().ismoothed).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::jacobian() const
    -> Covariance {
  assert(data().ijacobian != IndexData::kInvalid);
  return Covariance(m_traj->m_jac.col(data().ijacobian).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::projector() const
    -> Projector {
  assert(data().iprojector != IndexData::kInvalid);
  return bitsetToMatrix<Projector>(m_traj->m_projectors[data().iprojector]);
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::uncalibrated() const
    -> const SourceLink& {
  assert(data().iuncalibrated != IndexData::kInvalid);
  return m_traj->m_sourceLinks[data().iuncalibrated];
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::calibrated() const
    -> Measurement {
  assert(data().icalibrated != IndexData::kInvalid);
  return Measurement(m_traj->m_meas.col(data().icalibrated).data());
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::calibratedSourceLink() const
    -> const SourceLink& {
  assert(data().icalibratedsourcelink != IndexData::kInvalid);
  return m_traj->m_sourceLinks[data().icalibratedsourcelink];
}

template <typename SL, size_t N, size_t M, bool ReadOnly>
inline auto TrackStateProxy<SL, N, M, ReadOnly>::calibratedCovariance() const
    -> MeasurementCovariance {
  assert(data().icalibrated != IndexData::kInvalid);
  return MeasurementCovariance(
      m_traj->m_measCov.col(data().icalibrated).data());
}

}  // namespace detail_lt

template <typename SL>
inline size_t MultiTrajectory<SL>::addTrackState(TrackStatePropMask mask,
                                                 size_t iprevious) {
  using PropMask = TrackStatePropMask;

  m_index.emplace_back();
  detail_lt::IndexData& p = m_index.back();
  size_t index = m_index.size() - 1;

  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
  }

  // always set, but can be null
  m_referenceSurfaces.emplace_back(nullptr);
  p.irefsurface = m_referenceSurfaces.size() - 1;

  if (ACTS_CHECK_BIT(mask, PropMask::Predicted)) {
    m_params.addCol();
    m_cov.addCol();
    p.ipredicted = m_params.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Filtered)) {
    m_params.addCol();
    m_cov.addCol();
    p.ifiltered = m_params.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Smoothed)) {
    m_params.addCol();
    m_cov.addCol();
    p.ismoothed = m_params.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Jacobian)) {
    m_jac.addCol();
    p.ijacobian = m_jac.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Uncalibrated)) {
    m_sourceLinks.emplace_back();
    p.iuncalibrated = m_sourceLinks.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Calibrated)) {
    m_meas.addCol();
    m_measCov.addCol();
    p.icalibrated = m_meas.size() - 1;

    m_sourceLinks.emplace_back();
    p.icalibratedsourcelink = m_sourceLinks.size() - 1;

    m_projectors.emplace_back();
    p.iprojector = m_projectors.size() - 1;
  }

  return index;
}

template <typename SL>
template <typename F>
void MultiTrajectory<SL>::visitBackwards(size_t iendpoint, F&& callable) const {
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, ConstTrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid ||
          !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

template <typename SL>
template <typename F>
void MultiTrajectory<SL>::applyBackwards(size_t iendpoint, F&& callable) {
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid ||
          !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}
}  // namespace Acts
