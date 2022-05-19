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
template <size_t M, bool ReadOnly>
inline TrackStateProxy<M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory, ReadOnly>& trajectory, size_t istate)
    : m_traj(&trajectory), m_istate(istate) {}

template <size_t M, bool ReadOnly>
TrackStatePropMask TrackStateProxy<M, ReadOnly>::getMask() const {
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

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::parameters() const -> Parameters {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
  } else if (hasFiltered()) {
    idx = data().ifiltered;
  } else {
    idx = data().ipredicted;
  }

  return m_traj->parameters(idx);
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::covariance() const -> Covariance {
  IndexData::IndexType idx;
  if (hasSmoothed()) {
    idx = data().ismoothed;
  } else if (hasFiltered()) {
    idx = data().ifiltered;
  } else {
    idx = data().ipredicted;
  }
  return m_traj->covariance(idx);
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::predicted() const -> Parameters {
  assert(data().ipredicted != IndexData::kInvalid);
  return m_traj->parameters(m_traj->template component<IndexData::IndexType>(
      hashString("predicted"), m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::predictedCovariance() const
    -> Covariance {
  assert(data().ipredicted != IndexData::kInvalid);
  // return Covariance(m_traj->m_cov[data().ipredicted].data());
  return m_traj->covariance(m_traj->template component<IndexData::IndexType>(
      hashString("predicted"), m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::filtered() const -> Parameters {
  assert(data().ifiltered != IndexData::kInvalid);
  return m_traj->parameters(m_traj->template component<IndexData::IndexType>(
      hashString("filtered"), m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::filteredCovariance() const
    -> Covariance {
  assert(data().ifiltered != IndexData::kInvalid);
  return m_traj->covariance(m_traj->template component<IndexData::IndexType>(
      hashString("filtered"), m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::smoothed() const -> Parameters {
  assert(data().ismoothed != IndexData::kInvalid);
  return m_traj->parameters(
      m_traj->template component<IndexData::IndexType>("smoothed", m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::smoothedCovariance() const
    -> Covariance {
  assert(data().ismoothed != IndexData::kInvalid);
  return m_traj->covariance(
      m_traj->template component<IndexData::IndexType>("smoothed", m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::jacobian() const -> Covariance {
  assert(data().ijacobian != IndexData::kInvalid);
  // return Covariance(m_traj->m_jac[data().ijacobian].data());
  return m_traj->jacobian(
      m_traj->template component<IndexData::IndexType>("jacobian", m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::projector() const -> Projector {
  assert(data().iprojector != IndexData::kInvalid);
  // return bitsetToMatrix<Projector>(m_traj->m_projectors[data().iprojector]);
  return bitsetToMatrix<Projector>(
      m_traj->template component<ProjectorBitset>("projector", m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::uncalibrated() const
    -> const SourceLink& {
  assert(data().iuncalibrated != IndexData::kInvalid);
  // assert(m_traj->m_sourceLinks[data().iuncalibrated] != nullptr);
  // return *m_traj->m_sourceLinks[data().iuncalibrated];
  using T = const SourceLink*;
  const T& sl =
      m_traj->template component<const SourceLink*>("sourceLink", m_istate);
  assert(sl != nullptr);
  return *sl;
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibrated() const -> Measurement {
  assert(data().icalibrated != IndexData::kInvalid);
  // return Measurement(m_traj->m_meas[data().icalibrated].data());
  // return m_traj->template component<Measurement>(hashString("measurement"),
  // m_istate);
  return m_traj->measurement(m_traj->template component<IndexData::IndexType>(
      "measurement", m_istate));
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibratedSourceLink() const
    -> const SourceLink& {
  assert(data().icalibratedsourcelink != IndexData::kInvalid);
  // assert(m_traj->m_sourceLinks[data().icalibratedsourcelink] != nullptr);
  // return *m_traj->m_sourceLinks[data().icalibratedsourcelink];
  const SourceLink*& sl = m_traj->template component<const SourceLink*>(
      "calibratedSourceLink", m_istate);
  assert(sl != nullptr);
  return *sl;
}

template <size_t M, bool ReadOnly>
inline auto TrackStateProxy<M, ReadOnly>::calibratedCovariance() const
    -> MeasurementCovariance {
  assert(data().icalibrated != IndexData::kInvalid);
  // return MeasurementCovariance(m_traj->m_measCov[data().icalibrated].data());
  // return m_traj->template component<MeasurementCovariance>(
  // hashString("measurementCovariance"), m_istate);
  return m_traj->measurementCovariance(
      m_traj->template component<IndexData::IndexType>("measurement",
                                                       m_istate));
}

}  // namespace detail_lt

template <typename F>
void MultiTrajectory::visitBackwards(size_t iendpoint, F&& callable) const {
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, ConstTrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (previous(iendpoint) == detail_lt::IndexData::kInvalid || !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (previous(iendpoint) == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = previous(iendpoint);
  }
}

template <typename F>
void MultiTrajectory::applyBackwards(size_t iendpoint, F&& callable) {
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    if constexpr (std::is_same_v<std::invoke_result_t<F, TrackStateProxy>,
                                 bool>) {
      bool proceed = callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory, or a break was
      // requested
      if (previous(iendpoint) == detail_lt::IndexData::kInvalid || !proceed) {
        break;
      }
    } else {
      callable(getTrackState(iendpoint));
      // this point has no parent and ends the trajectory
      if (previous(iendpoint) == detail_lt::IndexData::kInvalid) {
        break;
      }
    }
    iendpoint = previous(iendpoint);
  }
}
}  // namespace Acts
