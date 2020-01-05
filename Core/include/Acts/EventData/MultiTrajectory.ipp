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

#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
namespace detail_lt {
template <typename SL, size_t N, size_t M, bool ReadOnly>
inline TrackStateProxy<SL, N, M, ReadOnly>::TrackStateProxy(
    ConstIf<MultiTrajectory<SL>, ReadOnly>& trajectory, size_t istate)
    : m_traj(&trajectory), m_istate(istate) {}

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
inline BoundParameters TrackStateProxy<SL, N, M, ReadOnly>::predictedParameters(
    const Acts::GeometryContext& gctx) const {
  return {gctx, predictedCovariance(), predicted(),
          m_traj->m_referenceSurfaces[data().irefsurface]};
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
inline BoundParameters TrackStateProxy<SL, N, M, ReadOnly>::filteredParameters(
    const Acts::GeometryContext& gctx) const {
  return {gctx, filteredCovariance(), filtered(),
          m_traj->m_referenceSurfaces[data().irefsurface]};
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
inline BoundParameters TrackStateProxy<SL, N, M, ReadOnly>::smoothedParameters(
    const Acts::GeometryContext& gctx) const {
  return {gctx, smoothedCovariance(), smoothed(),
          m_traj->m_referenceSurfaces[data().irefsurface]};
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
template <typename parameters_t>
inline size_t MultiTrajectory<SL>::addTrackState(
    const TrackState<SL, parameters_t>& ts, size_t iprevious) {
  using CovMap =
      typename detail_lt::Types<ParametersSize, false>::CovarianceMap;

  // use a TrackStateProxy to do the assignments
  m_index.emplace_back();
  detail_lt::IndexData& p = m_index.back();
  size_t index = m_index.size() - 1;

  TrackStateProxy nts = getTrackState(index);

  // make shared ownership held by this multi trajectory
  m_referenceSurfaces.push_back(ts.referenceSurface().getSharedPtr());
  p.irefsurface = m_referenceSurfaces.size() - 1;

  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
    TrackStateProxy previous = getTrackState(iprevious);
    previous.data().nchildren++;
  }

  if (ts.parameter.predicted) {
    const auto& predicted = *ts.parameter.predicted;
    m_params.addCol() = predicted.parameters();
    CovMap(m_cov.addCol().data()) = *predicted.covariance();
    p.ipredicted = m_params.size() - 1;
  }

  if (ts.parameter.filtered) {
    const auto& filtered = *ts.parameter.filtered;
    m_params.addCol() = filtered.parameters();
    CovMap(m_cov.addCol().data()) = *filtered.covariance();
    p.ifiltered = m_params.size() - 1;
  }

  if (ts.parameter.smoothed) {
    const auto& smoothed = *ts.parameter.smoothed;
    m_params.addCol() = smoothed.parameters();
    CovMap(m_cov.addCol().data()) = *smoothed.covariance();
    p.ismoothed = m_params.size() - 1;
  }

  // store jacobian
  if (ts.parameter.jacobian) {
    CovMap(m_jac.addCol().data()) = *ts.parameter.jacobian;
    p.ijacobian = m_jac.size() - 1;
  }

  // handle measurements
  if (ts.measurement.uncalibrated) {
    m_sourceLinks.push_back(*ts.measurement.uncalibrated);
    p.iuncalibrated = m_sourceLinks.size() - 1;
  }

  if (ts.measurement.calibrated) {
    std::visit([&](const auto& m) { nts.resetCalibrated(m); },
               *ts.measurement.calibrated);
  }

  nts.chi2() = ts.parameter.chi2;
  nts.pathLength() = ts.parameter.pathLength;
  nts.typeFlags() = ts.type();

  return index;
}

template <typename SL>
inline size_t MultiTrajectory<SL>::addTrackState(
    const TrackStatePropMask::Type& mask, size_t iprevious) {
  namespace PropMask = TrackStatePropMask;

  m_index.emplace_back();
  detail_lt::IndexData& p = m_index.back();
  size_t index = m_index.size() - 1;

  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
    TrackStateProxy previous = getTrackState(iprevious);
    previous.data().nchildren++;
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
