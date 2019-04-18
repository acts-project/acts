// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
namespace detail_lt {
  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline TrackStateProxy<SL, N, M, ReadOnly>::TrackStateProxy(
      ConstIf<MultiTrajectory<SL>, ReadOnly>& trajectory,
      size_t istate)
    : m_traj(trajectory), m_istate(istate), m_data(trajectory.m_index[istate])
  {
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::parameters() const -> Parameters
  {
    IndexData::IndexType idx;
    if (hasSmoothed()) {
      idx = m_data.ismoothed;
    } else if (hasFiltered()) {
      idx = m_data.ifiltered;
    } else {
      idx = m_data.ipredicted;
    }

    return Parameters(m_traj.m_params.data.col(idx).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::covariance() const -> Covariance
  {
    IndexData::IndexType idx;
    if (hasSmoothed()) {
      idx = m_data.ismoothed;
    } else if (hasFiltered()) {
      idx = m_data.ifiltered;
    } else {
      idx = m_data.ipredicted;
    }
    return Covariance(m_traj.m_cov.data.col(idx).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::predicted() const -> Parameters
  {
    assert(m_data.ipredicted != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::predictedCovariance() const -> Covariance
  {
    assert(m_data.ipredicted != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filtered() const -> Parameters
  {
    assert(m_data.ifiltered != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filteredCovariance() const -> Covariance
  {
    assert(m_data.ifiltered != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothed() const -> Parameters
  {
    assert(m_data.ismoothed != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothedCovariance() const -> Covariance
  {
    assert(m_data.ismoothed != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::jacobian() const -> Covariance
  {
    assert(m_data.ijacobian != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ijacobian).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::projector() const -> Projector
  {
    assert(m_data.iprojector != IndexData::kInvalid);
    return bitsetToMatrix<Projector>(m_traj.m_projectors[m_data.iprojector]);
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::uncalibrated() const -> const SourceLink&
  {
    assert(m_data.iuncalibrated != IndexData::kInvalid);
    return m_traj.m_sourceLinks[m_data.iuncalibrated];
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibrated() const -> Measurement
  {
    assert(m_data.icalibrated != IndexData::kInvalid);
    return Measurement(m_traj.m_meas.col(m_data.icalibrated).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibratedSourceLink() const
      -> const SourceLink&
  {
    assert(m_data.icalibratedsourcelink != IndexData::kInvalid);
    return m_traj.m_sourceLinks[m_data.icalibratedsourcelink];
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibratedCovariance() const
      -> MeasurementCovariance
  {
    assert(m_data.icalibrated != IndexData::kInvalid);
    return MeasurementCovariance(
        m_traj.m_measCov.col(m_data.icalibrated).data());
  }

}  // namespace detail_lt

template <typename SL>
template <typename parameters_t>
inline size_t
MultiTrajectory<SL>::addTrackState(const TrackState<SL, parameters_t>& ts,
                                   size_t iprevious)
{
  using CovMap =
      typename detail_lt::Types<ParametersSize, false>::CovarianceMap;

  detail_lt::IndexData p = {ts.referenceSurface()};
  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
  }

  if (ts.parameter.predicted) {
    const auto& predicted         = *ts.parameter.predicted;
    m_params.addCol()             = predicted.parameters();
    CovMap(m_cov.addCol().data()) = *predicted.covariance();
    p.ipredicted                  = m_params.size() - 1;
  }

  if (ts.parameter.filtered) {
    const auto& filtered          = *ts.parameter.filtered;
    m_params.addCol()             = filtered.parameters();
    CovMap(m_cov.addCol().data()) = *filtered.covariance();
    p.ifiltered                   = m_params.size() - 1;
  }

  if (ts.parameter.smoothed) {
    const auto& smoothed          = *ts.parameter.smoothed;
    m_params.addCol()             = smoothed.parameters();
    CovMap(m_cov.addCol().data()) = *smoothed.covariance();
    p.ismoothed                   = m_params.size() - 1;
  }

  // store jacobian
  if (ts.parameter.jacobian) {
    CovMap(m_cov.addCol().data()) = *ts.parameter.jacobian;
    p.ijacobian                   = m_cov.size() - 1;
  }

  // handle measurements
  if (ts.measurement.uncalibrated) {
    m_sourceLinks.push_back(*ts.measurement.uncalibrated);
    p.iuncalibrated = m_sourceLinks.size() - 1;
  }

  if (ts.measurement.calibrated) {
    auto meas    = m_meas.addCol();
    auto measCov = m_measCov.addCol();
    std::visit(
        [&meas, &measCov, &p, this](const auto& m) {
          using meas_t                         = std::decay_t<decltype(m)>;
          meas.template head<meas_t::size()>() = m.parameters();
          CovMap(measCov.data())
              .template topLeftCorner<meas_t::size(), meas_t::size()>()
              = m.covariance();

          // We can only store the projection if we have a calibrated
          // measurement. Place (possibly asymmetric) projector into
          // full size projector, padded with zeroes.
          // Convert to bitset before setting.
          typename TrackStateProxy::Projector fullProjector;
          fullProjector.setZero();
          fullProjector
              .template topLeftCorner<meas_t::size(), MeasurementSizeMax>()
              = m.projector();

          m_projectors.push_back(matrixToBitset(fullProjector));

          // we also need to set the measurement dimension
          p.measdim = meas_t::size();

          m_sourceLinks.push_back(m.sourceLink());
          p.icalibratedsourcelink = m_sourceLinks.size() - 1;
        },
        *ts.measurement.calibrated);
    p.icalibrated = m_meas.size() - 1;
    p.iprojector  = m_projectors.size() - 1;
  }

  m_index.push_back(std::move(p));
  return m_index.size() - 1;
}

template <typename SL>
template <typename F>
void
MultiTrajectory<SL>::visitBackwards(size_t iendpoint, F&& callable) const
{
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    callable(getTrackState(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

template <typename SL>
template <typename F>
void
MultiTrajectory<SL>::applyBackwards(size_t iendpoint, F&& callable)
{
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");
  while (true) {
    callable(getTrackState(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}
}
