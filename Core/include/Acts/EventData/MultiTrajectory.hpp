// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

// forward declarations
template <typename source_link_t>
class MultiTrajectory;
class Surface;

namespace detail_lt {
  /// Either type T or const T depending on the boolean.
  template <typename T, bool select>
  using ConstIf = std::conditional_t<select, const T, T>;
  /// wrapper for a dynamic Eigen type that adds support for automatic growth
  ///
  /// \warning Assumes the underlying storage has a fixed number of rows
  template <typename Storage, size_t kSizeIncrement>
  struct GrowableColumns
  {

    /// Make sure storage for @p n additional columns is allocated. Will update
    /// the size of the container accordingly. The indices added by this call
    /// can safely be written to.
    /// @param n Number of columns to add, defaults to 1.
    /// @return View into the last allocated column
    auto
    addCol(size_t n = 1)
    {
      size_t index = m_size + (n - 1);
      while (capacity() <= index) {
        data.conservativeResize(Eigen::NoChange, data.cols() + kSizeIncrement);
      }
      m_size = index + 1;
      return data.col(index);
    }

    /// Writable access to a column w/o checking its existence first.
    auto
    col(size_t index)
    {
      return data.col(index);
    }

    /// Read-only access to a column w/o checking its existence first.
    auto
    col(size_t index) const
    {
      return data.col(index);
    }

    /// Return the current allocated storage capacity
    size_t
    capacity() const
    {
      return static_cast<size_t>(data.cols());
    }

    size_t
    size() const
    {
      return m_size;
    }

  private:
    Storage data;
    size_t  m_size{0};
  };

  /// Type construction helper for coefficients and associated covariances.
  template <size_t Size, bool ReadOnlyMaps = true>
  struct Types
  {
    enum {
      Flags         = Eigen::ColMajor | Eigen::AutoAlign,
      SizeIncrement = 8,
    };
    using Scalar = double;
    // single items
    using Coefficients    = Eigen::Matrix<Scalar, Size, 1, Flags>;
    using Covariance      = Eigen::Matrix<Scalar, Size, Size, Flags>;
    using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
    using CovarianceMap   = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;
    // storage of multiple items in flat arrays
    using StorageCoefficients
        = GrowableColumns<Eigen::Array<Scalar, Size, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
    using StorageCovariance
        = GrowableColumns<Eigen::
                              Array<Scalar, Size * Size, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
  };
  struct IndexData
  {
    using IndexType = uint16_t;

    static constexpr IndexType kInvalid = UINT16_MAX;

    const Surface& surface;
    IndexType      iprevious     = kInvalid;
    IndexType      ipredicted    = kInvalid;
    IndexType      ifiltered     = kInvalid;
    IndexType      ismoothed     = kInvalid;
    IndexType      iuncalibrated = kInvalid;
    IndexType      icalibrated   = kInvalid;
    IndexType      measdim       = 0;
  };
  /// Proxy object to access a single point on the trajectory.
  ///
  /// @tparam source_link_t Type to link back to an original measurement
  /// @tparam N         Number of track parameters
  /// @tparam M         Maximum number of measurements
  /// @tparam ReadOnly  true for read-only access to underlying storage
  template <typename source_link_t, size_t N, size_t M, bool ReadOnly = true>
  class TrackStateProxy
  {
  public:
    using SourceLink            = source_link_t;
    using Parameters            = typename Types<N, ReadOnly>::CoefficientsMap;
    using Covariance            = typename Types<N, ReadOnly>::CovarianceMap;
    using Measurement           = typename Types<M, ReadOnly>::CoefficientsMap;
    using MeasurementCovariance = typename Types<M, ReadOnly>::CovarianceMap;

    /// Index within the trajectory.
    size_t
    index() const
    {
      return m_istate;
    }

    /// Reference surface.
    const Surface&
    referenceSurface() const
    {
      return m_data.surface;
    }

    /// Track parameters vector.
    Parameters
    parameters() const;

    /// Track parameters covariance matrix.
    Covariance
    covariance() const;

    /// Track parameters vector.
    Parameters
    predicted() const;

    /// Track parameters covariance matrix.
    Covariance
    predictedCovariance() const;

    bool
    hasPredicted() const
    {
      return m_data.ipredicted != IndexData::kInvalid;
    }

    /// Track parameters vector.
    Parameters
    filtered() const;

    /// Track parameters covariance matrix.
    Covariance
    filteredCovariance() const;

    bool
    hasFiltered() const
    {
      return m_data.ifiltered != IndexData::kInvalid;
    }

    /// Track parameters vector.
    Parameters
    smoothed() const;

    /// Track parameters covariance matrix.
    Covariance
    smoothedCovariance() const;

    bool
    hasSmoothed() const
    {
      return m_data.ismoothed != IndexData::kInvalid;
    }

    bool
    hasUncalibrated() const
    {
      return m_data.iuncalibrated != IndexData::kInvalid;
    }

    /// Uncalibrated measurement in the form of a source link
    const SourceLink&
    uncalibrated() const;

    /// Check if the point has an associated measurement.
    bool
    hasCalibrated() const
    {
      return m_data.icalibrated != IndexData::kInvalid;
    }

    /// Full measurement vector. Might contain additional zeroed dimensions.
    Measurement
    calibrated() const;

    /// Full measurement covariance matrix.
    MeasurementCovariance
    calibratedCovariance() const;

    /// Dynamic measurement vector with only the valid dimensions.
    auto
    effectiveCalibrated() const
    {
      return calibrated().head(m_data.measdim);
    }

    /// Dynamic measurement covariance matrix with only the valid dimensions.
    auto
    effectiveCalibratedCovariance() const
    {
      return calibratedCovariance().topLeftCorner(m_data.measdim,
                                                  m_data.measdim);
    }

  private:
    // Private since it can only be created by the trajectory.
    TrackStateProxy(ConstIf<MultiTrajectory<SourceLink>, ReadOnly>& trajectory,
                    size_t istate);

    ConstIf<MultiTrajectory<SourceLink>, ReadOnly>& m_traj;
    size_t    m_istate;
    IndexData m_data;

    friend class Acts::MultiTrajectory<SourceLink>;
  };
}  // namespace detail_lt

/// Store a trajectory of track states with multiple components.
///
/// This container supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over specific sub-components.
/// @tparam source_link_t Type to link back to an original measurement
template <typename source_link_t>
class MultiTrajectory
{
public:
  enum {
    ParametersSize     = NGlobalPars,
    MeasurementSizeMax = 2,
  };
  using SourceLink           = source_link_t;
  using ConstTrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  using TrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @param trackParameters  at the local point
  /// @param iprevious        index of the previous state, SIZE_MAX if first
  template <typename parameters_t>
  size_t
  addTrackState(const TrackState<SourceLink, parameters_t>& ts,
                size_t iprevious = SIZE_MAX);

  /// Access a read-only point on the trajectory by index.
  ConstTrackStateProxy
  getTrackState(size_t istate) const
  {
    return {*this, istate};
  }

  /// Access a writable point on the trajectory by index.
  TrackStateProxy
  getTrackState(size_t istate)
  {
    return {*this, istate};
  }

  /// Visit all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   non-modifying functor to be called with each point
  template <typename F>
  void
  visitBackwards(size_t iendpoint, F&& callable) const;
  /// Apply a function to all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F>
  void
  applyBackwards(size_t iendpoint, F&& callable);

private:
  /// index to map track states to the corresponding
  std::vector<detail_lt::IndexData>                                  m_index;
  typename detail_lt::Types<ParametersSize>::StorageCoefficients     m_params;
  typename detail_lt::Types<ParametersSize>::StorageCovariance       m_cov;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCoefficients m_meas;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCovariance   m_measCov;
  std::vector<SourceLink> m_sourceLinks;

  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;
};

// implementations

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
    return Parameters(m_traj.m_params.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::predictedCovariance() const -> Covariance
  {
    return Covariance(m_traj.m_cov.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filtered() const -> Parameters
  {
    return Parameters(m_traj.m_params.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filteredCovariance() const -> Covariance
  {
    return Covariance(m_traj.m_cov.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothed() const -> Parameters
  {
    return Parameters(m_traj.m_params.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothedCovariance() const -> Covariance
  {
    return Covariance(m_traj.m_cov.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::uncalibrated() const -> const SourceLink&
  {
    return m_traj.m_sourceLinks[m_data.iuncalibrated];
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibrated() const -> Measurement
  {
    return Measurement(m_traj.m_meas.col(m_data.icalibrated).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibratedCovariance() const
      -> MeasurementCovariance
  {
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

  // handle measurements
  if (ts.measurement.uncalibrated) {
    m_sourceLinks.push_back(*ts.measurement.uncalibrated);
    p.iuncalibrated = m_sourceLinks.size() - 1;
  }

  if (ts.measurement.calibrated) {
    auto meas    = m_meas.addCol();
    auto measCov = m_measCov.addCol();
    std::visit(
        [&meas, &measCov](const auto& m) {
          using meas_t                         = std::decay_t<decltype(m)>;
          meas.template head<meas_t::size()>() = m.parameters();
          measCov.template topLeftCorner<meas_t::size(), meas_t::size()>()
              = m.covariance();
        },
        *ts.measurement.calibrated);
    p.icalibrated = meas.size() - 1;
  }

  m_index.push_back(std::move(p));
  return m_index.size() - 1;
}

template <typename SL>
template <typename F>
void
MultiTrajectory<SL>::visitBackwards(size_t iendpoint, F&& callable) const
{
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
  while (true) {
    callable(getTrackState(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

}  // namespace Acts
