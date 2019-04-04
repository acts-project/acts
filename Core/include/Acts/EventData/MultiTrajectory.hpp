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

namespace Acts {

// forward declarations
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
    Storage data;

    /// Access a column after ensuring the underlying storage is large enough.
    auto
    ensureCol(size_t index)
    {
      while (static_cast<size_t>(data.cols()) <= index) {
        data.conservativeResize(Eigen::NoChange, data.cols() + kSizeIncrement);
      }
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
    static constexpr uint16_t kInvalid = UINT16_MAX;

    const Surface& surface;
    uint16_t       iprevious     = kInvalid;
    uint16_t       iparams       = kInvalid;
    uint16_t       iuncalibrated = kInvalid;
    uint16_t       icalibrated   = kInvalid;
    uint16_t       measdim       = 0;
  };
  /// Proxy object to access a single point on the trajectory.
  ///
  /// @tparam ReadOnly  true for read-only access to underlying storage
  /// @tparam N         Number of track parameters
  /// @tparam M         Maximum number of measurements
  template <size_t N, size_t M, bool ReadOnly = true>
  class TrackStateProxy
  {
  public:
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
    bool
    hasUncalibrated() const
    {
      return m_data.iuncalibrated != IndexData::kInvalid;
    }

    /// Full measurement vector. Might contain additional zeroed dimensions.
    Measurement
    uncalibrated() const;

    /// Full measurement covariance matrix.
    MeasurementCovariance
    uncalibratedCovariance() const;

    /// Dynamic measurement vector with only the valid dimensions.
    auto
    effectiveUncalibrated() const
    {
      return uncalibrated().head(m_data.measdim);
    }

    /// Dynamic measurement covariance matrix with only the valid dimensions.
    auto
    effectiveUncalibratedCovariance() const
    {
      return uncalibratedCovariance().topLeftCorner(m_data.measdim,
                                                    m_data.measdim);
    }

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
    TrackStateProxy(ConstIf<MultiTrajectory, ReadOnly>& trajectory,
                    size_t istate);

    ConstIf<MultiTrajectory, ReadOnly>& m_traj;
    size_t    m_istate;
    IndexData m_data;

    friend class Acts::MultiTrajectory;
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
class MultiTrajectory
{
public:
  enum {
    ParametersSize     = 6,
    MeasurementSizeMax = 2,
  };
  using ConstTrackStateProxy
      = detail_lt::TrackStateProxy<ParametersSize, MeasurementSizeMax, true>;
  using TrackStateProxy
      = detail_lt::TrackStateProxy<ParametersSize, MeasurementSizeMax, false>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @param trackParameters  at the local point
  /// @param iprevious        index of the previous state, SIZE_MAX if first
  size_t
  addPoint(const TrackParametersBase& trackParameters,
           size_t                     iprevious = SIZE_MAX);
  /// Access a read-only point on the trajectory by index.
  ConstTrackStateProxy
  getPoint(size_t istate) const
  {
    return {*this, istate};
  }
  /// Access a writable point on the trajectory by index.
  TrackStateProxy
  getPoint(size_t istate)
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
  std::vector<detail_lt::IndexData>                         m_index;
  detail_lt::Types<ParametersSize>::StorageCoefficients     m_params;
  detail_lt::Types<ParametersSize>::StorageCovariance       m_cov;
  detail_lt::Types<MeasurementSizeMax>::StorageCoefficients m_meas;
  detail_lt::Types<MeasurementSizeMax>::StorageCovariance   m_measCov;

  friend class detail_lt::
      TrackStateProxy<ParametersSize, MeasurementSizeMax, true>;
  friend class detail_lt::
      TrackStateProxy<ParametersSize, MeasurementSizeMax, false>;
};

// implementations

namespace detail_lt {
  template <size_t N, size_t M, bool ReadOnly>
  inline TrackStateProxy<N, M, ReadOnly>::TrackStateProxy(
      ConstIf<MultiTrajectory, ReadOnly>& trajectory,
      size_t istate)
    : m_traj(trajectory), m_istate(istate), m_data(trajectory.m_index[istate])
  {
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::parameters() const -> Parameters
  {
    return Parameters(m_traj.m_params.data.col(m_data.iparams).data());
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::covariance() const -> Covariance
  {
    return Covariance(m_traj.m_cov.data.col(m_data.iparams).data());
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::uncalibrated() const -> Measurement
  {
    return Measurement(m_traj.m_meas.data.col(m_data.iuncalibrated).data());
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::uncalibratedCovariance() const
      -> MeasurementCovariance
  {
    return MeasurementCovariance(
        m_traj.m_measCov.data.col(m_data.iuncalibrated).data());
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::calibrated() const -> Measurement
  {
    return Measurement(m_traj.m_meas.data.col(m_data.icalibrated).data());
  }

  template <size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<N, M, ReadOnly>::calibratedCovariance() const
      -> MeasurementCovariance
  {
    return MeasurementCovariance(
        m_traj.m_measCov.data.col(m_data.icalibrated).data());
  }
}  // namespace detail_lt

inline size_t
MultiTrajectory::addPoint(const TrackParametersBase& trackParameters,
                          size_t                     iprevious)
{
  using Par    = TrackParametersBase::ParVector_t;
  using CovMap = detail_lt::Types<ParametersSize, false>::CovarianceMap;

  constexpr auto nparams = Par::RowsAtCompileTime;
  auto           iparams = m_index.size();

  m_params.ensureCol(iparams).setZero();
  m_params.col(iparams).head<nparams>() = trackParameters.parameters();
  m_cov.ensureCol(iparams).setZero();
  const auto* cov = trackParameters.covariance();
  if (cov != nullptr) {
    CovMap(m_cov.col(iparams).data()).topLeftCorner<nparams, nparams>() = *cov;
  }

  detail_lt::IndexData p = {trackParameters.referenceSurface()};
  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
  }
  p.iparams = iparams;
  m_index.push_back(std::move(p));

  return m_index.size() - 1;
}

template <typename F>
void
MultiTrajectory::visitBackwards(size_t iendpoint, F&& callable) const
{
  while (true) {
    callable(getPoint(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

template <typename F>
void
MultiTrajectory::applyBackwards(size_t iendpoint, F&& callable)
{
  while (true) {
    callable(getPoint(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

}  // namespace Acts
