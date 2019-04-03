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

class MultiTrajectory;

namespace detail_lt {
  /// wrapper for a dynamic Eigen type that adds support for automatic growth
  ///
  /// \warning Assumes the underlying storage has a fixed number of rows
  template <typename Storage, Eigen::Index kSizeIncrement>
  struct GrowableColumns
  {
    Storage data;

    /// Access a column after ensuring the underlying storage is large enough.
    auto
    ensureCol(Eigen::Index index)
    {
      while (data.cols() <= index) {
        data.conservativeResize(Eigen::NoChange, data.cols() + kSizeIncrement);
      }
      return data.col(index);
    }
    /// Writable access to a column w/o checking its existence first.
    auto
    col(Eigen::Index index)
    {
      return data.col(index);
    }
    /// Read-only access to a column w/o checking its existence first.
    auto
    col(Eigen::Index index) const
    {
      return data.col(index);
    }
  };
  /// Either type T or const T depending on the boolean.
  template <typename T, bool select>
  using ConstIf = std::conditional_t<select, const T, T>;
  /// Type construction helper for coefficients and associated covariances.
  template <Eigen::Index MaxSize,
            bool         ReadOnlyMaps  = true,
            Eigen::Index SizeIncrement = 8>
  struct Types
  {
    static constexpr int Flags = Eigen::ColMajor | Eigen::AutoAlign;
    using Scalar               = double;

    // storage of multiple items in flat arrays (with overallocation up to max)
    using StorageCoefficients
        = GrowableColumns<Eigen::Array<Scalar, MaxSize, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
    using StorageCovariance = GrowableColumns<
        Eigen::Array<Scalar, MaxSize * MaxSize, Eigen::Dynamic, Flags>,
        SizeIncrement>;

    // full single items
    using FullCoefficients = Eigen::Matrix<Scalar, MaxSize, 1, Flags>;
    using FullCovariance   = Eigen::Matrix<Scalar, MaxSize, MaxSize, Flags>;
    using FullCoefficientsMap
        = Eigen::Map<ConstIf<FullCoefficients, ReadOnlyMaps>>;
    using FullCovarianceMap = Eigen::Map<ConstIf<FullCovariance, ReadOnlyMaps>>;

    // sub-vector, sub-matrix single items
    using Coefficients
        = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Flags, MaxSize, 1>;
    using Covariance = Eigen::
        Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Flags, MaxSize, MaxSize>;
    using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
    // TODO can we determine the pointer alignment at compile time?
    using CovarianceMap = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>,
                                     Eigen::Unaligned,
                                     Eigen::Stride<MaxSize, 1>>;
  };
  struct PointData
  {
    static constexpr uint16_t kInvalid = UINT16_MAX;

    uint16_t iprevious = kInvalid;
    uint16_t iparams   = kInvalid;
    uint16_t imeas     = kInvalid;
    uint8_t  nparams   = 0;
    uint8_t  nmeas     = 0;
  };
  /// Proxy object to access a single point on the trajectory.
  template <bool ReadOnly>
  class TrackStateProxy
  {
  public:
    using FullParameters  = typename Types<8, ReadOnly>::FullCoefficientsMap;
    using FullCovariance  = typename Types<8, ReadOnly>::FullCovarianceMap;
    using FullMeasurement = typename Types<2, ReadOnly>::FullCoefficientsMap;
    using FullMeasurementCovariance =
        typename Types<2, ReadOnly>::FullCovarianceMap;
    using Parameters            = typename Types<8, ReadOnly>::CoefficientsMap;
    using Covariance            = typename Types<8, ReadOnly>::CovarianceMap;
    using Measurement           = typename Types<2, ReadOnly>::CoefficientsMap;
    using MeasurementCovariance = typename Types<2, ReadOnly>::CovarianceMap;

    size_t
    index() const
    {
      return m_index;
    }

    FullParameters
    fullParameters() const;
    FullCovariance
    fullCovariance() const;
    FullMeasurement
    fullMeasurement() const;
    FullMeasurementCovariance
    fullMeasurementCovariance() const;

    Parameters
    parameters() const;
    Covariance
    covariance() const;
    Measurement
    measurement() const;
    MeasurementCovariance
    measurementCovariance() const;

    /// Check if the point has an associated measurement.
    bool
    hasMeasurement() const;

  private:
    // Private since it can only be created by the trajectory.
    TrackStateProxy(ConstIf<MultiTrajectory, ReadOnly>& trajectory,
                    size_t                              index)
      : m_traj(trajectory), m_index(index)
    {
    }

    ConstIf<MultiTrajectory, ReadOnly>& m_traj;
    size_t                              m_index;

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
  using ConstTrackStateProxy = detail_lt::TrackStateProxy<true>;
  using TrackStateProxy      = detail_lt::TrackStateProxy<false>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @param trackParameters  at the local point
  /// @param previous         point index or SIZE_MAX if its the first
  size_t
  addPoint(const TrackParametersBase& trackParameters,
           size_t                     previous = SIZE_MAX);
  /// Access a read-only point on the trajectory by index.
  ConstTrackStateProxy
  getPoint(size_t index) const
  {
    return {*this, index};
  }
  /// Access a writable point on the trajectory by index.
  TrackStateProxy
  getPoint(size_t index)
  {
    return {*this, index};
  }

  /// Visit all previous points starting at a given endpoint.
  ///
  /// @param endPoint  index of the last track point
  /// @param callable  non-modifying functor to be called with each point
  template <typename F>
  void
  visitBackwards(size_t endpoint, F&& callable) const;
  /// Apply a function to all previous points starting at a given endpoint.
  ///
  /// @param endPoint  index of the last track point
  /// @param callable  modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F>
  void
  applyBackwards(size_t endpoint, F&& callable);

private:
  std::vector<detail_lt::PointData>        m_points;
  detail_lt::Types<8>::StorageCoefficients m_params;
  detail_lt::Types<8>::StorageCovariance   m_cov;
  detail_lt::Types<2>::StorageCoefficients m_meas;
  detail_lt::Types<2>::StorageCovariance   m_measCov;

  friend class detail_lt::TrackStateProxy<true>;
  friend class detail_lt::TrackStateProxy<false>;
};

// implementations

namespace detail_lt {
  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::FullParameters
  TrackStateProxy<ReadOnly>::fullParameters() const
  {
    const auto& point = m_traj.m_points[m_index];
    return TrackStateProxy<ReadOnly>::FullParameters(
        m_traj.m_params.data.col(point.iparams).data());
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::FullCovariance
  TrackStateProxy<ReadOnly>::fullCovariance() const
  {
    const auto& point = m_traj.m_points[m_index];
    return TrackStateProxy<ReadOnly>::FullCovariance(
        m_traj.m_cov.data.col(point.iparams).data());
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::FullMeasurement
  TrackStateProxy<ReadOnly>::fullMeasurement() const
  {
    const auto& point = m_traj.m_points[m_index];
    return TrackStateProxy<ReadOnly>::FullMeasurement(
        m_traj.m_meas.data.col(point.imeas).data());
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::FullMeasurementCovariance
  TrackStateProxy<ReadOnly>::fullMeasurementCovariance() const
  {
    const auto& point = m_traj.m_points[m_index];
    return TrackStateProxy<ReadOnly>::FullMeasurementCovariance(
        m_traj.m_measCov.data.col(point.imeas).data());
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::Parameters
  TrackStateProxy<ReadOnly>::parameters() const
  {
    const auto& point = m_traj.m_points[m_index];
    return {m_traj.m_params.data.col(point.iparams).data(), point.nparams};
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::Covariance
  TrackStateProxy<ReadOnly>::covariance() const
  {
    const auto& point = m_traj.m_points[m_index];
    return {m_traj.m_cov.data.col(point.iparams).data(),
            point.nparams,
            point.nparams};
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::Measurement
  TrackStateProxy<ReadOnly>::measurement() const
  {
    const auto& point = m_traj.m_points[m_index];
    return {m_traj.m_meas.data.col(point.imeas).data(), point.nmeas};
  }

  template <bool ReadOnly>
  inline typename TrackStateProxy<ReadOnly>::MeasurementCovariance
  TrackStateProxy<ReadOnly>::measurementCovariance() const
  {
    const auto& point = m_traj.m_points[m_index];
    return {m_traj.m_measCov.data.col(point.imeas).data(),
            point.nmeas,
            point.nmeas};
  }

  template <bool ReadOnly>
  inline bool
  TrackStateProxy<ReadOnly>::hasMeasurement() const
  {
    return m_traj.m_points[m_index].imeas != detail_lt::PointData::kInvalid;
  }
}  // namespace detail_lt

inline size_t
MultiTrajectory::addPoint(const TrackParametersBase& trackParameters,
                          size_t                     previous)
{
  using Par        = TrackParametersBase::ParVector_t;
  using FullCov    = detail_lt::Types<8>::FullCovariance;
  using FullCovMap = Eigen::Map<FullCov>;

  constexpr auto nparams = Par::RowsAtCompileTime;
  auto           iparams = m_points.size();

  m_params.ensureCol(iparams).setZero();
  m_params.col(iparams).head<nparams>() = trackParameters.parameters();
  m_cov.ensureCol(iparams).setZero();
  const auto* cov = trackParameters.covariance();
  if (cov != nullptr) {
    FullCovMap(m_cov.col(iparams).data()).topLeftCorner<nparams, nparams>()
        = *cov;
  }

  detail_lt::PointData p;
  if (previous != SIZE_MAX) { p.iprevious = static_cast<uint16_t>(previous); }
  p.iparams = iparams;
  p.nparams = nparams;

  m_points.push_back(std::move(p));

  return m_points.size() - 1;
}

template <typename F>
void
MultiTrajectory::visitBackwards(size_t endpoint, F&& callable) const
{
  while (true) {
    callable(getPoint(endpoint));
    // this point has no parent and ends the trajectory
    if (m_points[endpoint].iprevious == detail_lt::PointData::kInvalid) {
      break;
    }
    endpoint = m_points[endpoint].iprevious;
  }
}

template <typename F>
void
MultiTrajectory::applyBackwards(size_t endpoint, F&& callable)
{
  while (true) {
    callable(getPoint(endpoint));
    // this point has no parent and ends the trajectory
    if (m_points[endpoint].iprevious == detail_lt::PointData::kInvalid) {
      break;
    }
    endpoint = m_points[endpoint].iprevious;
  }
}

}  // namespace Acts
