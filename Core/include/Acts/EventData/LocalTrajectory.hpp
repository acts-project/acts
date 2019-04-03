// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackParametersBase.hpp"

namespace Acts {
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
  /// Type construction helper for coefficients and associated covariances.
  template <Eigen::Index MaxSize, Eigen::Index SizeIncrement = 8>
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
    using FullCoefficientsConstMap = Eigen::Map<const FullCoefficients>;
    using FullCovarianceConstMap   = Eigen::Map<const FullCovariance>;

    // sub-vector, sub-matrix single items
    using Coefficients
        = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Flags, MaxSize, 1>;
    using Covariance = Eigen::
        Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Flags, MaxSize, MaxSize>;
    using CoefficientsConstMap = Eigen::Map<const Coefficients>;
    // TODO can we determine the pointer alignment at compile time?
    using CovarianceConstMap = Eigen::
        Map<const Covariance, Eigen::Unaligned, Eigen::Stride<MaxSize, 1>>;
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
}  // namespace detail_lt

class LocalTrajectory;

/// Proxy object to access a single point on the trajectory.
class LocalTrajectoryPoint
{
public:
  using FullParameters  = detail_lt::Types<8>::FullCoefficientsConstMap;
  using FullCovariance  = detail_lt::Types<8>::FullCovarianceConstMap;
  using FullMeasurement = detail_lt::Types<2>::FullCoefficientsConstMap;
  using FullMeasurementCovariance = detail_lt::Types<2>::FullCovarianceConstMap;
  using Parameters                = detail_lt::Types<8>::CoefficientsConstMap;
  using Covariance                = detail_lt::Types<8>::CovarianceConstMap;
  using Measurement               = detail_lt::Types<2>::CoefficientsConstMap;
  using MeasurementCovariance     = detail_lt::Types<2>::CovarianceConstMap;

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
  LocalTrajectoryPoint(const LocalTrajectory& trajectory, size_t index)
    : m_traj(trajectory), m_index(index)
  {
  }

  const LocalTrajectory& m_traj;
  size_t                 m_index;

  friend class LocalTrajectory;
};

/// Store local states, covariances, measurements along a trajectory.
///
/// The trajectory supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over sub-components.
class LocalTrajectory
{
public:
  /// Create an empty trajectory.
  LocalTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @param trackParameters  at the local point
  /// @param previous         point index or SIZE_MAX if its the first
  size_t
  addPoint(const TrackParametersBase& trackParameters,
           size_t                     previous = SIZE_MAX);
  /// Access a point on the trajectory by index.
  LocalTrajectoryPoint
  getPoint(size_t index) const
  {
    return {*this, index};
  }

  /// Visit all previous trajectory points starting from a given endpoint
  ///
  /// @param endPoint  index of the last track point in the trajectory
  /// @param visit     Functor to be called w/ each track point
  template <typename Visitor>
  void
  traverseBackward(size_t endpoint, Visitor visit) const;

private:
  std::vector<detail_lt::PointData>        m_points;
  detail_lt::Types<8>::StorageCoefficients m_params;
  detail_lt::Types<8>::StorageCovariance   m_cov;
  detail_lt::Types<2>::StorageCoefficients m_meas;
  detail_lt::Types<2>::StorageCovariance   m_measCov;

  friend class LocalTrajectoryPoint;
};

// implementations

inline LocalTrajectoryPoint::FullParameters
LocalTrajectoryPoint::fullParameters() const
{
  const auto& point = m_traj.m_points[m_index];
  return LocalTrajectoryPoint::FullParameters(
      m_traj.m_params.data.col(point.iparams).data());
}

inline LocalTrajectoryPoint::FullCovariance
LocalTrajectoryPoint::fullCovariance() const
{
  const auto& point = m_traj.m_points[m_index];
  return LocalTrajectoryPoint::FullCovariance(
      m_traj.m_cov.data.col(point.iparams).data());
}

inline LocalTrajectoryPoint::FullMeasurement
LocalTrajectoryPoint::fullMeasurement() const
{
  const auto& point = m_traj.m_points[m_index];
  return LocalTrajectoryPoint::FullMeasurement(
      m_traj.m_meas.data.col(point.imeas).data());
}

inline LocalTrajectoryPoint::FullMeasurementCovariance
LocalTrajectoryPoint::fullMeasurementCovariance() const
{
  const auto& point = m_traj.m_points[m_index];
  return LocalTrajectoryPoint::FullMeasurementCovariance(
      m_traj.m_measCov.data.col(point.imeas).data());
}

inline LocalTrajectoryPoint::Parameters
LocalTrajectoryPoint::parameters() const
{
  const auto& point = m_traj.m_points[m_index];
  return {m_traj.m_params.data.col(point.iparams).data(), point.nparams};
}

inline LocalTrajectoryPoint::Covariance
LocalTrajectoryPoint::covariance() const
{
  const auto& point = m_traj.m_points[m_index];
  return {m_traj.m_cov.data.col(point.iparams).data(),
          point.nparams,
          point.nparams};
}

inline LocalTrajectoryPoint::Measurement
LocalTrajectoryPoint::measurement() const
{
  const auto& point = m_traj.m_points[m_index];
  return {m_traj.m_meas.data.col(point.imeas).data(), point.nmeas};
}

inline LocalTrajectoryPoint::MeasurementCovariance
LocalTrajectoryPoint::measurementCovariance() const
{
  const auto& point = m_traj.m_points[m_index];
  return {
      m_traj.m_measCov.data.col(point.imeas).data(), point.nmeas, point.nmeas};
}

inline bool
LocalTrajectoryPoint::hasMeasurement() const
{
  return m_traj.m_points[m_index].imeas != detail_lt::PointData::kInvalid;
}

inline size_t
LocalTrajectory::addPoint(const TrackParametersBase& trackParameters,
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

template <typename Visitor>
void
LocalTrajectory::traverseBackward(size_t endpoint, Visitor visit) const
{
  if (m_points.size() <= endpoint) {
    // TODO to fail or not to fail here?
    return;
  }
  // TODO check input valididty
  while (true) {
    visit(LocalTrajectoryPoint(*this, endpoint));
    // this point has no parent and ends the trajectory
    if (m_points[endpoint].iprevious == detail_lt::PointData::kInvalid) {
      break;
    }
    endpoint = m_points[endpoint].iprevious;
  }
}

}  // namespace Acts
