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
  ///
  /// All types have an upper size bound but the
  template <Eigen::Index MaxSize, Eigen::Index SizeIncrement = 8>
  struct Types
  {
    static constexpr int Flags = Eigen::ColMajor | Eigen::AutoAlign;
    using Scalar               = double;
    using StorageCoefficients
        = GrowableColumns<Eigen::Array<Scalar, MaxSize, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
    using StorageCovariance = GrowableColumns<
        Eigen::Array<Scalar, MaxSize * MaxSize, Eigen::Dynamic, Flags>,
        SizeIncrement>;
    using FullCoefficients = Eigen::Matrix<Scalar, MaxSize, 1, Flags>;
    using FullCovariance   = Eigen::Matrix<Scalar, MaxSize, MaxSize, Flags>;
    using Coefficients
        = Eigen::Matrix<Scalar, Eigen::Dynamic, 1, Flags, MaxSize, 1>;
    using Covariance = Eigen::
        Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Flags, MaxSize, MaxSize>;
    using CovarianceStride = Eigen::Stride<MaxSize, 1>;
  };
  struct PointData
  {
    static constexpr uint16_t kInvalid = UINT16_MAX;

    uint16_t iparams = kInvalid;
    uint16_t nparams = 0;
    uint16_t imeas   = kInvalid;
    uint16_t nmeas   = 0;
  };
}  // namespace detail_lt

class LocalTrajectory;

class LocalTrajectoryPoint
{
public:
  // underlying storage types w/ overallocation
  using FullParameters        = detail_lt::Types<8>::FullCoefficients;
  using FullCovariance        = detail_lt::Types<8>::FullCovariance;
  using Parameters            = detail_lt::Types<8>::Coefficients;
  using Covariance            = detail_lt::Types<8>::Covariance;
  using CovarianceStride      = detail_lt::Types<8>::CovarianceStride;
  using Measurement           = detail_lt::Types<2>::Coefficients;
  using MeasurementCovariance = detail_lt::Types<2>::Covariance;

  LocalTrajectoryPoint(const LocalTrajectory& trajectory, size_t index);

  Eigen::Map<const FullParameters>
  fullParameters() const;
  Eigen::Map<const FullCovariance>
  fullCovariance() const;
  Eigen::Map<const Parameters>
  parameters() const;
  Eigen::Map<const Covariance, Eigen::Unaligned, CovarianceStride>
  covariance() const;

  /// Check if the point has an associated measurement.
  bool
  hasMeasurement() const;

private:
  const LocalTrajectory& m_traj;
  detail_lt::PointData   m_point;
};

/// @brief Store local states, covariances, measurements along a trajectory
class LocalTrajectory
{
public:
  /// Create an empty trajectory.
  LocalTrajectory() = default;

  /// Add a point without measurement and return its index.
  size_t
  addPoint(const TrackParametersBase& trackParameters);
  /// Access a point on the trajectory by index.
  LocalTrajectoryPoint
  getPoint(size_t index) const
  {
    return {*this, index};
  }

private:
  std::vector<detail_lt::PointData>        m_points;
  detail_lt::Types<8>::StorageCoefficients m_params;
  detail_lt::Types<8>::StorageCovariance   m_cov;
  detail_lt::Types<2>::StorageCoefficients m_meas;
  detail_lt::Types<2>::StorageCovariance   m_measCov;

  friend class LocalTrajectoryPoint;
};

// implementations

inline LocalTrajectoryPoint::LocalTrajectoryPoint(
    const LocalTrajectory& trajectory,
    size_t                 index)
  : m_traj(trajectory), m_point(trajectory.m_points.at(index))
{
}

inline Eigen::Map<const LocalTrajectoryPoint::FullParameters>
LocalTrajectoryPoint::fullParameters() const
{
  return Eigen::Map<const LocalTrajectoryPoint::FullParameters>(
      m_traj.m_params.data.col(m_point.iparams).data());
}

inline Eigen::Map<const LocalTrajectoryPoint::FullCovariance>
LocalTrajectoryPoint::fullCovariance() const
{
  return Eigen::Map<const LocalTrajectoryPoint::FullCovariance>(
      m_traj.m_cov.data.col(m_point.iparams).data());
}

inline Eigen::Map<const LocalTrajectoryPoint::Parameters>
LocalTrajectoryPoint::parameters() const
{
  return {m_traj.m_params.data.col(m_point.iparams).data(), m_point.nparams};
}

inline Eigen::Map<const LocalTrajectoryPoint::Covariance,
                  Eigen::Unaligned,
                  LocalTrajectoryPoint::CovarianceStride>
LocalTrajectoryPoint::covariance() const
{
  return {m_traj.m_cov.data.col(m_point.iparams).data(),
          m_point.nparams,
          m_point.nparams};
}

inline bool
LocalTrajectoryPoint::hasMeasurement() const
{
  return m_point.imeas != decltype(m_point)::kInvalid;
}

inline size_t
LocalTrajectory::addPoint(const TrackParametersBase& trackParameters)
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
  p.iparams = iparams;
  p.nparams = nparams;
  m_points.push_back(std::move(p));

  return m_points.size() - 1;
}

}  // namespace Acts
