// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

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

      // @TODO: do this or not? If we assume this happens only when something is
      // written, the expectation is that everything is zero
      data.col(index).setZero();

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

    IndexType irefsurface = kInvalid;
    IndexType iprevious   = kInvalid;
    IndexType ipredicted  = kInvalid;
    IndexType ifiltered   = kInvalid;
    IndexType ismoothed   = kInvalid;
    IndexType ijacobian   = kInvalid;
    IndexType iprojector  = kInvalid;

    IndexType iuncalibrated         = kInvalid;
    IndexType icalibrated           = kInvalid;
    IndexType icalibratedsourcelink = kInvalid;
    IndexType measdim               = 0;
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

    // as opposed to the types above, this is an actual Matrix (rather than a
    // map)
    // @TODO: Does not copy flags, because this fails: can't have col major row
    // vector, but that's required for 1xN projection matrices below.
    constexpr static auto ProjectorFlags = Eigen::RowMajor | Eigen::AutoAlign;
    using Projector
        = Eigen::Matrix<typename Covariance::Scalar, M, N, ProjectorFlags>;
    using EffectiveProjector = Eigen::Matrix<typename Projector::Scalar,
                                             Eigen::Dynamic,
                                             Eigen::Dynamic,
                                             ProjectorFlags,
                                             M,
                                             N>;

    /// Index within the trajectory.
    /// @return the index
    size_t
    index() const
    {
      return m_istate;
    }

    /// Reference surface.
    /// @return the reference surface
    const Surface&
    referenceSurface() const
    {
      assert(m_data.irefsurface != IndexData::kInvalid);
      return *m_traj.m_referenceSurfaces[m_data.irefsurface];
    }

    /// Track parameters vector. This tries to be somewhat smart and return the
    /// first parameters that are set in this order: predicted -> filtered ->
    /// smoothed
    /// @return one of predicted, filtered or smoothed parameters
    Parameters
    parameters() const;

    /// Track parameters covariance matrix. This tries to be somewhat smart and
    /// return the
    /// first parameters that are set in this order: predicted -> filtered ->
    /// smoothed
    /// @return one of predicted, filtered or smoothed covariances
    Covariance
    covariance() const;

    /// Predicted track parameters vector
    /// @return The predicted parameters
    Parameters
    predicted() const;

    /// Predicted track parameters covariance matrix.
    /// @return The predicted track parameter covariance
    Covariance
    predictedCovariance() const;

    /// Check whether the predicted parameters+covariance is set
    /// @return Whether it is set or not
    bool
    hasPredicted() const
    {
      return m_data.ipredicted != IndexData::kInvalid;
    }

    /// Filtered track parameters vector
    /// @return The filtered parameters
    Parameters
    filtered() const;

    /// Filtered track parameters covariance matrix
    /// @return The filtered parameters covariance
    Covariance
    filteredCovariance() const;

    /// Return whether filtered parameters+covariance is set
    /// @return Whether it is set
    bool
    hasFiltered() const
    {
      return m_data.ifiltered != IndexData::kInvalid;
    }

    /// Smoothed track parameters vector
    /// @return the parameter vector
    Parameters
    smoothed() const;

    /// Smoothed track parameters covariance matrix
    /// @return the parameter covariance matrix
    Covariance
    smoothedCovariance() const;

    /// Return whether smoothed parameters+covariance is set
    /// @return Whether it is set
    bool
    hasSmoothed() const
    {
      return m_data.ismoothed != IndexData::kInvalid;
    }

    /// Returns the jacobian from the previous trackstate to this one
    /// @return The jacobian matrix
    Covariance
    jacobian() const;

    /// Returns whether a jacobian is set for this trackstate
    /// @return Whether it is set
    bool
    hasJacobian() const
    {
      return m_data.ijacobian != IndexData::kInvalid;
    }

    /// Returns the projector (measurement mapping function) for this track
    /// state. It is derived from the uncalibrated measurement
    /// @note This function returns the overallocated projector. This means it
    /// is of dimension MxM, where M is the maximum number of measurement
    /// dimensions. The NxM submatrix, where N is the actual dimension of the
    /// measurement, is located in the top left corner, everything else is zero.
    /// @return The overallocated projector
    Projector
    projector() const;

    /// Returns whether a projector is set
    /// @return Whether it is set
    bool
    hasProjector() const
    {
      return m_data.iprojector != IndexData::kInvalid;
    }

    /// Returns the projector (measurement mapping function) for this track
    /// state. It is derived from the uncalibrated measurement
    /// @note This function returns the effective projector. This means it
    /// is of dimension NxM, where N is the actual dimension of the
    /// measurement.
    /// @return The effective projector
    EffectiveProjector
    effectiveProjector() const
    {
      return projector().topLeftCorner(m_data.measdim, M);
    }

    /// Return whether an uncalibrated measurement (source link) is set
    /// @return Whether it is set
    bool
    hasUncalibrated() const
    {
      return m_data.iuncalibrated != IndexData::kInvalid;
    }

    /// Uncalibrated measurement in the form of a source link. Const version
    /// @return The uncalibrated measurement source link
    const SourceLink&
    uncalibrated() const;

    /// Uncalibrated measurement in the form of a source link. Mutable version
    /// @note This overload is only available if @c ReadOnly is false
    /// @return The uncalibrated measurement source link
    template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
    SourceLink&
    uncalibrated()
    {
      assert(m_data.iuncalibrated != IndexData::kInvalid);
      return m_traj.m_sourceLinks[m_data.iuncalibrated];
    }

    /// Check if the point has an associated calibrated measurement.
    /// @return Whether it is set
    bool
    hasCalibrated() const
    {
      return m_data.icalibrated != IndexData::kInvalid;
    }

    /// The source link of the calibrated measurement. Const version
    /// @note This does not necessarily have to be the uncalibrated source link.
    /// @return The source link
    const SourceLink&
    calibratedSourceLink() const;

    /// The source link of the calibrated measurement. Mutable version
    /// @note This does not necessarily have to be the uncalibrated source link.
    /// @note This overload is only available if @c ReadOnly is false
    /// @return The source link
    template <bool RO = ReadOnly, typename = std::enable_if_t<!RO>>
    SourceLink&
    calibratedSourceLink()
    {
      assert(m_data.icalibratedsourcelink != IndexData::kInvalid);
      return m_traj.m_sourceLinks[m_data.icalibratedsourcelink];
    }

    /// Full calibrated measurement vector. Might contain additional zeroed
    /// dimensions.
    /// @return The measurement vector
    Measurement
    calibrated() const;

    /// Full calibrated measurement covariance matrix. The effective covariance
    /// is located in the top left corner, everything else is zeroed.
    /// @return The measurement covariance matrix
    MeasurementCovariance
    calibratedCovariance() const;

    /// Dynamic measurement vector with only the valid dimensions.
    /// @return The effective calibrated measurement vector
    auto
    effectiveCalibrated() const
    {
      return calibrated().head(m_data.measdim);
    }

    /// Dynamic measurement covariance matrix with only the valid dimensions.
    /// @return The effective calibrated covariance matrix
    auto
    effectiveCalibratedCovariance() const
    {
      return calibratedCovariance().topLeftCorner(m_data.measdim,
                                                  m_data.measdim);
    }

    /// Return the (dynamic) number of dimensions stored for this measurement.
    /// @note The underlying storage is overallocated to MeasurementSizeMax
    /// regardless of this value
    /// @return The number of dimensions
    size_t
    calibratedSize() const
    {
      return m_data.measdim;
    }

    template <bool RO  = ReadOnly,
              typename = std::enable_if_t<!RO>,
              ParID_t... params>
    void
    setCalibrated(const Acts::Measurement<SourceLink, params...>& meas)
    {
      constexpr size_t measdim
          = Acts::Measurement<SourceLink, params...>::size();

      m_data.measdim = measdim;

      calibrated().setZero();
      calibrated().template head<measdim>() = meas.parameters();

      calibratedCovariance().setZero();
      calibratedCovariance().template topLeftCorner<measdim, measdim>()
          = meas.covariance();

      typename TrackStateProxy::Projector fullProjector;
      fullProjector.setZero();
      fullProjector.template topLeftCorner<measdim,
                                           MultiTrajectory<SourceLink>::
                                               MeasurementSizeMax>()
          = meas.projector();
      m_traj.m_projectors[m_data.iprojector] = matrixToBitset(fullProjector);

      calibratedSourceLink() = meas.sourceLink();
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

  // implement track state visitor concept
  template <typename T, typename TS>
  using call_operator_t = decltype(std::declval<T>()(std::declval<TS&>()));

  template <typename T, typename TS>
  constexpr bool VisitorConcept
      = concept::require<concept::exists<call_operator_t, T, TS>>;

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
    MeasurementSizeMax = NGlobalPars,
  };
  using SourceLink           = source_link_t;
  using ConstTrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  using TrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;

  using ProjectorBitset = std::bitset<ParametersSize * MeasurementSizeMax>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @tparam parameters_t The parameter type used for the trackstate
  /// @param trackParameters  at the local point
  /// @param iprevious        index of the previous state, SIZE_MAX if first
  /// @note The parameter type from @p parameters_t is not currently stored in
  /// MultiTrajectory.
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
  std::vector<SourceLink>      m_sourceLinks;
  std::vector<ProjectorBitset> m_projectors;
  // owning vector of shared pointers to surfaces
  // @TODO: This might be problematic when appending a large number of surfaces
  // / trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;
};

}  // namespace Acts

#include "Acts/EventData/MultiTrajectory.ipp"
