// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"

#include <unordered_map>

namespace Acts {

class VectorMultiTrajectory final
    : public MultiTrajectory<VectorMultiTrajectory> {
#ifndef DOXYGEN
  friend MultiTrajectory<VectorMultiTrajectory>;
#endif

  struct IndexData {
    IndexType iprevious = kInvalid;
    IndexType ipredicted = kInvalid;
    IndexType ifiltered = kInvalid;
    IndexType ismoothed = kInvalid;
    IndexType ijacobian = kInvalid;
    IndexType iprojector = kInvalid;

    double chi2 = 0;
    double pathLength;
    TrackStateType typeFlags;

    IndexType iuncalibrated = kInvalid;
    IndexType icalibrated = kInvalid;
    IndexType icalibratedsourcelink = kInvalid;
    IndexType measdim = 0;
  };

 public:
  VectorMultiTrajectory() = default;

  VectorMultiTrajectory(const VectorMultiTrajectory& other)
      : m_index{other.m_index},
        m_previous{other.m_previous},
        m_params{other.m_params},
        m_cov{other.m_cov},
        m_meas{other.m_meas},
        m_measCov{other.m_measCov},
        m_jac{other.m_jac},
        m_sourceLinks{other.m_sourceLinks},
        m_projectors{other.m_projectors},
        m_referenceSurfaces{other.m_referenceSurfaces} {
    for (const auto& [key, value] : other.m_dynamic) {
      m_dynamic.insert({key, value->clone()});
    }
  };

  VectorMultiTrajectory(VectorMultiTrajectory&&) = default;
  VectorMultiTrajectory& operator=(const VectorMultiTrajectory&) = default;
  VectorMultiTrajectory& operator=(VectorMultiTrajectory&&) = default;

 private:
  // BEGIN INTERFACE
  TrackStateProxy::Parameters parameters_impl(IndexType parIdx) {
    return TrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Parameters parameters_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  TrackStateProxy::Covariance covariance_impl(IndexType parIdx) {
    return TrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  TrackStateProxy::Covariance jacobian_impl(IndexType parIdx) {
    return TrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  TrackStateProxy::Measurement measurement_impl(IndexType parIdx) {
    return TrackStateProxy::Measurement{m_meas[parIdx].data()};
  }

  ConstTrackStateProxy::Measurement measurement_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Measurement{m_meas[parIdx].data()};
  }

  TrackStateProxy::MeasurementCovariance measurementCovariance_impl(
      IndexType parIdx) {
    return TrackStateProxy::MeasurementCovariance{m_measCov[parIdx].data()};
  }

  ConstTrackStateProxy::MeasurementCovariance measurementCovariance_impl(
      IndexType parIdx) const {
    return ConstTrackStateProxy::MeasurementCovariance{
        m_measCov[parIdx].data()};
  }

  IndexType addTrackState_impl(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid);

  void shareFrom_impl(IndexType iself, IndexType iother,
                      TrackStatePropMask shareSource,
                      TrackStatePropMask shareTarget);

  void unset_impl(TrackStatePropMask target, IndexType istate);

  constexpr bool has_impl(HashedString key, IndexType istate) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
        return m_index[istate].ipredicted != kInvalid;
      case "filtered"_hash:
        return m_index[istate].ifiltered != kInvalid;
      case "smoothed"_hash:
        return m_index[istate].ismoothed != kInvalid;
      case "calibrated"_hash:
        return m_index[istate].icalibrated != kInvalid;
      case "jacobian"_hash:
        return m_index[istate].ijacobian != kInvalid;
      case "projector"_hash:
        return m_index[istate].iprojector != kInvalid;
      case "uncalibrated"_hash: {
        const auto& sl = m_sourceLinks[m_index[istate].iuncalibrated];
        return sl != nullptr;
      }
      case "previous"_hash:
      case "calibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return m_dynamic.find(key) != m_dynamic.end();
    }
  }

  IndexType size_impl() const { return m_index.size(); }

  void clear_impl();

  std::any component_impl(HashedString key, IndexType istate) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "previous"_hash:
        return &m_index[istate].iprevious;
      case "predicted"_hash:
        return &m_index[istate].ipredicted;
      case "filtered"_hash:
        return &m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &m_index[istate].ismoothed;
      case "calibrated"_hash:
        return &m_index[istate].icalibrated;
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "uncalibrated"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[istate];
      case "measdim"_hash:
        return &m_index[istate].measdim;
      case "chi2"_hash:
        return &m_index[istate].chi2;
      case "pathLength"_hash:
        return &m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &m_index[istate].typeFlags;
      default:
        auto it = m_dynamic.find(key);
        if (it == m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        auto& col = it->second;
        assert(col && "Dynamic column is null");
        return col->get(istate);
    }
  }

  std::any component_impl(HashedString key, IndexType istate) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "previous"_hash:
        return &m_index[istate].iprevious;
      case "predicted"_hash:
        return &m_index[istate].ipredicted;
      case "filtered"_hash:
        return &m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &m_index[istate].ismoothed;
      case "calibrated"_hash:
        return &m_index[istate].icalibrated;
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "uncalibrated"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[istate];
      case "measdim"_hash:
        return &m_index[istate].measdim;
      case "chi2"_hash:
        return &m_index[istate].chi2;
      case "pathLength"_hash:
        return &m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &m_index[istate].typeFlags;
      default:
        auto it = m_dynamic.find(key);
        if (it == m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        auto& col = it->second;
        assert(col && "Dynamic column is null");
        return col->get(istate);
    }
  }

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {
    m_dynamic.insert({hashString(key), std::make_unique<DynamicColumn<T>>()});
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
      case "previous"_hash:
      case "uncalibrated"_hash:
      case "calibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return m_dynamic.find(key) != m_dynamic.end();
    }
  }

  // END INTERFACE

  struct DynamicColumnBase {
    virtual ~DynamicColumnBase() = 0;

    virtual std::any get(size_t i) = 0;
    virtual std::any get(size_t i) const = 0;

    virtual void add() = 0;
    virtual void clear() = 0;

    virtual std::unique_ptr<DynamicColumnBase> clone() const = 0;
  };

  template <typename T>
  struct DynamicColumn : public VectorMultiTrajectory::DynamicColumnBase {
    ~DynamicColumn() override = default;

    std::any get(size_t i) override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    std::any get(size_t i) const override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    void add() override { m_vector.emplace_back(); }
    void clear() override { m_vector.clear(); }

    std::unique_ptr<DynamicColumnBase> clone() const override {
      return std::make_unique<DynamicColumn<T>>(*this);
    }

    std::vector<T> m_vector;
  };

  /// index to map track states to the corresponding
  std::vector<IndexData> m_index;
  std::vector<IndexType> m_previous;
  std::vector<typename detail_lt::Types<eBoundSize>::Coefficients> m_params;
  std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_cov;
  std::vector<typename detail_lt::Types<MeasurementSizeMax>::Coefficients>
      m_meas;
  std::vector<typename detail_lt::Types<MeasurementSizeMax>::Covariance>
      m_measCov;
  std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_jac;
  std::vector<const SourceLink*> m_sourceLinks;
  std::vector<ProjectorBitset> m_projectors;

  // owning vector of shared pointers to surfaces
  //
  // This might be problematic when appending a large number of surfaces
  // trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  std::unordered_map<HashedString, std::unique_ptr<DynamicColumnBase>>
      m_dynamic;
};  // namespace Acts

}  // namespace Acts
