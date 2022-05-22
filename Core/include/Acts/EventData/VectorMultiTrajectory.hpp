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
  struct IndexData {
    size_t iprevious = kInvalid;
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
  TrackStateProxy::Parameters parameters(IndexType parIdx) override {
    return TrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Parameters parameters(IndexType parIdx) const override {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  TrackStateProxy::Covariance covariance(IndexType parIdx) override {
    return TrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance(IndexType parIdx) const override {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  TrackStateProxy::Covariance jacobian(IndexType parIdx) override {
    return TrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian(IndexType parIdx) const override {
    return ConstTrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  TrackStateProxy::Measurement measurement(IndexType parIdx) override {
    return TrackStateProxy::Measurement{m_meas[parIdx].data()};
  }

  ConstTrackStateProxy::Measurement measurement(
      IndexType parIdx) const override {
    return ConstTrackStateProxy::Measurement{m_meas[parIdx].data()};
  }

  TrackStateProxy::MeasurementCovariance measurementCovariance(
      IndexType parIdx) override {
    return TrackStateProxy::MeasurementCovariance{m_measCov[parIdx].data()};
  }

  ConstTrackStateProxy::MeasurementCovariance measurementCovariance(
      IndexType parIdx) const override {
    return ConstTrackStateProxy::MeasurementCovariance{
        m_measCov[parIdx].data()};
  }

  std::size_t addTrackState(TrackStatePropMask mask = TrackStatePropMask::All,
                            size_t iprevious = kNoPrevious) override;

  void shareFrom(IndexType iself, IndexType iother,
                 TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) override;

  void unset(TrackStatePropMask target, IndexType istate) override;

  bool has(HashedString key, IndexType istate) const override;

  std::size_t size() const override { return m_index.size(); }

  void clear() override;

  struct DynamicColumnBase {
    virtual ~DynamicColumnBase() = 0;

    virtual void* get(size_t i) = 0;
    virtual const void* get(size_t i) const = 0;

    virtual void add() = 0;
    virtual void clear() = 0;
  };

  template <typename T>
  struct DynamicColumn : public VectorMultiTrajectory::DynamicColumnBase {
    ~DynamicColumn() override = default;

    void* get(size_t i) override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    const void* get(size_t i) const override {
      assert(i < m_vector.size() && "DynamicColumn out of bounds");
      return &m_vector[i];
    }

    void add() override { m_vector.emplace_back(); }
    void clear() override { m_vector.clear(); }

    std::vector<T> m_vector;
  };

 protected:
  void* componentImpl(HashedString key, IndexType istate) override;

  const void* componentImpl(HashedString key, IndexType istate) const override;

  template <typename T>
  constexpr void addColumnImpl(HashedString key) {
    m_dynamic.insert({key, std::make_unique<DynamicColumn<T>>()});
  }

  constexpr bool hasColumnImpl(HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
      case "previous"_hash:
      case "sourceLink"_hash:
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

 private:
  /// index to map track states to the corresponding
  std::vector<IndexData> m_index;
  std::vector<size_t> m_previous;
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
  // @TODO: This might be problematic when appending a large number of surfaces
  // trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  std::unordered_map<HashedString, std::unique_ptr<DynamicColumnBase>>
      m_dynamic;
};  // namespace Acts

}  // namespace Acts
