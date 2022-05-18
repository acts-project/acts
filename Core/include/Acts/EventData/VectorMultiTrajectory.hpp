// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"

namespace Acts {

class VectorMultiTrajectory final : public MultiTrajectory {
 public:
  std::size_t previous(IndexType istate) const override {
    return m_index[istate].iprevious;
  }

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
                            size_t iprevious = SIZE_MAX) override {
    using PropMask = TrackStatePropMask;

    m_index.emplace_back();
    detail_lt::IndexData& p = m_index.back();
    size_t index = m_index.size() - 1;

    if (iprevious != SIZE_MAX) {
      p.iprevious = static_cast<uint16_t>(iprevious);
    }

    // always set, but can be null
    m_referenceSurfaces.emplace_back(nullptr);
    p.irefsurface = m_referenceSurfaces.size() - 1;

    assert(m_params.size() == m_cov.size());

    if (ACTS_CHECK_BIT(mask, PropMask::Predicted)) {
      m_params.emplace_back();
      m_cov.emplace_back();
      p.ipredicted = m_params.size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, PropMask::Filtered)) {
      m_params.emplace_back();
      m_cov.emplace_back();
      p.ifiltered = m_params.size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, PropMask::Smoothed)) {
      m_params.emplace_back();
      m_cov.emplace_back();
      p.ismoothed = m_params.size() - 1;
    }

    assert(m_params.size() == m_cov.size());

    if (ACTS_CHECK_BIT(mask, PropMask::Jacobian)) {
      m_jac.emplace_back();
      p.ijacobian = m_jac.size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, PropMask::Uncalibrated)) {
      m_sourceLinks.emplace_back();
      p.iuncalibrated = m_sourceLinks.size() - 1;
    }

    if (ACTS_CHECK_BIT(mask, PropMask::Calibrated)) {
      m_meas.emplace_back();
      m_measCov.emplace_back();
      p.icalibrated = m_meas.size() - 1;

      m_sourceLinks.emplace_back();
      p.icalibratedsourcelink = m_sourceLinks.size() - 1;

      m_projectors.emplace_back();
      p.iprojector = m_projectors.size() - 1;
    }

    return index;
  }

  std::size_t size() const override { return m_index.size(); }

  void clear() override {
    std::cout << "traj: " << this << " clear" << std::endl;
    m_index.clear();
    m_params.clear();
    m_cov.clear();
    m_meas.clear();
    m_measCov.clear();
    m_jac.clear();
    m_sourceLinks.clear();
    m_projectors.clear();
    m_referenceSurfaces.clear();
  }

 protected:
  void* componentImpl(HashedString key, IndexType istate) override {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
        return &m_index[istate].ipredicted;
      case "filtered"_hash:
        return &m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &m_index[istate].ismoothed;
      case "measurement"_hash:
        return &m_index[istate].icalibrated;
      // case "measurementCovariance"_hash:
      // return &m_index[istate].icalibrated;
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "sourceLink"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[m_index[istate].irefsurface];
      default:
        assert(false);
    }
  }

  const void* componentImpl(HashedString key, IndexType istate) const override {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
        return &m_index[istate].ipredicted;
      case "filtered"_hash:
        return &m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &m_index[istate].ismoothed;
      case "measurement"_hash:
        return &m_index[istate].icalibrated;
      // case "measurementCovariance"_hash:
      // return &m_measCov[m_index[istate].icalibrated];
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "sourceLink"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[m_index[istate].irefsurface];
      default:
        assert(false);
    }
  }

  const detail_lt::IndexData& data(IndexType istate) const override {
    return m_index[istate];
  }

  detail_lt::IndexData& data(IndexType istate) override {
    return m_index[istate];
  }

 private:
  /// index to map track states to the corresponding
  std::vector<detail_lt::IndexData> m_index;
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
};

}  // namespace Acts
