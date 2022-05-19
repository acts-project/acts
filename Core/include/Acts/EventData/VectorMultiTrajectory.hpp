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

    m_sourceLinks.emplace_back();
    p.iuncalibrated = m_sourceLinks.size() - 1;

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

  void shareFrom(IndexType iself, IndexType iother,
                 TrackStatePropMask shareSource,
                 TrackStatePropMask shareTarget) override {
    // auto other = getTrackState(iother);
    detail_lt::IndexData& self = m_index[iself];
    detail_lt::IndexData& other = m_index[iother];

    assert(ACTS_CHECK_BIT(getTrackState(iother).getMask(), shareSource) &&
           "Source has incompatible allocation");

    // @TODO: Push behind interface somehow
    using PM = TrackStatePropMask;

    IndexType sourceIndex{kInvalid};
    switch (shareSource) {
      case PM::Predicted:
        sourceIndex = other.ipredicted;
        break;
      case PM::Filtered:
        sourceIndex = other.ifiltered;
        break;
      case PM::Smoothed:
        sourceIndex = other.ismoothed;
        break;
      case PM::Jacobian:
        sourceIndex = other.ijacobian;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }

    assert(sourceIndex != kInvalid);

    switch (shareTarget) {
      case PM::Predicted:
        assert(shareSource != PM::Jacobian);
        self.ipredicted = sourceIndex;
        break;
      case PM::Filtered:
        assert(shareSource != PM::Jacobian);
        self.ifiltered = sourceIndex;
        break;
      case PM::Smoothed:
        assert(shareSource != PM::Jacobian);
        self.ismoothed = sourceIndex;
        break;
      case PM::Jacobian:
        assert(shareSource == PM::Jacobian);
        self.ijacobian = sourceIndex;
        break;
      default:
        throw std::domain_error{"Unable to share this component"};
    }
  }

  void unset(TrackStatePropMask target, IndexType istate) override {
    using PM = TrackStatePropMask;

    switch (target) {
      case PM::Predicted:
        m_index[istate].ipredicted = kInvalid;
        break;
      case PM::Filtered:
        m_index[istate].ifiltered = kInvalid;
        break;
      case PM::Smoothed:
        m_index[istate].ismoothed = kInvalid;
        break;
      case PM::Jacobian:
        m_index[istate].ijacobian = kInvalid;
        break;
      case PM::Calibrated:
        m_index[istate].icalibrated = kInvalid;
        break;
      default:
        throw std::domain_error{"Unable to unset this component"};
    }
  }

  bool has(HashedString key, IndexType istate) const override {
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
      case "sourceLink"_hash:
      case "calibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        assert(false && "Unable to handle this component");
    }
  }

  std::size_t size() const override { return m_index.size(); }

  void clear() override {
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
      case "calibrated"_hash:
        return &m_index[istate].icalibrated;
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "sourceLink"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[m_index[istate].irefsurface];
      case "measdim"_hash:
        return &m_index[istate].measdim;
      case "chi2"_hash:
        return &m_index[istate].chi2;
      case "pathLength"_hash:
        return &m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &m_index[istate].typeFlags;
      default:
        assert(false && "Unable to handle this component");
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
      case "calibrated"_hash:
        return &m_index[istate].icalibrated;
      case "jacobian"_hash:
        return &m_index[istate].ijacobian;
      case "projector"_hash:
        return &m_projectors[m_index[istate].iprojector];
      case "sourceLink"_hash:
        return &m_sourceLinks[m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &m_sourceLinks[m_index[istate].icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &m_referenceSurfaces[m_index[istate].irefsurface];
      case "measdim"_hash:
        return &m_index[istate].measdim;
      case "chi2"_hash:
        return &m_index[istate].chi2;
      case "pathLength"_hash:
        return &m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &m_index[istate].typeFlags;
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
