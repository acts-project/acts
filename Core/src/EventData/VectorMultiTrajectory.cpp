// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/VectorMultiTrajectory.hpp"

namespace Acts {

VectorMultiTrajectory::DynamicColumnBase::~DynamicColumnBase() = default;

auto VectorMultiTrajectory::addTrackState_impl(TrackStatePropMask mask,
                                               IndexType iprevious)
    -> IndexType {
  using PropMask = TrackStatePropMask;

  m_index.emplace_back();
  IndexData& p = m_index.back();
  size_t index = m_index.size() - 1;

  if (iprevious != kInvalid) {
    p.iprevious = iprevious;
  }

  // always set, but can be null
  m_referenceSurfaces.emplace_back(nullptr);

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

  m_sourceLinks.push_back(nullptr);
  p.iuncalibrated = m_sourceLinks.size() - 1;

  if (ACTS_CHECK_BIT(mask, PropMask::Calibrated)) {
    m_meas.emplace_back();
    m_measCov.emplace_back();
    p.icalibrated = m_meas.size() - 1;

    m_sourceLinks.push_back(nullptr);
    p.icalibratedsourcelink = m_sourceLinks.size() - 1;

    m_projectors.emplace_back();
    p.iprojector = m_projectors.size() - 1;
  }

  // dynamic columns
  for (auto& [key, vec] : m_dynamic) {
    vec->add();
  }

  return index;
}

void VectorMultiTrajectory::shareFrom_impl(IndexType iself, IndexType iother,
                                           TrackStatePropMask shareSource,
                                           TrackStatePropMask shareTarget) {
  // auto other = getTrackState(iother);
  IndexData& self = m_index[iself];
  IndexData& other = m_index[iother];

  assert(ACTS_CHECK_BIT(getTrackState(iother).getMask(), shareSource) &&
         "Source has incompatible allocation");

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

void VectorMultiTrajectory::unset_impl(TrackStatePropMask target,
                                       IndexType istate) {
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

void VectorMultiTrajectory::clear_impl() {
  m_index.clear();
  m_params.clear();
  m_cov.clear();
  m_meas.clear();
  m_measCov.clear();
  m_jac.clear();
  m_sourceLinks.clear();
  m_projectors.clear();
  m_referenceSurfaces.clear();
  for (auto& [key, vec] : m_dynamic) {
    vec->clear();
  }
}

}  // namespace Acts
