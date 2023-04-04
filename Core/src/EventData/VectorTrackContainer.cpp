// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/VectorTrackContainer.hpp"

namespace Acts {

namespace detail_vtc {

VectorTrackContainerBase::VectorTrackContainerBase(
    const VectorTrackContainerBase& other)
    : m_tipIndex{other.m_tipIndex},
      m_params{other.m_params},
      m_cov{other.m_cov},
      m_referenceSurfaces{other.m_referenceSurfaces},
      m_nMeasurements{other.m_nMeasurements},
      m_nHoles{other.m_nHoles} {
  for (const auto& [key, value] : other.m_dynamic) {
    m_dynamic.insert({key, value->clone()});
  }
  assert(checkConsistency());
}
}  // namespace detail_vtc

VectorTrackContainer::IndexType VectorTrackContainer::addTrack_impl() {
  assert(checkConsistency());

  m_tipIndex.emplace_back(kInvalid);

  m_params.emplace_back();
  m_cov.emplace_back();
  m_referenceSurfaces.emplace_back();

  m_nMeasurements.emplace_back();
  m_nHoles.emplace_back();

  // dynamic columns
  for (auto& [key, vec] : m_dynamic) {
    vec->add();
  }

  assert(checkConsistency());

  return m_tipIndex.size() - 1;
}

void VectorTrackContainer::removeTrack_impl(IndexType itrack) {
  auto erase = [&](auto& vec) {
    assert(itrack < vec.size() && "Index is out of range");
    auto it = vec.begin();
    std::advance(it, itrack);
    vec.erase(it);
  };

  erase(m_tipIndex);

  erase(m_params);
  erase(m_cov);
  erase(m_referenceSurfaces);

  erase(m_nMeasurements);
  erase(m_nHoles);

  for (auto& [key, vec] : m_dynamic) {
    vec->erase(itrack);
  }
}

void VectorTrackContainer::copyDynamicFrom_impl(
    IndexType dstIdx, const VectorTrackContainerBase& src, IndexType srcIdx) {
  for (const auto& [key, value] : src.m_dynamic) {
    auto it = m_dynamic.find(key);
    if (it == m_dynamic.end()) {
      throw std::invalid_argument{
          "Destination container does not have matching dynamic column"};
    }

    it->second->copyFrom(dstIdx, *value, srcIdx);
  }
}

void VectorTrackContainer::ensureDynamicColumns_impl(
    const detail_vtc::VectorTrackContainerBase& other) {
  for (auto& [key, value] : other.m_dynamic) {
    if (m_dynamic.find(key) == m_dynamic.end()) {
      m_dynamic[key] = value->clone(true);
    }
  }
}

void VectorTrackContainer::reserve(IndexType size) {
  m_tipIndex.reserve(size);

  m_params.reserve(size);
  m_cov.reserve(size);
  m_referenceSurfaces.reserve(size);

  m_nMeasurements.reserve(size);
  m_nHoles.reserve(size);

  for (auto& [key, vec] : m_dynamic) {
    vec->reserve(size);
  }
}

}  // namespace Acts
