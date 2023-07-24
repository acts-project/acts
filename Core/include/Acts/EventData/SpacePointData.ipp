// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

namespace Acts {

inline float SpacePointData::x(const std::size_t idx) const {
  return m_x[idx];
}

inline float SpacePointData::y(const std::size_t idx) const {
  return m_y[idx];
}

inline float SpacePointData::z(const std::size_t idx) const {
  return m_z[idx];
}

inline float SpacePointData::radius(const std::size_t idx) const {
  return m_radius[idx];
}

inline float SpacePointData::phi(const std::size_t idx) const {
  return m_phi[idx];
}

inline float SpacePointData::varianceZ(const std::size_t idx) const {
  return m_varianceZ[idx];
}

inline float SpacePointData::varianceR(const std::size_t idx) const {
  return m_varianceR[idx];
}

inline void SpacePointData::setX(const std::size_t idx, const float value) {
  m_x[idx] = value;
}

inline void SpacePointData::setY(const std::size_t idx, const float value) {
  m_y[idx] = value;
}

inline void SpacePointData::setZ(const std::size_t idx, const float value) {
  m_z[idx] = value;
}

inline void SpacePointData::setRadius(const std::size_t idx,
                                      const float value) {
  m_radius[idx] = value;
}

inline void SpacePointData::setPhi(const std::size_t idx, const float value) {
  m_phi[idx] = value;
}

inline void SpacePointData::setVarianceZ(const std::size_t idx,
                                         const float value) {
  m_varianceZ[idx] = value;
}

inline void SpacePointData::setVarianceR(const std::size_t idx,
                                         const float value) {
  m_varianceR[idx] = value;
}

inline float SpacePointData::quality(const std::size_t idx) const {
  return m_quality[idx];
}

inline float SpacePointData::deltaR(const std::size_t idx) const {
  return m_deltaR[idx];
}

inline void SpacePointData::setQuality(const std::size_t idx,
                                       const float value) const {
  if (value > m_quality[idx]) {
    m_quality[idx] = value;
  }
}

inline void SpacePointData::setDeltaR(const std::size_t idx,
                                      const float value) const {
  m_deltaR[idx] = value;
}

inline bool SpacePointData::hasDynamicVariable() const {
  return not m_topStripVector.empty();
}

inline const Acts::Vector3& SpacePointData::topStripVector(
    const std::size_t idx) const {
  return m_topStripVector[idx];
}

inline const Acts::Vector3& SpacePointData::bottomStripVector(
    const std::size_t idx) const {
  return m_bottomStripVector[idx];
}

inline const Acts::Vector3& SpacePointData::stripCenterDistance(
    const std::size_t idx) const {
  return m_stripCenterDistance[idx];
}

inline const Acts::Vector3& SpacePointData::topStripCenterPosition(
    const std::size_t idx) const {
  return m_topStripCenterPosition[idx];
}

inline void SpacePointData::setTopStripVector(const std::size_t idx,
                                              const Acts::Vector3& value) {
  m_topStripVector[idx] = value;
}

inline void SpacePointData::setBottomStripVector(const std::size_t idx,
                                                 const Acts::Vector3& value) {
  m_bottomStripVector[idx] = value;
}

inline void SpacePointData::setStripCenterDistance(const std::size_t idx,
                                                   const Acts::Vector3& value) {
  m_stripCenterDistance[idx] = value;
}

inline void SpacePointData::setTopStripCenterPosition(
    const std::size_t idx, const Acts::Vector3& value) {
  m_topStripCenterPosition[idx] = value;
}

inline void SpacePointData::resize(const std::size_t n, bool resizeDynamic) {
  clear();

  m_x.resize(n, 0.);
  m_y.resize(n, 0.);
  m_z.resize(n, 0.);
  m_radius.resize(n, 0.);
  m_phi.resize(n, 0.);
  m_varianceZ.resize(n, 0.);
  m_varianceR.resize(n, 0.);

  m_quality.resize(n, -std::numeric_limits<float>::infinity());
  m_deltaR.resize(n, 0.);

  if (resizeDynamic) {
    m_topStripVector.resize(n, {0, 0, 0});
    m_bottomStripVector.resize(n, {0, 0, 0});
    m_stripCenterDistance.resize(n, {0, 0, 0});
    m_topStripCenterPosition.resize(n, {0, 0, 0});
  }
}

inline void SpacePointData::clear() {
  //
  m_x.clear();
  m_y.clear();
  m_z.clear();
  m_radius.clear();
  m_phi.clear();
  m_varianceZ.clear();
  m_varianceR.clear();
  // mutable variables
  m_quality.clear();
  m_deltaR.clear();
  // dynamicvariables
  m_topStripVector.clear();
  m_bottomStripVector.clear();
  m_stripCenterDistance.clear();
  m_topStripCenterPosition.clear();
}

}  // namespace Acts
