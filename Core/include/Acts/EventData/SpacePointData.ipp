// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

inline const float& SpacePointData::x(const std::size_t& idx) const {
  return m_x[idx];
}
  
inline const float& SpacePointData::y(const std::size_t& idx) const {
  return m_y[idx];
}
  
inline const float& SpacePointData::z(const std::size_t& idx) const {
  return m_z[idx];
}

inline const float& SpacePointData::radius(const std::size_t& idx) const {
  return m_radius[idx];
}
  
inline const float& SpacePointData::phi(const std::size_t& idx) const {
  return m_phi[idx];
}

inline const float& SpacePointData::varianceZ(const std::size_t& idx) const {
  return m_varianceZ[idx];
}

inline const float& SpacePointData::varianceR(const std::size_t& idx) const {
  return m_varianceR[idx];
}

  
inline void SpacePointData::setX(const std::size_t& idx, const float& value) {
  m_x[idx] = value;
}

inline void SpacePointData::setY(const std::size_t& idx, const float& value) {
  m_y[idx] = value;
}

inline void SpacePointData::setZ(const std::size_t& idx, const float& value) {
  m_z[idx] = value;
}

inline void SpacePointData::setRadius(const std::size_t& idx, const float& value) {
  m_radius[idx] = value;
}

inline void SpacePointData::setPhi(const std::size_t& idx, const float& value) {
  m_phi[idx] = value;
}

inline void SpacePointData::setVarianceZ(const std::size_t& idx, const float& value) {
  m_varianceZ[idx] = value;
}
  
inline void SpacePointData::setVarianceR(const std::size_t& idx, const float& value) {
  m_varianceR[idx] = value;
}
  
inline const float& SpacePointData::quality(const std::size_t& idx) const {
  return m_quality[idx];
}

inline const float& SpacePointData::deltaR(const std::size_t& idx) const {
  return m_deltaR[idx];
}

inline void SpacePointData::setQuality(const std::size_t& idx, const float& value) {
  if (value > m_quality[idx]) {
    m_quality[idx] = value;
  }
}

inline void SpacePointData::setDeltaR(const std::size_t& idx, const float& value) {
  m_deltaR[idx] = value;
}

inline bool SpacePointData::hasDynamicVariable() const {
  return not m_topHalfStripLength.empty();
}

inline std::any SpacePointData::component(Acts::HashedString key,
					  const std::size_t& n) const {
  using namespace Acts::HashedStringLiteral;
  switch (key) {
  case "TopHalfStripLength"_hash:
    return m_topHalfStripLength[n];
  case "BottomHalfStripLength"_hash:
    return m_bottomHalfStripLength[n];
  case "TopStripDirection"_hash:
    return m_topStripDirection[n];
  case "BottomStripDirection"_hash:
    return m_bottomStripDirection[n];
  case "StripCenterDistance"_hash:
    return m_stripCenterDistance[n];
  case "TopStripCenterPosition"_hash:
    return m_topStripCenterPosition[n];
  default:
    throw std::runtime_error("no such component " + std::to_string(key));
  }
}
  
inline const float& SpacePointData::getTopHalfStripLength(
    const std::size_t& idx) const {
  return m_topHalfStripLength[idx];
}

inline const float& SpacePointData::getBottomHalfStripLength(
    const std::size_t& idx) const {
  return m_bottomHalfStripLength[idx];
}

inline const Acts::Vector3& SpacePointData::getTopStripDirection(
    const std::size_t& idx) const {
  return m_topStripDirection[idx];
}

inline const Acts::Vector3& SpacePointData::getBottomStripDirection(
    const std::size_t& idx) const {
  return m_bottomStripDirection[idx];
}

inline const Acts::Vector3& SpacePointData::getStripCenterDistance(
    const std::size_t& idx) const {
  return m_stripCenterDistance[idx];
}

inline const Acts::Vector3& SpacePointData::getTopStripCenterPosition(
    const std::size_t& idx) const {
  return m_topStripCenterPosition[idx];
}

inline void SpacePointData::setTopHalfStripLength(const std::size_t& idx,
                                                  const float& value) {
  m_topHalfStripLength[idx] = value;
}

inline void SpacePointData::setBottomHalfStripLength(const std::size_t& idx,
                                                     const float& value) {
  m_bottomHalfStripLength[idx] = value;
}

inline void SpacePointData::setTopStripDirection(const std::size_t& idx,
                                                 const Acts::Vector3& value) {
  m_topStripDirection[idx] = value;
}

inline void SpacePointData::setBottomStripDirection(
    const std::size_t& idx, const Acts::Vector3& value) {
  m_bottomStripDirection[idx] = value;
}

inline void SpacePointData::setStripCenterDistance(const std::size_t& idx,
                                                   const Acts::Vector3& value) {
  m_stripCenterDistance[idx] = value;
}

inline void SpacePointData::setTopStripCenterPosition(
    const std::size_t& idx, const Acts::Vector3& value) {
  m_topStripCenterPosition[idx] = value;
}

inline void SpacePointData::resize(const std::size_t& n, bool resizeDynamic) {
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
    m_topHalfStripLength.resize(n, 0.);
    m_bottomHalfStripLength.resize(n, 0.);
    m_topStripDirection.resize(n, {0, 0, 0});
    m_bottomStripDirection.resize(n, {0, 0, 0});
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
  m_topHalfStripLength.clear();
  m_bottomHalfStripLength.clear();
  m_topStripDirection.clear();
  m_bottomStripDirection.clear();
  m_stripCenterDistance.clear();
  m_topStripCenterPosition.clear();
}

}  // namespace Acts
