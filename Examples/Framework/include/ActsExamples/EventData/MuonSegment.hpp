// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <cmath>
#include <vector>

namespace ActsExamples {
class MuonSegment {
 public:
  using MuonId = MuonSpacePoint::MuonId;
  /** @brief Empty default constructor */
  MuonSegment() = default;
  /** @brief Standard copy constructor */
  MuonSegment(const MuonSegment& other) = default;
  /** @brief Standard move constructor */
  MuonSegment(MuonSegment&& other) = default;
  /** @brief Standard copy assignment */
  MuonSegment& operator=(const MuonSegment& other) = default;

  /** @brief Returns the reconstructed segment position in global coordinates */
  const Acts::Vector3& globalPosition() const { return m_globPos; }
  /** @brief Returns the reconstructed segment direction in global coordinates */
  const Acts::Vector3& globalDirection() const { return m_globDir; }
  /** @brief Returns the reconstructed segment position expressed in the local space point frame*/
  const Acts::Vector3& localPosition() const { return m_localPos; }
  /** @brief Returns the reconstructed segment direction expressed in the local space point frame*/
  const Acts::Vector3& localDirection() const { return m_localDir; }
  /** @brief Returns the associated MS sector identifier */
  const MuonId& id() const { return m_id; }
  /** @brief returns the fitted segment time & uncertainty */
  double time() const { return m_time; }
  double timeUncert() const { return m_timeError; }
  /** @brief Returns the chiSquared & degreed of freedom */
  double chiSquared() const { return m_chiSquared; }
  double nDoF() const { return m_nDoF; }
  /** @brief Returns the number of precision hits building the segment */
  unsigned nPrecisionHits() const { return m_precisionHits; }
  /** @brief Returns the number of complementary trigger hits in the bending plane */
  unsigned nTrigEtaLayers() const { return m_triEtaLayers; }
  /** @brief Returns the number of complementary trigger hits in the non-bending plane */
  unsigned nTrigPhiLayers() const { return m_phiLayers; }
  /** @brief Set the global segment position & direction */
  void setGlobalCoords(const Acts::Vector3& pos, const Acts::Vector3& dir) {
    m_globPos = pos;
    m_globDir = dir;
  }
  /** @brief Set the local segment position & direction */
  void setLocalCoords(const Acts::Vector3& pos, const Acts::Vector3& dir) {
    m_localPos = pos;
    m_localDir = dir;
  }
  /** @brief Set the Identifier */
  void setId(const MuonId& id) { m_id = id; }
  /** @brief Set the time & unceratinty */
  void setTime(const double time, const double timeError) {
    m_time = time;
    m_timeError = timeError;
  }
  /** @brief Set the chi2 & the number of degrees of freedom */
  void setFitQuality(const double chi2, const unsigned nDoF) {
    m_chiSquared = chi2;
    m_nDoF = nDoF;
  }
  /** @brief Set the segment hit summary in terms of precision hits (Straw/ Mm /sTgc)
   *         & trigger eta (Rpc /Tgc) & trigger phi (Rpc / Tgc) hits  */
  void setHitSummary(unsigned nPrec, unsigned nTrigEta, unsigned nTrigPhi) {
    m_precisionHits = nPrec;
    m_triEtaLayers = nTrigEta;
    m_phiLayers = nTrigPhi;
  }

 private:
  MuonId m_id{};
  Acts::Vector3 m_globPos{Acts::Vector3::Zero()};
  Acts::Vector3 m_globDir{Acts::Vector3::Zero()};
  Acts::Vector3 m_localPos{Acts::Vector3::Zero()};
  Acts::Vector3 m_localDir{Acts::Vector3::Zero()};
  /** @brief Fitted segment time of arrival */
  double m_time{0.f};
  /** @brief Uncertainty on the time of arrival */
  double m_timeError{0.f};
  /** @brief segment chi2 & number of degrees of freedom */
  double m_chiSquared{0.f};
  unsigned m_nDoF{0u};
  /** @brief how many precision hits are on the segment (Straw tubes or Mm) */
  unsigned m_precisionHits{0u};
  /** @brief  Complementary hits in the non-bending direction (Rpc / Tgc / sTgc) */
  unsigned m_phiLayers{0u};
  /** @brief  Complementary hits in the bending direction (Rpc / Tgc) */
  unsigned m_triEtaLayers{0u};
};
using MuonSegmentContainer = std::vector<MuonSegment>;

}  // namespace ActsExamples
