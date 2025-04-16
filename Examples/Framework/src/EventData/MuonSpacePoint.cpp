// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include "Acts/Utilities/StringHelpers.hpp"

#include "TString.h"

std::ostream& operator<<(std::ostream& ostr,
                         const ActsExamples::MuonSpacePoint::MuonId& id) {
  ostr << Form("%s in %2d on %s",
               ActsExamples::to_string(id.msStation()).c_str(), id.sector(),
               ActsExamples::to_string(id.side()).c_str());
  return ostr;
}
std::ostream& operator<<(std::ostream& ostr,
                         const ActsExamples::MuonSpacePoint& sp) {
  ostr << "Id: " << sp.id() << ", pos: " << Acts::toString(sp.localPosition())
       << ", dir: " << Acts::toString(sp.sensorDirection())
       << ", covariance: " << Acts::toString(sp.covariance());
  return ostr;
}

namespace ActsExamples {
using TechField = MuonSpacePoint::MuonId::TechField;
using StationName = MuonSpacePoint::MuonId::StationName;
using DetSide = MuonSpacePoint::MuonId::DetSide;

std::string to_string(const StationName st) {
  switch (st) {
    case StationName::BIS:
      return "BIS";
    case StationName::BIL:
      return "BIL";
    case StationName::BMS:
      return "BMS";
    case StationName::BML:
      return "BML";
    case StationName::BOS:
      return "BOS";
    case StationName::BOL:
      return "BOL";
    case StationName::BEE:
      return "BEE";
    case StationName::EIL:
      return "EIL";
    case StationName::EIS:
      return "EIS";
    case StationName::EMS:
      return "EMS";
    case StationName::EML:
      return "EML";
    case StationName::EOS:
      return "EOS";
    case StationName::EOL:
      return "EOL";
    case StationName::EES:
      return "EES";
    case StationName::EEL:
      return "EEL";
    default:
      return "Unknown";
  }
}
std::string to_string(const TechField tech) {
  switch (tech) {
    case TechField::Mdt:
      return "Mdt";
    case TechField::Tgc:
      return "Tgc";
    case TechField::Rpc:
      return "Rpc";
    case TechField::sTgc:
      return "sTgc";
    case TechField::Mm:
      return "Mm";
    default:
      return "Unknown";
  }
}
std::string to_string(const DetSide side) {
  switch (side) {
    case DetSide::A:
      return "A-side";
    case DetSide::C:
      return "C-side";
    default:
      return "Unknown";
  }
}
void MuonSpacePoint::MuonId::setChamber(StationName stName, DetSide side,
                                        int sector, TechField tech) {
  m_stName = stName;
  m_side = side;
  m_sector = sector;
  m_tech = tech;
}
void MuonSpacePoint::MuonId::setLayAndCh(std::uint8_t layer, std::uint16_t ch) {
  m_layer = layer;
  m_channel = ch;
}
void MuonSpacePoint::MuonId::setCoordFlags(bool measEta, bool measPhi) {
  m_measEta = measEta;
  m_measPhi = measPhi;
}
void MuonSpacePoint::defineCoordinates(Acts::Vector3&& pos,
                                       Acts::Vector3&& sensorDir) {
  m_pos = std::move(pos);
  m_dir = std::move(sensorDir);
}
void MuonSpacePoint::defineNormal(Acts::Vector3&& norm) {
  m_norm = std::move(norm);
}
void MuonSpacePoint::setRadius(const double r) {
  m_radius = r;
}
void MuonSpacePoint::setTime(const double t) {
  m_time = t;
}
void MuonSpacePoint::setSpatialCov(const double xx, const double xy,
                                   const double yx, const double yy) {
  m_cov(Acts::eX, Acts::eX) = xx;
  m_cov(Acts::eX, Acts::eY) = xy;
  m_cov(Acts::eY, Acts::eX) = yx;
  m_cov(Acts::eY, Acts::eY) = yy;
}
void MuonSpacePoint::setId(const MuonId& id) {
  m_id = id;
}
}  // namespace ActsExamples
