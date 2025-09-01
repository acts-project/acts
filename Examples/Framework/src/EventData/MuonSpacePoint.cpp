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

#include <format>

namespace ActsExamples {

std::ostream& operator<<(std::ostream& ostr,
                         const ActsExamples::MuonSpacePoint::MuonId& id) {
  ostr << id.toString();
  return ostr;
}
std::ostream& operator<<(std::ostream& ostr,
                         const ActsExamples::MuonSpacePoint& sp) {
  ostr << std::format(
      "Muon-ID: {}, position: {}, orientation (d/n/o): {}/{}/{}, covariance: "
      "({:.2f}, {:.2f}, {:.2f})",
      sp.id().toString(), Acts::toString(sp.localPosition()),
      Acts::toString(sp.sensorDirection()), Acts::toString(sp.toNextSensor()),
      Acts::toString(sp.planeNormal()), sp.covariance()[0], sp.covariance()[1],
      sp.covariance()[2]);
  return ostr;
}

using TechField = MuonSpacePoint::MuonId::TechField;
using StationName = MuonSpacePoint::MuonId::StationName;
using DetSide = MuonSpacePoint::MuonId::DetSide;

std::string MuonSpacePoint::MuonId::toString(const StationName st) {
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
std::string MuonSpacePoint::MuonId::toString(const TechField tech) {
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
std::string MuonSpacePoint::MuonId::toString(const DetSide side) {
  switch (side) {
    case DetSide::A:
      return "A-side";
    case DetSide::C:
      return "C-side";
    default:
      return "Unknown";
  }
}
std::string MuonSpacePoint::MuonId::toString() const {
  return std::format("{:} in {:2d} on {:}", toString(msStation()), sector(),
                     toString(side()));
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
                                       Acts::Vector3&& sensorDir,
                                       Acts::Vector3&& toNextSensor) {
  m_pos = std::move(pos);
  m_dir = std::move(sensorDir);
  m_toNext = std::move(toNextSensor);
  m_norm = m_dir.cross(m_toNext).normalized();
}
void MuonSpacePoint::setRadius(const double r) {
  m_radius = r;
}
void MuonSpacePoint::setTime(const double t) {
  m_time = t;
}
void MuonSpacePoint::setCovariance(const double covX, const double covY,
                                   const double covT) {
  m_cov[0] = covX;
  m_cov[1] = covY;
  m_cov[2] = covT;
}
void MuonSpacePoint::setId(const MuonId& id) {
  m_id = id;
}
void MuonSpacePoint::setGeometryId(const Acts::GeometryIdentifier& geoId) {
  m_geoId = geoId;
}
}  // namespace ActsExamples
