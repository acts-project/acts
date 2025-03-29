
// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include <format>
namespace ActsExamples{
    using TechField = MuonSpacePoint::MuonId::TechField;
    using StationName = MuonSpacePoint::MuonId::StationName;
    using DetSide = MuonSpacePoint::MuonId::DetSide;

    std::string to_string(const StationName st) {
        switch (st) {
            case StationName::BIS: return "BIS";
            case StationName::BIL: return "BIL";
            case StationName::BMS: return "BMS";
            case StationName::BML: return "BML";
            case StationName::BOS: return "BOS";
            case StationName::BOL: return "BOL";
            case StationName::BEE: return "BEE";
            case StationName::EIL: return "EIL";
            case StationName::EIS: return "EIS";
            case StationName::EMS: return "EMS";
            case StationName::EML: return "EML";
            case StationName::EOS: return "EOS";
            case StationName::EOL: return "EOL";
            case StationName::EES: return "EES";
            case StationName::EEL: return "EEL";
            default:
                return "Unknown";
        };
    }
    std::string to_string(const TechField tech) {
        switch (tech) {
            case TechField::Mdt: return "Mdt";
            case TechField::Tgc: return "Tgc";
            case TechField::Rpc: return "Rpc";
            case TechField::sTgc: return "sTgc";
            case TechField::Mm: return "Mm";
            default:
                return "Unknown";
        };
    }
    std::string to_string(const DetSide side) {
        switch (side) {
            case DetSide::A:
                return "A-side";
            case DetSide::C:
                return "C-side";
            default:
                return "Unknown";
        };
    }

    std::ostream& operator<<(std::ostream& ostr, const MuonSpacePoint::MuonId& id){
        ostr<<std::format("{:} in {:2d} on {:}", to_string(id.msStation()), id.sector(), to_string(id.side()));
        return ostr;
    }
    void MuonSpacePoint::MuonId::setChamber(StationName stName, DetSide side, int sector, TechField tech) {
        m_stName = stName;
        m_side = side;
        m_sector = sector;
        m_tech = tech;
    }
    void MuonSpacePoint::MuonId::setLayAndCh(uint8_t layer, uint16_t ch) {
        m_layer = layer;
        m_channel = ch;
    }
    void MuonSpacePoint::MuonId::setCoordFlags(bool measEta, bool measPhi) {
        m_measEta = measEta;
        m_measPhi = measPhi;
    }

    std::ostream& operator<<(const MuonSpacePoint& sp, std::ostream& ostr) {
        ostr<<"Id: "<<sp.id()<<", pos: "<<Acts::toString(sp.localPosition())
            <<", dir: "<<Acts::toString(sp.sensorDirection());
        return ostr;
    }
    void MuonSpacePoint::defineCoordinates(Acts::Vector3&& pos, Acts::Vector3&& sensorDir) {
        m_pos = std::move(pos);
        m_dir = std::move(sensorDir);
    }
    void MuonSpacePoint::defineNormal(Acts::Vector3&& norm) {
        m_norm = std::move(norm);
    }
    void MuonSpacePoint::setRadius(const double r) {
        m_radius = r;
    }
    void MuonSpacePoint::setTime(const double t){
        m_time = t;
    }
    void MuonSpacePoint::setSpatialCov(const double xx, const double xy, 
                                       const double yx, const double yy) {
        m_cov(Acts::eX, Acts::eX) = xx;
        m_cov(Acts::eX, Acts::eY) = xy;
        m_cov(Acts::eY, Acts::eX) = yx;
        m_cov(Acts::eY, Acts::eY) = yy;
    }
    void MuonSpacePoint::setId(const MuonId & id) {
        m_id = id;
    }

    MuonSpacePointSorter::MuonSpacePointSorter(const SpVec_t& spacePoints) {
        m_strawHits.reserve(spacePoints.size());
        m_stripHits.reserve(spacePoints.size());
        for (const MuonSpacePoint* sp : spacePoints){
            LayerVec& pushMe{sp->id().technology() == MuonId::TechField::Mdt ? m_strawHits : m_stripHits};
            assert(sp->id().detLayer() >=1);
            const unsigned idx = sp->id().detLayer() - 1;
            assert(idx >=0);
            /// Ensure that there's enough space
            if (idx >= pushMe.size()) {
                pushMe.resize(idx +1);
            }
            pushMe[idx].push_back(sp);
        }
        /** The lowest strip gasGap number is the max +1 straw hit gasGap number. Erase all the hit vectors which are empty */
        m_stripHits.erase(m_stripHits.begin(), std::remove_if(m_stripHits.begin(), m_stripHits.end(),[](const SpVec_t& lay){ 
                         return lay.empty();}));
    }



}