// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/EventData/StationSpacePoint.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <iostream>
#include <memory>

namespace ActsExamples{
    /** @brief Example implementation of a StationSpacePoint concept inspired by the ATLAS Muon::SpacePoint EDM. 
     *         The space points are expressed in a local frame such that the x-axis is parallel to the ATLAS 
     *         Monitored Drift Tubes (Mdt), the y-axis points within the tube layer & the z-axis outside of the plane */
    class MuonSpacePoint {
        public:   
            /** @brief Identifier helper class to distinguish different measurements */         
            class MuonId {
                public:
                    /** @brief Technology encoding of the measurement*/
                    enum class TechField: int8_t {
                        UnDef = -1, Mdt = 0, Rpc = 2, Tgc = 3, sTgc = 4, Mm = 5
                    };
                    /** @brief ATLAS-MS station encoding */
                    enum class StationName: int8_t {
                        UnDef = -1, BIS, BIL, BMS, BML, BOS, BOL, BEE,
                        EIS, EIL, EMS, EML, EOS, EOL, EES, EEL
                    };
                    /** @brief Detector side encoding */
                    enum class DetSide: int8_t {
                       UnDef = 0, A = 1, C = -1
                    };
                    /** @brief Empty default Identifier constructor */
                    explicit MuonId() = default;
                    /** @brief Default copy constructor */
                    MuonId(const MuonId& other) = default;
                    /** @brief Default move constructor */
                    MuonId(MuonId && other) = default;
                    /** @brief Default copy assignment */
                    MuonId& operator=(const MuonId& other) = default;
                    /** @brief Default move assignment */
                    MuonId& operator=(MuonId&& other) = default;

                    /** @brief Returns the technology of the measurement */
                    TechField technology() const {
                        return m_tech;
                    }
                    /** @brief Returns the MS station in which the measurement was recorded */
                    StationName msStation() const {
                        return m_stName;
                    }
                    /** @brief Returns the sector in which the measurement was recorded */
                    uint8_t sector() const {
                        return m_sector;
                    }
                    /** @brief Returns the detector side */
                    DetSide side() const {
                        return m_side;
                    }
                    /** @brief Returns the layer */
                    uint8_t detLayer() const {
                        return m_layer;
                    }
                    /** @brief Returns the channel number */
                    uint16_t channel() const { 
                        return m_channel;
                    }
                    /** @brief Returns whether the id corresponds to a precision coordinate (eta) measurement */
                    bool measuresEta() const { return m_measEta; }
                    /**  @brief Returns whether the id corresponds to a non-precision coordinate (phi) measurement */
                    bool measuresPhi() const { return m_measPhi; }
                    /** @brief Returns whether two Identifiers belong to the same station, which
                     *         is characterized that both share the same msStation, sector & side field. */
                    bool sameStation(const MuonId& other) const {
                        return msStation() == other.msStation() &&
                               sector() == other.sector() &&
                               side() == other.side();
                    }
                    /** @brief Set the fields neede to Identify the detector in the system
                     *  @param stName: StationName where the muon station is located
                     *  @param side: Positive or negative side
                     *  @param sector: Phi sector in which the chamber is installed
                     *  @param tech: Technology of the chamber within the chamber */
                    void setChamber(StationName stName, DetSide side, int sector, TechField tech);
                    /** @brief Set the measurement layer & channel */
                    void setLayAndCh(uint8_t layer, uint16_t ch);
                    /** @brief Define the measurement type of the space point
                     *  @param measEta: Flag stating whether the space point measures the precision (eta) coordinate
                     *  @param measPhi: Flag stating whether the space point measures the non-precsion (phi) coordinate */
                    void setCoordFlags(bool measEta, bool measPhi);
                private:
                    TechField m_tech{TechField::UnDef};
                    StationName m_stName{StationName::UnDef};
                    DetSide m_side{DetSide::UnDef};
                    uint8_t m_sector{0};
                    uint8_t m_layer{0};
                    uint16_t m_channel{0}; 
                    bool m_measEta{false};
                    bool m_measPhi{false};
            };
            /** @brief Empty default constructor */
            MuonSpacePoint() = default;
            /** @brief Copy constructor */
            MuonSpacePoint(const MuonSpacePoint& other) = default;
            /** @brief Move constructor */
            MuonSpacePoint(MuonSpacePoint&& other) = default;
            /** @brief Returns the Identifier of the space point */
            const MuonId& id() const {
                return m_id;
            }
            /** @brief Returns the local meaurement position */
            const Acts::Vector3& localPosition() const{
                return m_pos;
            }
            /** @brief Returns the local sensor direction */
            const Acts::Vector3& sensorDirection() const{
                return m_dir;
            }
            /** @brief Returns the normal vector to the plane */
            const Acts::Vector3& stripPlaneNormal() const {
                return m_pos;
            }
            /** @brief Returns the space point covariance */
            const Acts::ActsSquareMatrix<3>& covariance() const {
                return m_cov;
            }
            /** @brief Returns the drift radius */
            double driftRadius() const{
                return m_radius;
            }
            /** @brief Returns the measurement time */
            double time() const {
                return m_time;
            }
            /** @brief Define the space point's identifier */
            void setId(const MuonId & id);
            /** @brief Define the space point coordinates.
             *  @param pos: Space point position
             *  @param sensorDir: Direction of the sensor */
            void defineCoordinates(Acts::Vector3&& pos, Acts::Vector3&& sensorDir);
            /** @brief Define the space point normal*/
            void defineNormal(Acts::Vector3&& norm);
            /** @brief Define the space point radius */
            void setRadius(const double r);
            /** @brief Define the time of the space point measurement */
            void setTime(const double t);
            /** @brief Define the spatial components of the covariance */
            void setSpatialCov(const double xx, const double xy, 
                               const double yx, const double yy);
        private:
            MuonId m_id{};
            Acts::Vector3 m_pos{Acts::Vector3::Zero()};
            Acts::Vector3 m_dir{Acts::Vector3::Zero()};
            Acts::Vector3 m_norm{Acts::Vector3::Zero()};
            Acts::ActsSquareMatrix<3> m_cov{Acts::ActsSquareMatrix<3>::Identity()};
            double m_radius{0.};
            double m_time{0.};
    };
    using MuonSpacePointContainer = std::vector<std::vector<MuonSpacePoint>>;

    /** @brief Helper class to sort the MuonspacePoints whether they're Mdt hits, i.e. straw hits or whether they're ordinary strip hits.\
      *        Hits in each category are additionally sorted by straw or tube layer. Each layer is represented by a dedicated hit vector. 
      *        The sorter then holds the sorted vectors in two jagged vectors and returns them via its interface.      */
    class MuonSpacePointSorter{
        public:
            /** @brief abbreviation of the space point vector */
            using SpVec_t = std::vector<const MuonSpacePoint*>;
            /** @brief abbreviation of the layer vector. Each element inside the vector represents
             *         all hits in a given detector layer */
            using LayerVec = std::vector<SpVec_t>;

            using MuonId = MuonSpacePoint::MuonId;
            /** @brief Constructor taking a list of hits in the same muon station. */
            MuonSpacePointSorter(const SpVec_t& spacePoints);
            /** @brief Returns all straw hits sorted by straw layer */
            const LayerVec& strawHits() const{
                return m_strawHits;
            }
            /** @brief Returns all strip hits sorted by strip layer */
            const LayerVec& stripHits() const { 
                return m_stripHits;
            }
        private:
            LayerVec m_strawHits{};
            LayerVec m_stripHits{};
    };

    


    std::ostream& operator<<(std::ostream&ostr, const MuonSpacePoint::MuonId& id);
    std::ostream& operator<<(std::ostream&ostr, const MuonSpacePoint& sp);

    std::string to_string(const MuonSpacePoint::MuonId::StationName st);
    std::string to_string(const MuonSpacePoint::MuonId::TechField tech);
    std::string to_string(const MuonSpacePoint::MuonId::DetSide side);

}