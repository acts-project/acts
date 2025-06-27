// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <mutex>

#include <GeoModelRead/ReadGeoModel.h>

class GeoVPhysVol;

namespace Acts {


struct GeoModelTree {
  /// @brief Helper struct to manage the list of full physical volumes 
  ///        per detector system.
  struct VolumePublisher{
    /// @brief Abrivation of the published full physical volumes
    using VolumeMap_t = std::map<std::string, const GeoVFullPhysVol*>;
    /// @brief Abrivation of the published map per system
    using PublisherMap_t = std::unordered_map<std::string, VolumeMap_t>;
    /// @brief Default constructor
    VolumePublisher() = default;
    /// @brief Constructor taking an instance to the GeoModelReader.
    ///        If the publisher does not know yet about the system,
    ///        the published volumes are queried from the database
    VolumePublisher(std::shared_ptr<GeoModelIO::ReadGeoModel> geoReader);

    VolumePublisher(const VolumePublisher& other);
    VolumePublisher(VolumePublisher&& other);
    VolumePublisher& operator=(const VolumePublisher& other);
    VolumePublisher& operator=(VolumePublisher&& other);
    /// @brief Returns the 
    const VolumeMap_t& getPublishedVol(const std::string& systemName) const;

    void publishVolumes(const std::string& systemName,
                         VolumeMap_t&&  publishedVols);
    
    const std::shared_ptr<GeoModelIO::ReadGeoModel>& reader () const{
        return m_reader;
    }
    private:
        std::shared_ptr<GeoModelIO::ReadGeoModel> m_reader{};
        mutable PublisherMap_t m_publishedVols{};
        mutable std::mutex m_mutex{};
  };
  /// @brief Empty default constructor
  GeoModelTree() = default;
  /// @brief Copy constructor
  GeoModelTree(const GeoModelTree& other) = default;
  /// @brief Move constructor
  GeoModelTree(GeoModelTree& other) = default;
  /// @brief Copy assignment 
  GeoModelTree& operator=(const GeoModelTree& other) = default;
  /// @brief Move assignment
  GeoModelTree& operator=(GeoModelTree&& other) = default;
  /// @brief Constructor taking a 
  GeoModelTree(std::shared_ptr<GMDBManager> db);


  std::shared_ptr<GMDBManager> dbMgr{};
  VolumePublisher publisher{};
  PVConstLink worldVolume{};
};

}  // namespace Acts
