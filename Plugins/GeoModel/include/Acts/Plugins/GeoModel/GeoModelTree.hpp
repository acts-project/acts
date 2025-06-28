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

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelRead/ReadGeoModel.h>

namespace Acts {

struct GeoModelTree {
  using FPVLink = GeoIntrusivePtr<GeoFullPhysVol>;
  using FPVConstLink = GeoIntrusivePtr<const GeoFullPhysVol>;
  /// @brief Helper struct to manage the of published full physical volumes per
  ///        sub detector. The volumes are stored in std::map
  ///
  struct VolumePublisher {
    /// @brief Abrivation of the published full physical volumes
    using VolumeMap_t = std::map<std::string, FPVConstLink>;
    /// @brief Abrivation of the published map per system
    using PublisherMap_t = std::unordered_map<std::string, VolumeMap_t>;
    /// @brief Default constructor
    VolumePublisher() = default;
    /// @brief Constructor taking an instance to the GeoModelReader.
    ///        If the publisher does not know yet about the system,
    ///        the published volumes are queried from the database
    VolumePublisher(std::shared_ptr<GeoModelIO::ReadGeoModel> geoReader);
    /// @brief Copy constructor
    VolumePublisher(const VolumePublisher& other);
    /// @brief Move constructor
    VolumePublisher(VolumePublisher&& other);
    /// @brief Copy assignment operator
    VolumePublisher& operator=(const VolumePublisher& other);
    /// @brief Move assignment operator
    VolumePublisher& operator=(VolumePublisher&& other);
    /// @brief Returns the
    const VolumeMap_t& getPublishedVol(const std::string& systemName) const;
    /// @brief  Declare a list of published
    /// @param systemName
    /// @param publishedVols
    void publishVolumes(const std::string& systemName,
                        VolumeMap_t&& publishedVols);

    const std::shared_ptr<GeoModelIO::ReadGeoModel>& reader() const {
      return m_reader;
    }

   private:
    std::shared_ptr<GeoModelIO::ReadGeoModel> m_reader{};
    mutable PublisherMap_t m_publishedVols{};
    mutable std::mutex m_mutex{};
  };
  /// @brief Empty default constructor
  GeoModelTree() = default;
  /// @brief Constructor taking a shared pointer to an opened
  ///        database session.
  /// @param db: Valid pointer to a GeoModel database manager
  GeoModelTree(std::shared_ptr<GMDBManager> db);

  std::shared_ptr<GMDBManager> dbMgr{};
  VolumePublisher publisher{};
  PVConstLink worldVolume{};
  std::string worldVolumeName{"World"};
};

}  // namespace Acts
