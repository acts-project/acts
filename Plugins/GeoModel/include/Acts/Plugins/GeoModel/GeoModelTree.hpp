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
/// @brief Holder struct to manage a GeoModel world. It holds the pointer to the
///        root volume of the world and provides additional links to dedicated
///        volumes inside the tree representing the sensors and chambers in a
///        detector.
struct GeoModelTree {
  using FpvLink = GeoIntrusivePtr<GeoFullPhysVol>;
  using FpvConstLink = GeoIntrusivePtr<const GeoFullPhysVol>;
  /// @brief Helper struct to manage the list of published full physical volumes.
  ///        These volumes have unique name tags and seed the construction of
  ///        tracking volumes. The list of published full physical volumes are
  ///        the name tags associated with a particular node in the GeoModel
  ///        tree. The lists are usually separated by subdetector and are
  ///        shipped together with the SQLite database. The latter is queried,
  ///        if the geometry is constructed from an external file. If the
  ///        geometry is built dynamically in memory, the VolumePublisher
  ///        provides the mechanism to store and return the published volumes.
  struct VolumePublisher {
    /// @brief Abrivation of the published full physical volumes
    using VolumeMap_t = std::map<std::string, FpvConstLink>;
    /// @brief Abrivation of the published map per system
    using PublisherMap_t = std::unordered_map<std::string, VolumeMap_t>;
    /// @brief Default constructor
    VolumePublisher() = default;
    /// @brief Constructor taking an instance to the GeoModelReader.
    ///        If the publisher does not know yet about the system,
    ///        the published volumes are queried from the database
    explicit VolumePublisher(
        const std::shared_ptr<GeoModelIO::ReadGeoModel>& geoReader) noexcept;

    /// @brief Returns the list of published full physical volume for a given subsystem
    /// @param systemName: System of interest (e.g. Tracker, MS)
    const VolumeMap_t& getPublishedVol(const std::string& systemName);
    /// @brief  Publish a list of full physical volumes for a given subsystem
    /// @param systemName: Subsystem name under which the list will be registered
    /// @param publishedVols: List of volumes to publish.
    void publishVolumes(const std::string& systemName,
                        VolumeMap_t&& publishedVols);
    /// @brief Returns the pointer to the ReadGeoModel instance (might be invalid)
    const std::shared_ptr<GeoModelIO::ReadGeoModel>& reader() const {
      return m_reader;
    }

   private:
    /// @brief Pointer to the ReadGeoModel instance
    std::shared_ptr<GeoModelIO::ReadGeoModel> m_reader{};
    /// @brief Map of all published full phyical volumes per subsytstem
    PublisherMap_t m_publishedVols{};
  };

  /// @brief Empty default constructor
  GeoModelTree() = default;
  /// @brief Constructor taking a shared pointer to an opened
  ///        SQLite DB session.
  /// @param db: Valid pointer to a GeoModel database manager
  explicit GeoModelTree(const std::shared_ptr<GMDBManager>& db);

  /// @brief Database manager used to construct the Root node and
  ///        also to query information from additional tables
  std::shared_ptr<GMDBManager> dbMgr{};
  /// @brief Manager class of the published volumes per subsystem
  std::shared_ptr<VolumePublisher> publisher{
      std::make_shared<VolumePublisher>()};
  /// @brief Root node of the GeoModel world
  PVConstLink worldVolume{};
  /// @brief Name of the Root node
  std::string worldVolumeName{"World"};
};

}  // namespace Acts
