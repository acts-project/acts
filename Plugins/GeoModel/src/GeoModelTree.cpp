// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"

namespace Acts {
GeoModelTree::VolumePublisher::VolumePublisher(
    const std::shared_ptr<GeoModelIO::ReadGeoModel>& geoReader) noexcept
    : m_reader{geoReader} {}

const GeoModelTree::VolumePublisher::VolumeMap_t&
GeoModelTree::VolumePublisher::getPublishedVol(const std::string& systemName) {
  PublisherMap_t::const_iterator find_itr = m_publishedVols.find(systemName);
  if (find_itr != m_publishedVols.end()) {
    return find_itr->second;
  }
  if (m_reader) {
    auto qFPV =
        m_reader->getPublishedNodes<std::string, GeoFullPhysVol*>(systemName);
    VolumeMap_t newPublished{};
    newPublished.insert(qFPV.begin(), qFPV.end());
    return m_publishedVols
        .insert(std::make_pair(systemName, std::move(newPublished)))
        .first->second;
  }
  throw std::domain_error(
      "GeoModelTree::VolumePublisher() - No volumes have been published for " +
      systemName);
}
void GeoModelTree::VolumePublisher::publishVolumes(
    const std::string& systemName, VolumeMap_t&& publishedVols) {
  if (systemName.empty()) {
    throw std::logic_error(
        "GeoModelTree::VolumePublisher() - System name cannot be empty");
  }
  if (!m_publishedVols
           .insert(std::make_pair(systemName, std::move(publishedVols)))
           .second) {
    throw std::domain_error("GeoModelTree::VolumePublisher - System " +
                            systemName + " has already published volumes.");
  }
}
GeoModelTree::GeoModelTree(const std::shared_ptr<GMDBManager>& db)
    : dbMgr{db},
      publisher{std::make_shared<VolumePublisher>(
          std::make_shared<GeoModelIO::ReadGeoModel>(db.get()))},
      worldVolume{publisher->reader()->buildGeoModel()} {}

}  // namespace Acts
