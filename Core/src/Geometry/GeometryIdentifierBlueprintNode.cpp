// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifierBlueprintNode.hpp"

#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>

#include <boost/algorithm/string/join.hpp>

namespace Acts::Experimental {

namespace {
class Configuration {
 public:
  virtual ~Configuration() = default;

  virtual void apply(const std::string& prefix, TrackingVolume& volume,
                     const Logger& logger) = 0;
  virtual const std::string& name() const = 0;
};

struct FixedLayerConfiguration : public Configuration {
  explicit FixedLayerConfiguration(GeometryIdentifier::Value layer)
      : m_layer(layer) {
    m_name = "GeoIdFixLayer(lay=" + std::to_string(m_layer) + ")";
  }

  void apply(const std::string& prefix, TrackingVolume& volume,
             const Logger& logger) override {
    ACTS_DEBUG(prefix << "~> Setting layer to " << m_layer
                      << " for volume with ID " << volume.geometryId());
    volume.assignGeometryId(volume.geometryId().withLayer(m_layer));
  }

  const std::string& name() const override { return m_name; }

 private:
  GeometryIdentifier::Value m_layer;
  std::string m_name;
};

struct IncrementLayerConfiguration : public Configuration {
  explicit IncrementLayerConfiguration(GeometryIdentifier::Value start)
      : m_value(start) {
    m_name = "GeoIdIncLay(start=" + std::to_string(m_value) + ")";
  }

  void apply(const std::string& prefix, TrackingVolume& volume,
             const Logger& logger) override {
    ACTS_DEBUG(prefix << "Incrementing layer component for volume with ID "
                      << volume.geometryId());
    if (volume.geometryId().layer() != 0) {
      ACTS_ERROR("Volume " << volume.volumeName() << " already has layer ID "
                           << volume.geometryId().layer() << ". Please check "
                           << "your geometry configuration.");
      throw std::logic_error("Volume already has a layer ID");
    }
    GeometryIdentifier id = volume.geometryId().withLayer(m_value);
    ACTS_DEBUG(prefix << "~> Setting layer to " << m_value
                      << " for volume with ID " << id);
    volume.assignGeometryId(id);
    m_value++;
  }

  const std::string& name() const override { return m_name; }

 private:
  GeometryIdentifier::Value m_value;
  std::string m_name;
};

struct FixedVolumeConfiguration : public Configuration {
  explicit FixedVolumeConfiguration(GeometryIdentifier::Value volumeId)
      : m_volumeId(volumeId) {
    m_name = "GeoIdFixVol(vol=" + std::to_string(m_volumeId) + ")";
  }

  void apply(const std::string& prefix, TrackingVolume& volume,
             const Logger& logger) override {
    ACTS_DEBUG(prefix << "~> Setting volume ID to " << m_volumeId
                      << " for volume " << volume.volumeName()
                      << " and all descendents");
    volume.apply([&](TrackingVolume& v) {
      if (v.geometryId().volume() != 0) {
        ACTS_ERROR("Volume " << v.volumeName() << " already has volume ID "
                             << v.geometryId().volume()
                             << ". Please check your geometry configuration.");
        throw std::logic_error("Volume already has a volume ID");
      }
      ACTS_DEBUG(prefix << "~> Setting volume ID to " << m_volumeId
                        << " for volume " << v.volumeName());
      v.assignGeometryId(v.geometryId().withVolume(m_volumeId));
    });
  }

  const std::string& name() const override { return m_name; }

 private:
  GeometryIdentifier::Value m_volumeId;
  std::string m_name;
};

}  // namespace

struct GeometryIdentifierBlueprintNodeImpl {
  void add(std::unique_ptr<Configuration> configuration) {
    m_configurations.push_back(std::move(configuration));

    std::vector<std::string> names;
    for (auto& conf : m_configurations) {
      names.push_back(conf->name());
    }

    m_name = boost::algorithm::join(names, ", ");
  }

  std::vector<std::unique_ptr<Configuration>> m_configurations;
  std::string m_name;

  GeometryIdentifierBlueprintNode::CompareVolumes m_sortBy;
};

GeometryIdentifierBlueprintNode::GeometryIdentifierBlueprintNode()
    : m_impl(std::make_unique<GeometryIdentifierBlueprintNodeImpl>()) {}

GeometryIdentifierBlueprintNode::~GeometryIdentifierBlueprintNode() = default;

Volume& GeometryIdentifierBlueprintNode::build(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               const Logger& logger) {
  if (children().size() != 1) {
    throw std::invalid_argument(
        "GeometryIdentifierBlueprintNode must have exactly one child");
  }

  if (m_impl->m_configurations.empty()) {
    throw std::invalid_argument(
        "GeometryIdentifierBlueprintNode has no configuration");
  }

  return children().at(0).build(options, gctx, logger);
}

PortalShellBase& GeometryIdentifierBlueprintNode::connect(
    const BlueprintOptions& options, const GeometryContext& gctx,
    const Logger& logger) {
  return children().at(0).connect(options, gctx, logger);
}

void GeometryIdentifierBlueprintNode::finalize(const BlueprintOptions& options,
                                               const GeometryContext& gctx,
                                               TrackingVolume& parent,
                                               const Logger& logger) {
  ACTS_DEBUG(prefix() << "Finalizing geo id " << name() << " with parent "
                      << parent.volumeName());
  std::set<const TrackingVolume*> previous;
  std::ranges::for_each(parent.volumes(),
                        [&](const auto& v) { previous.insert(&v); });

  children().at(0).finalize(options, gctx, parent, logger);

  std::vector<TrackingVolume*> volumes;
  for (auto& v : parent.volumes()) {
    // Skip volumes that were already in the parent before the subtree
    // was processed
    if (previous.contains(&v)) {
      continue;
    }

    volumes.push_back(&v);
  }

  if (m_impl->m_sortBy) {
    std::ranges::sort(volumes, m_impl->m_sortBy,
                      [](TrackingVolume* v) -> TrackingVolume& { return *v; });
  }

  for (auto* volumePtr : volumes) {
    auto& volume = *volumePtr;
    ACTS_VERBOSE(
        prefix() << " Applying " << m_impl->m_configurations.size()
                 << " geometry ID configuration(s) on subtree starting from "
                 << volume.volumeName());
    for (auto& configuration : m_impl->m_configurations) {
      ACTS_VERBOSE(prefix()
                   << "~> Applying configuration " << configuration->name());
      configuration->apply(prefix(), volume, logger);
    }
    ACTS_DEBUG(prefix() << "~> Final volume ID for " << volume.volumeName()
                        << ": " << volume.geometryId());
  }
}

const std::string& GeometryIdentifierBlueprintNode::name() const {
  return m_impl->m_name;
}

GeometryIdentifierBlueprintNode& GeometryIdentifierBlueprintNode::setLayerIdTo(
    GeometryIdentifier::Value layer) {
  m_impl->add(std::make_unique<FixedLayerConfiguration>(layer));
  return *this;
}

GeometryIdentifierBlueprintNode&
GeometryIdentifierBlueprintNode::incrementLayerIds(
    GeometryIdentifier::Value start) {
  m_impl->add(std::make_unique<IncrementLayerConfiguration>(start));
  return *this;
}

GeometryIdentifierBlueprintNode&
GeometryIdentifierBlueprintNode::setAllVolumeIdsTo(
    GeometryIdentifier::Value volumeId) {
  m_impl->add(std::make_unique<FixedVolumeConfiguration>(volumeId));
  return *this;
}

GeometryIdentifierBlueprintNode& GeometryIdentifierBlueprintNode::sortBy(
    const CompareVolumes& compare) {
  if (!compare) {
    throw std::invalid_argument("Invalid sorting function");
  }
  m_impl->m_sortBy = compare;
  return *this;
}

}  // namespace Acts::Experimental
