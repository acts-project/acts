// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryDelegates.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace Acts {
namespace Experimental {

namespace detail {

/// @brief This is the simplest of all GeometryIdenfier generators
/// it simply counts objects and sets the associated geometry identifiers
///
/// if `processSurfacesPortals` is set to true, it will also
/// recursively set surface and portal geometry ids, more complicated
/// geometry id setting can be done by chaining different prescriptors
class VolumeCounter {
 public:
  bool processSurfacesPortals = false;
  size_t volumeCounter = 0;

  /// @brief Generate geometry Ids for the objects here
  /// @param volume is the volume for which the geometry Ids
  /// are set
  void generateIds(DetectorVolume& volume);
};

/// @brief This emulates a layer like detector, i.e. it keeps the
/// volumeId fixed and increases the layer counting instead of
/// increasing each time the volume id.
///
/// if `processSurfacesPortals` is set to true, it will also
/// recursively set surface and portal geometry ids, more complicated
/// geometry id setting can be done by chaining different prescriptors
///
/// This allows the geometry id in all layer type volumes of a
/// dedicated sub detector have the same volume id but differ
/// in layer id.
class LayerCounter {
 public:
  GeometryIdentifier baseID{0};
  size_t layerCounter = 0;

  /// @brief Generate geometry Ids for the objects here
  /// @param volume is the volume for which the geometry Ids
  /// are set
  void generateIds(DetectorVolume& volume);
};

/// @brief A portal counter per volume
class PortalCounter {
 public:
  /// @brief Generate geometry Ids for the objects here
  /// @param volume is the volume for which the geometry Ids
  /// are set
  void generateIds(DetectorVolume& volume);
};

/// @brief A sensitive counter per volume
class SensitiveCounter {
 public:
  /// @brief Generate geometry Ids for the objects here
  /// @param volume is the volume for which the geometry Ids
  /// are set
  void generateIds(DetectorVolume& volume);
};

/// @brief A passive counter per volume
class PassiveCounter {
 public:
  /// @brief Generate geometry Ids for the objects here
  /// @param volume is the volume for which the geometry Ids
  /// are set
  void generateIds(DetectorVolume& volume);
};

/// @brief checks for duplicate geometry ids and
/// throws exception
class DuplicateIdChecker {
 public:
  // The already 'seen' geometry Ids
  std::unordered_map<GeometryIdentifier, bool> m_geometryIdMap;
  // The already 'seen' volumes
  std::unordered_map<const DetectorVolume*, bool> m_processedVolumes;
  // The already 'seen' portals
  std::unordered_map<const Portal*, bool> m_processedPortals;

  /// @brief Checks the geometry ids of all objects and fills
  /// them into a map, throws an exception if a geometry id is
  /// assigned twice.
  ///
  /// @param volume is the volume for which the geometry Ids
  /// are checked
  void generateIds(DetectorVolume& volume) noexcept(false);
};

/// @brief Checks for unset Geometry Ids and throws exception
/// in case on eis found
class UnsetIdChecker {
 public:
  /// @brief Checks the geometry ids of all objects and fills
  /// them into a map, throws an exception if a geometry id is
  /// assigned twice.
  ///
  /// @param volume is the volume for which the geometry Ids
  /// are checked
  void generateIds(DetectorVolume& volume) noexcept(false);
};

/// Volume restricted generator - this is a helper class that
/// allows to restrict Id Generators to a specific volume
/// by provind a list of detector volume pointers.
///
/// @tparam generator_t the generator type restricted tot his volume(s)
///
/// This can be used, e.g. for a layer id generation for
/// volumes representing different layer volumes
template <typename generator_t>
class VolumeRestrictedIdGenerator {
 public:
  /// Any type of Id generator
  generator_t generator;
  /// The volume restriction for this generator
  std::vector<const DetectorVolume*> restrictedVolumes = {};

  /// @brief  Convenience constructor with a single generator
  ///
  /// @param gen the generator in question
  VolumeRestrictedIdGenerator(
      const generator_t& gen,
      const std::vector<const DetectorVolume*>& rvolumes)
      : generator(gen), restrictedVolumes(rvolumes) {}

  /// A restriced generator, it checks first if the volume is part of the
  /// restriction list
  ///
  /// @param volume is the detector volume for which (and for its sub objects) the detector
  /// geometry id is generated and set
  ///
  void generateIds(DetectorVolume& volume) {
    if (std::find(restrictedVolumes.begin(), restrictedVolumes.end(),
                  &volume) != restrictedVolumes.end()) {
      generator.generateIds(volume);
    }
  }
};

/// String restricted generator - this is a helper class that
/// allows to restrict Id Generators to volumes having/containing
/// a certain name string.
///
/// @tparam generator_t the generator type restricted tot his volume(s)
///
/// This can be used, e.g. for a layer id generation for
/// volumes representing different layer volumes
template <typename generator_t>
class NameRestrictedIdGenerator {
 public:
  /// Any type of Id generator
  generator_t generator;
  /// The name restriction for this generator
  std::string restrictedName;
  /// Indicate if you need an exact match or a contains
  bool exactMatch = false;

  /// @brief  Convenience constructor with a single generator
  ///
  /// @param gen the generator in question
  NameRestrictedIdGenerator(const generator_t& gen, const std::string& rname)
      : generator(gen), restrictedName(rname) {}

  /// A restriced generator, it checks first if the volume is part of the
  /// restriction list
  ///
  /// @param volume is the detector volume for which (and for its sub objects) the detector
  /// geometry id is generated and set
  ///
  void generateIds(DetectorVolume& volume) {
    std::string volumeName = volume.name();
    if ((exactMatch and volumeName == restrictedName) or
        (not exactMatch and
         volumeName.find(restrictedName) != std::string::npos)) {
      generator.generateIds(volume);
    }
  }
};

/// This is a helper class for chained generators, it allows
/// to set up a set of rules how to generate goeometry Ids
///
/// @tparam generators_t the generators that will be called in sequence
///
template <typename... generators_t>
class ChainedGeometryIdGenerator {
 public:
  // The stored generators
  std::tuple<generators_t...> generators;

  /// Constructor for chained generators in a tuple, this will unroll
  /// the tuple and call them in sequence
  ///
  /// @param gens the generators to be called in chain
  ChainedGeometryIdGenerator(const std::tuple<generators_t...>& gens)
      : generators(gens) {}

  /// A combined chain generator
  ///
  /// @param volume is the detector volume for which (and for its sub objects) the detector
  /// geometry id is generated and set
  ///
  void generateIds(DetectorVolume& volume) {
    // Unfold the tuple and add the attachers
    std::apply(
        [&](auto&&... generator) { ((generator.generateIds(volume)), ...); },
        generators);
  }
};

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
