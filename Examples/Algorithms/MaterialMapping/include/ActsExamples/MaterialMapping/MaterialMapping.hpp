// This file is part of the Acts project.
//
// Copyright (C) 2017-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/Material/interface/IMaterialMapper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/MaterialMapping/IMaterialWriter.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ActsExamples {
class IMaterialWriter;
struct AlgorithmContext;
}  // namespace ActsExamples

namespace Acts {

class TrackingGeometry;
class ISurfaceMaterial;
class IVolumeMaterial;

using SurfaceMaterialMap =
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;

using VolumeMaterialMap =
    std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;

using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;

namespace Experimental {
class Detector;
}

}  // namespace Acts

namespace ActsExamples {

/// @class MaterialMapping
///
/// @brief Initiates and executes material mapping
///
/// The MaterialMapping reads in the MaterialTrack with a dedicated
/// reader and uses the material mapper to project the material onto
/// either surfaces or volumes of the input geometry
///
class MaterialMapping final : public IAlgorithm {
 public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config {
    // Geometry context for the state creation
    std::reference_wrapper<const Acts::GeometryContext> geoContext;

    // MagneticField  context for the state creation
    std::reference_wrapper<const Acts::MagneticFieldContext> magFieldContext;

    /// Input collection
    std::string collection = "material_tracks";

    /// The material collection to be stored - mapped
    std::string mappedCollection = "mapped_material_tracks";

    /// The material collection to be stored - unmapped
    std::string unmappedCollection = "unmapped_material_tracks";

    /// The material mappers - they are executed in a sequential order,
    /// i.e. the subsequent mapper runs only on the remaining hits of the
    /// prior mappper
    std::vector<std::shared_ptr<Acts::IMaterialMapper>> materialMappers = {};

    /// The writer of the material
    std::vector<std::shared_ptr<IMaterialWriter>> materialWriters{};

    /// The TrackingGeometry to be mapped on
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;

    /// The Detector to be mapped on
    std::shared_ptr<const Acts::Experimental::Detector> detector = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  MaterialMapping(const Config& cfg,
                  Acts::Logging::Level level = Acts::Logging::INFO);

  /// Destructor
  ///
  /// @note the destructor also invokes the material writer(s)
  ~MaterialMapping() final;

  /// Framework execute method
  ///
  /// @param context The algorithm context for event consistency
  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;  //!< internal config object

  std::vector<std::unique_ptr<Acts::MaterialMappingState>>
      m_mappingStates;  //!< Material mapping states, one for each mapper

  ReadDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};
  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_mappedMaterialTracks{this, "MappedMaterialTracks"};
  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_unmappedMaterialTracks{this, "UnmappedMaterialTracks"};
};

}  // namespace ActsExamples
