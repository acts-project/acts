// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/CosmicSeedCuts.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

class SimHitSpacePointMaker final : public IAlgorithm {
 public:
  struct Config {
    // sim hits input collection (ActsExamples::SimHitContainer)
    std::string inputSimHits;
    // sim hits output collection (Acts::SpacePointContainer)
    std::string outputSpacePoints;
    // optional: if not empty, writes a std::vector<SimBarcode> with the barcode
    // of each sim hit
    //  (same index as the space point)
    std::string outputParticleIds;
    // assigned standard deviation to r [mm]; only defines the variance
    double sigmaR = 0.1 * Acts::UnitConstants::mm;
    // assigned standard deviation to z [mm]; only defines the variance
    double sigmaZ = 0.1 * Acts::UnitConstants::mm;
    // optional: if not empty, writes a std::vector<Acts::CosmicSpacePointClass>
    // with the class of each space point (same index as the space point)
    std::string outputSpacePointClasses;
    // optional: if not empty, writes a std::vector<Acts::GeometryIdentifier>
    // with the geometryID() of SimHit of each space point
    std::string outputSpacePointSurfaces;

    // volume classification
    std::vector<std::uint32_t> pixelBarrelVolumes;
    std::vector<std::uint32_t> pixelEndcapVolumes;
    std::vector<std::uint32_t> stripBarrelVolumes;
    std::vector<std::uint32_t> stripEndcapVolumes;

    // geometric classification
    double stripMinRadius = std::numeric_limits<double>::quiet_NaN();
    double pixelBarrelMaxAbsZ = std::numeric_limits<double>::quiet_NaN();
    double stripBarrelMaxAbsZ = std::numeric_limits<double>::quiet_NaN();
  };

  explicit SimHitSpacePointMaker(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }
  // conversion
  /// @param simHits     input hits
  /// @param sigmaR      standard deviation assigned to r [mm]
  /// @param sigmaZ      standard deviation assigned to z [mm]
  /// @param particleIds if not nullptr, receives the barcode for each hit
  /// @return container with one space point per hit
  static Acts::SpacePointContainer makeSpacePoints(
      const SimHitContainer& simHits, double sigmaR, double sigmaZ,
      std::vector<SimBarcode>* particleIds = nullptr);

  // class of a SimHit
  static Acts::CosmicSpacePointClass classifyHit(const Config& cfg,
                                                 const SimHit& hit);

  // class of each hit
  static std::vector<Acts::CosmicSpacePointClass> classifyHits(
      const Config& cfg, const SimHitContainer& simHits);

 private:
  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  WriteDataHandle<Acts::SpacePointContainer> m_outputSpacePoints{
      this, "OutputSpacePoints"};
  WriteDataHandle<std::vector<SimBarcode>> m_outputParticleIds{
      this, "OutputParticleIds"};
  WriteDataHandle<std::vector<Acts::CosmicSpacePointClass>>
      m_outputSpacePointClasses{this, "OutputSpacePointClasses"};
  WriteDataHandle<std::vector<Acts::GeometryIdentifier>>
      m_outputSpacePointSurfaces{this, "OutputSpacePointSurfaces"};
};

}  // namespace ActsExamples
