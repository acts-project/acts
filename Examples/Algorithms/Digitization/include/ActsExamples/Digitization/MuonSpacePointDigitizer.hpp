// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once


#include "Acts/Utilities/Logger.hpp"

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"


namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class MuonSpacePointDigitizer final : public IAlgorithm {
 public:
    struct Config{
      /// @brief Pointer to the tracking geometry to fetch the surfaces
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{};
      /// @brief Name of the input simulated hits collection
      std::string inputSimHits{"simHits"};
      /// @brief Name of the input simulated particles collection
      std::string inputParticles{"particles_simulated"};
      /// @brief Name of the output space points collection
      std::string outputSpacePoints{"space_points"};
      /// @brief Random number generator service
      std::shared_ptr<const RandomNumbers> randomNumbers{};

    };
    /// @brief Constructor  
    MuonSpacePointDigitizer(const Config& cfg, Acts::Logging::Level lvl);
    /// @brief Destructor
    ~MuonSpacePointDigitizer() = default;
    
    /// @brief Initialize the digitizer
    ProcessCode initialize() override;
    /// @brief Execute the digitization 
    ProcessCode execute(const AlgorithmContext& ctx) const override;
    
 private:
    /// @brief Configuration of the digitizer
    Config m_cfg;
    /// @brief Logger for the digitizer
    std::unique_ptr<const Acts::Logger> m_logger;
    /// @brief Data handle for the input simulated hits
   ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};

       ReadDataHandle<SimParticleContainer> m_inputParticles{this,
                                                          "G4Particles"};
      /// @brief Data handle for the output space points
      WriteDataHandle<MuonSpacePointContainer> m_outputSpacePoints{this,"SpacePoints"};
};
}