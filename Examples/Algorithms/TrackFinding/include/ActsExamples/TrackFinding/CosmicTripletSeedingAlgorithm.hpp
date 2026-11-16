#pragma once

#include "Acts/EventData/SeedContainer.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Seeding/CosmicSeedCuts.hpp"
#include "Acts/Seeding/CosmicTripletSeedFinder.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <memory>
#include "Acts/Utilities/Logger.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class CosmicTripletSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSpacePoints;
    std::string outputSeeds;
    // optional, with cuts.requireSameSubdetector or cuts.useCovarianceScale:
    // a vector of CosmicSpacePointClass, one per space point.
    std::string inputSpacePointClasses;
    // optional, std::vector<Acts::GeometryIdentifier>, the surface of each
    // space point. Required if cuts.rejectSameSurface.
    std::string inputSpacePointSurfaces;
    // optional, with cuts.useStripCalibration: a copy of the space points (same indices)
    // with each strip SP used in a seed in the position recalculated with the direction of its best seed.
    // The output seeds always use this copy as their container,
    // so that the CKF parameter estimation starts from the recalculated positions.
    std::string outputSpacePoints;
    // optional, with cuts.useStripCalibration: seeds with two space points (bottom, top),
    // the Athena find2Sp strategy (see Acts::CosmicTripletSeedFinder::createTwoSpacePointSeeds).
    // They index the same container as the three-point seeds. Empty = not produced.
    std::string outputTwoSpacePointSeeds;
    // with outputTwoSpacePointSeeds: doublets only with the SPs that are not used
    // by any three-point seed.
    bool twoSpacePointSeedsFromUnused = false;
    Acts::CosmicSeedCuts cuts;
  };

  // constructor with configuration and optional logger
  CosmicTripletSeedingAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);
  // execute() runs the algorithm: reads the input space points, applies the cosmic triplet seeding finder with the configured cuts,
  // and writes the output seeds and optional recalibrated space points
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  Acts::CosmicTripletSeedFinder m_finder;
  // data handles for input and output
  ReadDataHandle<Acts::SpacePointContainer> m_inputSpacePoints{this, "InputSpacePoints"};
  ReadDataHandle<std::vector<Acts::CosmicSpacePointClass>> m_inputSpacePointClasses{this, "InputSpacePointClasses"};
  ReadDataHandle<std::vector<Acts::GeometryIdentifier>> m_inputSpacePointSurfaces{this, "InputSpacePointSurfaces"};
  WriteDataHandle<Acts::SeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<Acts::SpacePointContainer> m_outputSpacePoints{this, "OutputSpacePoints"};
  WriteDataHandle<Acts::SeedContainer> m_outputTwoSpacePointSeeds{this, "OutputTwoSpacePointSeeds"};
};

}
