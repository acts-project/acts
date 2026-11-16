#include "ActsExamples/TrackFinding/CosmicTripletSeedingAlgorithm.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <memory>
#include "Acts/Utilities/Logger.hpp"

namespace ActsExamples {

CosmicTripletSeedingAlgorithm::CosmicTripletSeedingAlgorithm(Config cfg,std::unique_ptr<const Acts::Logger> logger):
                                IAlgorithm("CosmicTripletSeedingAlgorithm", std::move(logger)), m_cfg(std::move(cfg)),m_finder(m_cfg.cuts) {
  
  // requirements on the configuration
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing input space point collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output seed collection");
  }
  if ((m_cfg.cuts.requireSameSubdetector || m_cfg.cuts.useCovarianceScale || m_cfg.cuts.useSubdetectorAngularCuts) &&
      m_cfg.inputSpacePointClasses.empty()) {
    throw std::invalid_argument(
        "cuts.requireSameSubdetector / cuts.useCovarianceScale / "
        "cuts.useSubdetectorAngularCuts require inputSpacePointClasses");
  }
  if (m_cfg.cuts.rejectSameSurface && m_cfg.inputSpacePointSurfaces.empty()) {
    throw std::invalid_argument(
        "cuts.rejectSameSurface requires inputSpacePointSurfaces");
  }
  if (!m_cfg.outputSpacePoints.empty() && !m_cfg.cuts.useStripCalibration) {
    throw std::invalid_argument(
        "outputSpacePoints requires cuts.useStripCalibration");
  }
  if (m_cfg.twoSpacePointSeedsFromUnused &&
      m_cfg.outputTwoSpacePointSeeds.empty()) {
    throw std::invalid_argument(
        "twoSpacePointSeedsFromUnused requires outputTwoSpacePointSeeds");
  }

  m_inputSpacePoints.initialize(m_cfg.inputSpacePoints);
  // handle opciotional, which initialize only if the corresponding input is not empty
  if (!m_cfg.inputSpacePointClasses.empty()) {
    m_inputSpacePointClasses.initialize(m_cfg.inputSpacePointClasses);
  }
  if (!m_cfg.inputSpacePointSurfaces.empty()) {
    m_inputSpacePointSurfaces.initialize(m_cfg.inputSpacePointSurfaces);
  }
  m_outputSeeds.initialize(m_cfg.outputSeeds);
  if (!m_cfg.outputSpacePoints.empty()) {
    m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  }
  if (!m_cfg.outputTwoSpacePointSeeds.empty()) {
    m_outputTwoSpacePointSeeds.initialize(m_cfg.outputTwoSpacePointSeeds);
  }
}

ProcessCode CosmicTripletSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  const Acts::SpacePointContainer& spacePoints = m_inputSpacePoints(ctx);

  Acts::SeedContainer seeds;

  const std::vector<Acts::CosmicSpacePointClass>* classes = nullptr;
  if (m_inputSpacePointClasses.isInitialized()) {
    classes = &m_inputSpacePointClasses(ctx);
  }
  const std::vector<Acts::GeometryIdentifier>* surfaces = nullptr;
  if (m_inputSpacePointSurfaces.isInitialized()) {
    surfaces = &m_inputSpacePointSurfaces(ctx);
  }
  // if the input space points have strip calibration details,
  // we can use them to recalibrate the positions of strip space points in the seeds
  std::shared_ptr<Acts::SpacePointContainer> calibrated;
  if (m_cfg.cuts.useStripCalibration && spacePoints.hasColumns(Acts::SpacePointColumns::StripCalibrationDetails)) {
    calibrated = std::make_shared<Acts::SpacePointContainer>(spacePoints);
    seeds.assignSpacePointContainer(std::shared_ptr<const Acts::SpacePointContainer>(calibrated));
  }

  // createSeeds() generates triplet seeds from the given space points,
  // applying the configured cuts and storing the results in outputSeeds
  try {m_finder.createSeeds(spacePoints, seeds, classes, calibrated.get(),surfaces);
  } catch (const std::invalid_argument& e) {
    ACTS_ERROR("CosmicTripletSeedingAlgorithm: " << e.what());
    return ProcessCode::ABORT;
  }

  ACTS_DEBUG("CosmicTripletSeedingAlgorithm: " << spacePoints.size() << " space points -> "<< seeds.size() << " seeds");

  if (m_outputTwoSpacePointSeeds.isInitialized()) {
    Acts::SeedContainer doublets;
    if (calibrated != nullptr) {
      doublets.assignSpacePointContainer(
          std::shared_ptr<const Acts::SpacePointContainer>(calibrated));
    }
    // SPs used by any triplet seed are marked as used,
    // and only the unused SPs are considered for doublet seeds
    std::vector<std::uint8_t> used;
    if (m_cfg.twoSpacePointSeedsFromUnused) {
      used.assign(spacePoints.size(), 0);
      const Acts::SeedContainer& triplets = seeds;
      for (Acts::ConstSeedProxy seed : triplets) {
        for (Acts::SpacePointIndex i : seed.spacePointIndices()) {
          used[i] = 1;
        }
      }
    }
    try {
      m_finder.createTwoSpacePointSeeds(spacePoints, doublets, classes,
                                        used.empty() ? nullptr : &used,
                                        surfaces);
    } catch (const std::invalid_argument& e) {
      ACTS_ERROR("CosmicTripletSeedingAlgorithm: " << e.what());
      return ProcessCode::ABORT;
    }
    ACTS_DEBUG("CosmicTripletSeedingAlgorithm: " << doublets.size() << " two space point seeds");
    m_outputTwoSpacePointSeeds(ctx, std::move(doublets));
  }

  if (m_outputSpacePoints.isInitialized()) {
    m_outputSpacePoints(ctx, calibrated != nullptr
                                 ? Acts::SpacePointContainer(*calibrated)
                                 : Acts::SpacePointContainer(spacePoints));
  }
  m_outputSeeds(ctx, std::move(seeds));
  return ProcessCode::SUCCESS;
}

}