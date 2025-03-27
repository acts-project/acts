

#include "Acts/Seeding/StrawChamberLineSeeder.hpp"

#include "ActsExamples/TrackFinding/MuonSegmentFinder.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"


namespace ActsExamples{



    MuonSegmentFinder::MuonSegmentFinder(Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("MuonSegmentFinder", lvl),
      m_cfg(std::move(cfg)),
      m_logger(Acts::getDefaultLogger("MuonSegmentFinder", lvl)) {

    if (m_cfg.inHoughSeeds.empty()) {
      throw std::invalid_argument("MuonSegmentFinder: Missing maxima collection");
    }
    if (m_cfg.outSegments.empty()) {
      throw std::invalid_argument("MuonSegmentFinder: Missing out segment container");
    }

}

ProcessCode MuonSegmentFinder::execute( const AlgorithmContext& ctx) const {
  
  const MuonHoughMaxContainer& gotMaxima = m_inputMax(ctx);
  
  using UnCalibSpCont_t = MuonSpacePointCalibrator::UnCalibSpVec_t;
  using CalibSpCont_t = MuonSpacePointCalibrator::CalibSpCont_t;

  using StrawSeeder_t = Acts::StrawChamberLineSeeder<UnCalibSpCont_t, MuonSpacePointSorter, CalibSpCont_t,
                                                     MuonSpacePointCalibrator>;
  using SeedCandidate_t = StrawSeeder_t::DriftCircleSeed;
  
  
  Acts::CalibrationContext calibCtx{ctx};
  for (const MuonHoughMaximum& max : gotMaxima){
    StrawSeeder_t::Config seedCfg{};
    StrawSeeder_t seeder{max.hits(), std::move(seedCfg), m_logger->clone()};
    while (std::optional<SeedCandidate_t> cand = seeder.generateSeed(calibCtx));
  }
  MuonSegmentContainer outSegments{};


  m_outSegments(ctx, std::move(outSegments));

  return ProcessCode::SUCCESS;
}

ProcessCode MuonSegmentFinder::initialize() {
  m_inputMax.initialize(m_cfg.inHoughSeeds);
  m_outSegments.initialize(m_cfg.outSegments);
  return ProcessCode::SUCCESS;
}
}