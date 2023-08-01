#include "ActsExamples/Utilities/MakeMeasurementParticlesMap.hpp"

using namespace ActsExamples;

ActsExamples::MakeMeasurementParticlesMap::MakeMeasurementParticlesMap(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("MakeMeasurementParticlesMap", lvl), m_cfg(cfg) {
  m_inputHitMap.initialize(m_cfg.inputMeasurementSimhitMap);
  m_inputHits.initialize(m_cfg.inputSimHits);
  m_outputParticleMap.initialize(m_cfg.outputMeasurementParticlesMap);
}

ProcessCode ActsExamples::MakeMeasurementParticlesMap::execute(
    const AlgorithmContext &ctx) const {
  const auto hits = m_inputHits(ctx);
  const auto hitMeasMap = m_inputHitMap(ctx);

  IndexMultimap<ActsFatras::Barcode> outputMap;

  for (const auto &[measIdx, hitIdx] : hitMeasMap) {
    const auto &hit = hits.nth(hitIdx);
    outputMap.emplace(measIdx, hit->particleId());
  }

  m_outputParticleMap(ctx, std::move(outputMap));

  return ProcessCode::SUCCESS;
}
