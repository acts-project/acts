#include "ActsExamples/Utilities/MakeMeasurmentParticlesMap.hpp"

using namespace ActsExamples;

ActsExamples::MakeMeasurmentParticlesMap::MakeMeasurmentParticlesMap(
    Config cfg, Acts::Logging::Level lvl)
    : IAlgorithm("MakeMeasurmentParticlesMap", lvl), m_cfg(cfg) {}

ProcessCode ActsExamples::MakeMeasurmentParticlesMap::execute(
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
