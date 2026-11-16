#include "ActsExamples/TrackFinding/SimHitSpacePointMaker.hpp"

#include "Acts/EventData/SpacePointColumns.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace ActsExamples {

namespace {

bool contains(const std::vector<std::uint32_t>& v, std::uint32_t x) {
  return std::find(v.begin(), v.end(), x) != v.end();
}

bool hasVolumeLists(const SimHitSpacePointMaker::Config& cfg) {
  return !cfg.pixelBarrelVolumes.empty() || !cfg.pixelEndcapVolumes.empty() ||
         !cfg.stripBarrelVolumes.empty() || !cfg.stripEndcapVolumes.empty();
}

// classify a hit by its volume ID, using the lists in the config (returns false if the volume is not in any list)
bool classByVolume(const SimHitSpacePointMaker::Config& cfg,
                   std::uint32_t volume, Acts::CosmicSpacePointClass& out) {
  using Acts::CosmicSpacePointClass;
  if (contains(cfg.pixelBarrelVolumes, volume)) {
    out = CosmicSpacePointClass::BarrelPixel;
  } else if (contains(cfg.pixelEndcapVolumes, volume)) {
    out = CosmicSpacePointClass::EndcapPixel;
  } else if (contains(cfg.stripBarrelVolumes, volume)) {
    out = CosmicSpacePointClass::BarrelStrip;
  } else if (contains(cfg.stripEndcapVolumes, volume)) {
    out = CosmicSpacePointClass::EndcapStrip;
  } else {
    return false;
  }
  return true;
}

}

SimHitSpacePointMaker::SimHitSpacePointMaker(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("SimHitSpacePointMaker", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input sim hit collection");
  }
  if (m_cfg.outputSpacePoints.empty()) {
    throw std::invalid_argument("Missing output space point collection");
  }
  if (!(m_cfg.sigmaR >= 0.0) || !(m_cfg.sigmaZ >= 0.0)) {
    throw std::invalid_argument("sigmaR and sigmaZ must be non-negative");
  }

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_outputSpacePoints.initialize(m_cfg.outputSpacePoints);
  // optional output: particle IDs of the sim hits
  if (!m_cfg.outputParticleIds.empty()) {
    m_outputParticleIds.initialize(m_cfg.outputParticleIds);
  }
  // optional output: classes of space points for CosmicTripletSeedFinder
  if (!m_cfg.outputSpacePointClasses.empty()) {
    if (hasVolumeLists(m_cfg)) {
      // a volume cannot be in two different lists
      std::map<std::uint32_t, int> listOf;
      const std::array<const std::vector<std::uint32_t>*, 4> lists = {
          &m_cfg.pixelBarrelVolumes, &m_cfg.pixelEndcapVolumes,
          &m_cfg.stripBarrelVolumes, &m_cfg.stripEndcapVolumes};
      for (int i = 0; i < 4; ++i) {
        for (std::uint32_t v : *lists[i]) {
          const auto it = listOf.emplace(v, i).first;
          if (it->second != i) {
            throw std::invalid_argument(
                "volume " + std::to_string(v) +
                " appears in more than one class list");
          }
        }
      }
    } else if (!std::isfinite(m_cfg.stripMinRadius) ||
               !std::isfinite(m_cfg.pixelBarrelMaxAbsZ) ||
               !std::isfinite(m_cfg.stripBarrelMaxAbsZ)) {
      throw std::invalid_argument(
          "outputSpacePointClasses needs either the four *Volumes lists or "
          "stripMinRadius, pixelBarrelMaxAbsZ and stripBarrelMaxAbsZ");
    }
    m_outputSpacePointClasses.initialize(m_cfg.outputSpacePointClasses);
  }
  // optional output: surfaces of each space point (rejectSameSurface)
  if (!m_cfg.outputSpacePointSurfaces.empty()) {
    m_outputSpacePointSurfaces.initialize(m_cfg.outputSpacePointSurfaces);
  }
}

Acts::SpacePointContainer SimHitSpacePointMaker::makeSpacePoints(
    const SimHitContainer& simHits, double sigmaR, double sigmaZ,
    std::vector<SimBarcode>* particleIds) {
  using Acts::SpacePointColumns;

  // all columns that the seeders read (x, y, z, r, variances) plus time
  Acts::SpacePointContainer spacePoints(
      SpacePointColumns::SourceLinks | SpacePointColumns::X |
      SpacePointColumns::Y | SpacePointColumns::Z | SpacePointColumns::R |
      SpacePointColumns::Time | SpacePointColumns::VarianceZ |
      SpacePointColumns::VarianceR);
  spacePoints.reserve(static_cast<std::uint32_t>(simHits.size()));

  if (particleIds != nullptr) {
    particleIds->clear();
    particleIds->reserve(simHits.size());
  }

  const auto varianceR = static_cast<float>(sigmaR * sigmaR);
  const auto varianceZ = static_cast<float>(sigmaZ * sigmaZ);

  for (const SimHit& hit : simHits) {
    const Acts::Vector3 pos = hit.position();

    auto sp = spacePoints.createSpacePoint();
    sp.x() = static_cast<float>(pos.x());
    sp.y() = static_cast<float>(pos.y());
    sp.z() = static_cast<float>(pos.z());
    sp.r() = static_cast<float>(std::hypot(pos.x(), pos.y()));
    sp.time() = static_cast<float>(hit.time());
    sp.varianceR() = varianceR;
    sp.varianceZ() = varianceZ;

    if (particleIds != nullptr) {
      particleIds->push_back(hit.particleId());
    }
  }

  return spacePoints;
}

Acts::CosmicSpacePointClass SimHitSpacePointMaker::classifyHit(
    const Config& cfg, const SimHit& hit) {
  using Acts::CosmicSpacePointClass;

  if (hasVolumeLists(cfg)) {
    const auto volume = static_cast<std::uint32_t>(hit.geometryId().volume());
    CosmicSpacePointClass cls = CosmicSpacePointClass::BarrelPixel;
    if (!classByVolume(cfg, volume, cls)) {
      throw std::runtime_error("SimHit volume " + std::to_string(volume) +
                               " is not in any class list");
    }
    return cls;
  }

  const Acts::Vector3 pos = hit.position();
  const bool strip = std::hypot(pos.x(), pos.y()) >= cfg.stripMinRadius;
  const bool barrel = std::abs(pos.z()) < (strip ? cfg.stripBarrelMaxAbsZ
                                                 : cfg.pixelBarrelMaxAbsZ);
  if (strip) {
    return barrel ? CosmicSpacePointClass::BarrelStrip
                  : CosmicSpacePointClass::EndcapStrip;
  }
  return barrel ? CosmicSpacePointClass::BarrelPixel
                : CosmicSpacePointClass::EndcapPixel;
}

std::vector<Acts::CosmicSpacePointClass> SimHitSpacePointMaker::classifyHits(
    const Config& cfg, const SimHitContainer& simHits) {
  std::vector<Acts::CosmicSpacePointClass> classes;
  classes.reserve(simHits.size());
  // volumes lists are used to classify hits; if a volume is not in any list, it is an error
  std::set<std::uint32_t> unmapped;
  const bool byVolume = hasVolumeLists(cfg);
  for (const SimHit& hit : simHits) {
    if (byVolume) {
      const auto volume =
          static_cast<std::uint32_t>(hit.geometryId().volume());
      Acts::CosmicSpacePointClass cls = Acts::CosmicSpacePointClass::BarrelPixel;
      if (!classByVolume(cfg, volume, cls)) {
        unmapped.insert(volume);
      }
      classes.push_back(cls);
    } else {
      classes.push_back(classifyHit(cfg, hit));
    }
  }
  if (!unmapped.empty()) {
    std::string list;
    for (std::uint32_t v : unmapped) {
      list += (list.empty() ? "" : ", ") + std::to_string(v);
    }
    throw std::runtime_error("SimHit volumes not in any class list: " + list);
  }
  return classes;
}

ProcessCode SimHitSpacePointMaker::execute(const AlgorithmContext& ctx) const {
  const SimHitContainer& simHits = m_inputSimHits(ctx);

  // classify hits if the outputSpacePointClasses handle is initialized
  std::vector<Acts::CosmicSpacePointClass> classes;
  if (m_outputSpacePointClasses.isInitialized()) {
    try {
      classes = classifyHits(m_cfg, simHits);
    } catch (const std::exception& e) {
      ACTS_ERROR("SimHitSpacePointMaker: " << e.what());
      return ProcessCode::ABORT;
    }
    std::array<std::size_t, 4> counts{};
    for (const auto cls : classes) {
      ++counts[static_cast<std::size_t>(cls)];
    }
    ACTS_DEBUG("SimHitSpacePointMaker classes: pixel barrel=" << counts[0] << " pixel endcap=" << counts[1] << " strip barrel=" << counts[2] << " strip endcap=" << counts[3]);
  }

  std::vector<SimBarcode> particleIds;
  Acts::SpacePointContainer spacePoints = makeSpacePoints(
      simHits, m_cfg.sigmaR, m_cfg.sigmaZ,
      m_outputParticleIds.isInitialized() ? &particleIds : nullptr);

  ACTS_DEBUG("SimHitSpacePointMaker: " << simHits.size() << " sim hits -> " << spacePoints.size() << " space points");

  m_outputSpacePoints(ctx, std::move(spacePoints));
  if (m_outputParticleIds.isInitialized()) {
    m_outputParticleIds(ctx, std::move(particleIds));
  }
  if (m_outputSpacePointClasses.isInitialized()) {
    m_outputSpacePointClasses(ctx, std::move(classes));
  }
  if (m_outputSpacePointSurfaces.isInitialized()) {
    // space point index it is the position of the hit in the SimHitContainer
    std::vector<Acts::GeometryIdentifier> surfaces;
    surfaces.reserve(simHits.size());
    for (const SimHit& hit : simHits) {
      surfaces.push_back(hit.geometryId());
    }
    m_outputSpacePointSurfaces(ctx, std::move(surfaces));
  }

  return ProcessCode::SUCCESS;
}

}