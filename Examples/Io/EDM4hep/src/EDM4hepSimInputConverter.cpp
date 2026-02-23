// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepSimInputConverter.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <map>
#include <numeric>
#include <ranges>
#include <stdexcept>

#include <DD4hep/Detector.h>
#include <boost/histogram.hpp>
#include <boost/histogram/ostream.hpp>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>
#include <podio/ObjectID.h>

namespace bh = boost::histogram;

namespace ActsExamples {
namespace detail {
struct ParticleInfo {
  std::size_t particleIndex;
  // std::uint16_t numHits;
};

}  // namespace detail

EDM4hepSimInputConverter::EDM4hepSimInputConverter(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : PodioInputConverter("EDM4hepSimInputConverter", config.inputFrame,
                          std::move(logger)),
      m_cfg(config) {
  if (m_cfg.outputParticlesGenerator.empty()) {
    throw std::invalid_argument(
        "Missing output collection generator particles");
  }
  if (m_cfg.outputParticlesSimulation.empty()) {
    throw std::invalid_argument(
        "Missing output collection simulated particles");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing output collection sim hits");
  }

  if (m_cfg.outputSimVertices.empty()) {
    throw std::invalid_argument("Missing output collection sim vertices");
  }

  m_outputParticlesGenerator.initialize(m_cfg.outputParticlesGenerator);
  m_outputParticlesSimulation.initialize(m_cfg.outputParticlesSimulation);
  m_outputSimHits.initialize(m_cfg.outputSimHits);
  m_outputSimHitAssociation.maybeInitialize(m_cfg.outputSimHitAssociation);
  m_outputSimVertices.initialize(m_cfg.outputSimVertices);

  ACTS_LOG_WITH_LOGGER(*m_logger, Acts::Logging::INFO,
                       "Configured EDM4hepSimInputConverter:");
  auto printCut = [](std::optional<double> opt) {
    if (opt.has_value()) {
      return std::to_string(opt.value());
    } else {
      return std::string{"<none>"};
    }
  };
  ACTS_LOG_WITH_LOGGER(*m_logger, Acts::Logging::INFO,
                       "- particle r: [" << printCut(m_cfg.particleRMin) << ", "
                                         << printCut(m_cfg.particleRMax)
                                         << "] mm");
  ACTS_LOG_WITH_LOGGER(*m_logger, Acts::Logging::INFO,
                       "- particle z: [" << printCut(m_cfg.particleZMin) << ", "
                                         << printCut(m_cfg.particleZMax)
                                         << "] mm");

  m_cfg.trackingGeometry->visitSurfaces([&](const auto* surface) {
    const auto* detElement =
        dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
            surface->surfacePlacement());

    if (detElement == nullptr) {
      ACTS_LOG_WITH_LOGGER(*m_logger, Acts::Logging::ERROR,
                           "Surface has no associated detector element");
      return;
    }

    const auto translation = detElement->sourceElement()
                                 .nominal()
                                 .worldTransformation()
                                 .GetTranslation();
    Acts::Vector3 position;
    position << translation[0], translation[1], translation[2];
    position *= Acts::UnitConstants::cm;

    m_surfaceMap.insert({detElement->sourceElement().key(), surface});
  });
}

bool EDM4hepSimInputConverter::acceptParticle(
    const edm4hep::MCParticle& particle) const {
  double z = particle.getVertex().z * Acts::UnitConstants::mm;
  if (m_cfg.particleZMin.has_value() && z < m_cfg.particleZMin.value()) {
    ACTS_VERBOSE("Skipping particle with z=" << z << " < zMin="
                                             << m_cfg.particleZMin.value());
    return false;
  }

  if (m_cfg.particleZMax.has_value() && z > m_cfg.particleZMax.value()) {
    ACTS_VERBOSE("Skipping particle with z=" << z << " > zMax="
                                             << m_cfg.particleZMax.value());
    return false;
  }

  double r = std::hypot(particle.getVertex()[0], particle.getVertex()[1]) *
             Acts::UnitConstants::mm;

  if (m_cfg.particleRMin.has_value() && r < m_cfg.particleRMin.value()) {
    ACTS_VERBOSE("Skipping particle with r=" << r << " < rMin="
                                             << m_cfg.particleRMin.value());
    return false;
  }

  if (m_cfg.particleRMax.has_value() && r > m_cfg.particleRMax.value()) {
    ACTS_VERBOSE("Skipping particle with r=" << r << " > rMax="
                                             << m_cfg.particleRMax.value());
    return false;
  }

  double pt = std::hypot(particle.getMomentum()[0], particle.getMomentum()[1]);

  if (m_cfg.particlePtMin.has_value() && pt < m_cfg.particlePtMin.value()) {
    ACTS_VERBOSE("Skipping particle with pt=" << pt << " < ptMin="
                                              << m_cfg.particlePtMin.value());
    return false;
  }

  if (m_cfg.particlePtMax.has_value() && pt > m_cfg.particlePtMax.value()) {
    ACTS_VERBOSE("Skipping particle with pt=" << pt << " > ptMax="
                                              << m_cfg.particlePtMax.value());
    return false;
  }

  return true;
}

bool EDM4hepSimInputConverter::particleOrDescendantsHaveHits(
    const edm4hep::MCParticle& particle,
    const std::function<std::uint8_t(const edm4hep::MCParticle&)>& getNumHits)
    const {
  // Check if this particle has hits
  if (getNumHits(particle) > 0) {
    return true;
  }

  // Recursively check all descendants
  for (const auto& daughter : particle.getDaughters()) {
    if (particleOrDescendantsHaveHits(daughter, getNumHits)) {
      return true;
    }
  }

  return false;
}

namespace {
bool isGeneratorStable(const edm4hep::MCParticle& particle) {
  // https://arxiv.org/pdf/1912.08005#subsection.1.A.1
  constexpr int kUndecayedPhysicalParticleStatus = 1;
  constexpr int kDecayedPhysicalParticleStatus = 2;
  int status = particle.getGeneratorStatus();
  return status == kUndecayedPhysicalParticleStatus ||
         status == kDecayedPhysicalParticleStatus;
}

void findGeneratorStableParticles(
    const edm4hep::MCParticle& particle,
    std::vector<edm4hep::MCParticle>& outputParticles) {
  // Theoretically we shouldn't get here because we shouldn't descend this far
  if (particle.isCreatedInSimulation()) {
    return;
  }

  if (isGeneratorStable(particle)) {
    // This is a generator stable particle, record and do not descend further
    if (std::ranges::find(outputParticles, particle) == outputParticles.end()) {
      outputParticles.push_back(particle);
    }
    return;
  }

  for (const auto& daughter : particle.getDaughters()) {
    findGeneratorStableParticles(daughter, outputParticles);
  }
}
}  // namespace

ProcessCode EDM4hepSimInputConverter::convert(const AlgorithmContext& ctx,
                                              const podio::Frame& frame) const {
  ACTS_DEBUG("Reading EDM4hep inputs");

  const auto& mcParticleCollection =
      frame.get<edm4hep::MCParticleCollection>(m_cfg.inputParticles);

  ACTS_DEBUG("Total input particles: " << mcParticleCollection.size()
                                       << " particles");

  std::vector<SimBarcode> unorderedParticlesInitial;

  // Read particles from the input file
  // Find particles without parents and group them by vtx position to find
  // primary vertices
  std::vector<std::pair<Acts::Vector3, std::vector<int>>> primaryVertices;

  std::size_t nVertexParticles = 0;

  {
    Acts::ScopedTimer timer("Finding primary vertices", logger(),
                            Acts::Logging::DEBUG);
    for (const auto& particle : mcParticleCollection) {
      if (particle.parents_size() > 0) {
        // not a primary vertex
        continue;
      }
      const auto& vtx = particle.getEndpoint();

      // @TODO: Might have to use the time here as well
      Acts::Vector3 vtxPos = {vtx[0], vtx[1], vtx[2]};
      vtxPos *= Acts::UnitConstants::mm;

      // linear search for vector
      auto it = std::ranges::find_if(primaryVertices, [&vtxPos](const auto& v) {
        return v.first == vtxPos;
      });

      std::vector<int>* vertexParticles = nullptr;

      if (it == primaryVertices.end()) {
        ACTS_VERBOSE("Found primary vertex at " << vtx.x << ", " << vtx.y
                                                << ", " << vtx.z);
        primaryVertices.emplace_back(vtxPos, std::vector<int>{});
        vertexParticles = &primaryVertices.back().second;
      } else {
        vertexParticles = &it->second;
      }

      assert(vertexParticles != nullptr);

      if (isGeneratorStable(particle)) {
        vertexParticles->push_back(particle.getObjectID().index);
        nVertexParticles += 1;
      }

      nVertexParticles += particle.getDaughters().size();
      std::ranges::copy(
          particle.getDaughters() | std::views::transform([](const auto& p) {
            return p.getObjectID().index;
          }),
          std::back_inserter(*vertexParticles));
    }
  }

  ACTS_DEBUG("Found " << primaryVertices.size() << " primary vertices with "
                      << nVertexParticles << " outgoing particles in total");

  // key: child, value: parent
  ParentRelationship parentRelationship;

  // key: input particle index, value: index in the unordered particle
  // container
  std::unordered_map<int, detail::ParticleInfo> edm4hepParticleMap;

  std::vector<std::uint16_t> numSimHits;
  numSimHits.resize(mcParticleCollection.size());

  std::size_t nGeneratorParticles = 0;
  for (const auto& particle : mcParticleCollection) {
    if (!particle.isCreatedInSimulation()) {
      nGeneratorParticles += 1;
    }
  }

  std::vector<const edm4hep::SimTrackerHitCollection*> simHitCollections;
  for (const auto& name : m_cfg.inputSimHits) {
    simHitCollections.push_back(
        &frame.get<edm4hep::SimTrackerHitCollection>(name));
  }

  // Let's figure out first how many hits each particle has:
  for (const auto* inputHits : simHitCollections) {
    for (const auto& hit : *inputHits) {
      auto particle = ActsPlugins::EDM4hepUtil::getParticle(hit);

      std::size_t index = particle.getObjectID().index;

      auto& num = numSimHits.at(index);
      constexpr unsigned int maxNum =
          (1 << (sizeof(decltype(numSimHits)::value_type) * 8)) - 1;

      if (num == maxNum) {
        throw std::runtime_error{"Hit count " + std::to_string(num) +
                                 " is at the limit of " +
                                 std::to_string(maxNum)};
      }

      num += 1;
    }
  }

  std::function<std::uint16_t(const edm4hep::MCParticle&)> getNumHits =
      [&numSimHits](const edm4hep::MCParticle& p) {
        return numSimHits.at(p.getObjectID().index);
      };

  std::optional particlesGeneratorUnordered = std::vector<SimParticle>{};

  std::size_t nPrimaryVertices = 0;
  // Walk the particle tree
  {
    Acts::ScopedTimer timer("Walking particle tree", logger(),
                            Acts::Logging::DEBUG);

    Acts::AveragingScopedTimer timerFindStable(
        "Finding generator-stable particles", logger(), Acts::Logging::DEBUG);

    std::vector<edm4hep::MCParticle> generatorStableParticles;

    for (const auto& [vtxPos, particles] : primaryVertices) {
      nPrimaryVertices += 1;
      ACTS_VERBOSE("Walking particle tree for primary vertex at "
                   << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z());
      std::size_t nParticles = 0;
      std::size_t nSecondaryVertices = 0;
      std::size_t maxGen = 0;

      auto startSize = unorderedParticlesInitial.size();

      // Find all GENERATOR STABLE particles (i.e. particles that were handed
      // over to the simulation)
      generatorStableParticles.clear();

      ACTS_VERBOSE("Finding generator stable particles in "
                   << particles.size()
                   << " particles recorded for this primary vertex");

      {
        auto s = timerFindStable.sample();
        for (const auto& index : particles) {
          const auto& inParticle = mcParticleCollection.at(index);
          findGeneratorStableParticles(inParticle, generatorStableParticles);
        }
      }

      ACTS_VERBOSE(
          "Have " << generatorStableParticles.size()
                  << " generator stable particles for this primary vertex");

      particlesGeneratorUnordered->reserve(particlesGeneratorUnordered->size() +
                                           generatorStableParticles.size());

      for (const auto& genParticle : generatorStableParticles) {
        nParticles += 1;
        const auto particleId = SimBarcode()
                                    .withParticle(nParticles)
                                    .withVertexPrimary(nPrimaryVertices);
        SimParticle particle =
            EDM4hepUtil::readParticle(genParticle).withParticleId(particleId);
        particlesGeneratorUnordered->push_back(particle);
        ACTS_VERBOSE("+ add GEN particle " << particle);
        ACTS_VERBOSE("  - at " << particle.position().transpose());

        const auto pid = particle.particleId();
        unorderedParticlesInitial.push_back(particle.particleId());
        edm4hepParticleMap[genParticle.getObjectID().index] =
            detail::ParticleInfo{.particleIndex =
                                     unorderedParticlesInitial.size() - 1};
        processChildren(genParticle, pid, unorderedParticlesInitial,
                        parentRelationship, edm4hepParticleMap,
                        nSecondaryVertices, maxGen, getNumHits);
      }

      ACTS_VERBOSE("Primary vertex complete, produced "
                   << (unorderedParticlesInitial.size() - startSize)
                   << " particles and " << nSecondaryVertices
                   << " secondary vertices in " << maxGen << " generations");
      std::span particleSpan{
          std::next(unorderedParticlesInitial.begin(), startSize),
          unorderedParticlesInitial.end()};
      setSubParticleIds(particleSpan);
    }
  }

  ACTS_DEBUG("Found " << nGeneratorParticles << " generator particles and "
                      << unorderedParticlesInitial.size() - nGeneratorParticles
                      << " particles from simulation");
  ACTS_DEBUG("Found " << unorderedParticlesInitial.size()
                      << " particles in total");

  if (unorderedParticlesInitial.empty()) {
    ACTS_WARNING(
        "No particles found in the input file, this will likely cause issues "
        "if sim hits are converted");
  }

  std::optional particlesSimulatedUnordered = std::vector<SimParticle>{};

  auto convertPosition4 = [&](const edm4hep::MCParticle& inParticle) {
    Acts::Vector4 vtxPos4 = {inParticle.getVertex()[0],
                             inParticle.getVertex()[1],
                             inParticle.getVertex()[2], inParticle.getTime()};
    vtxPos4.head<3>() *= Acts::UnitConstants::mm;
    vtxPos4[3] *= Acts::UnitConstants::ns;
    return vtxPos4;
  };

  std::size_t nSkipped = 0;
  {
    Acts::ScopedTimer timer("Converting particles", logger(),
                            Acts::Logging::DEBUG);

    for (const auto& inParticle : mcParticleCollection) {
      ACTS_VERBOSE("Converting particle:\n" << inParticle);

      auto particleIt = edm4hepParticleMap.find(inParticle.getObjectID().index);
      if (particleIt == edm4hepParticleMap.end()) {
        nSkipped += 1;
        continue;
      }

      const auto& [index] = particleIt->second;

      if (!acceptParticle(inParticle) && getNumHits(inParticle) == 0) {
        // Only reject particles if they don't have simhits
        ACTS_VERBOSE(" - skipping particle (no hits AND not accepted)");
        continue;
      }

      const auto& pid = unorderedParticlesInitial.at(index);

      // Copy the particle to the simulated particle container, because we'll
      // make modified version for the "final" state (i.e. after simulation)
      SimParticle particleSimulated = EDM4hepUtil::readParticle(inParticle);
      particleSimulated.setParticleId(pid);
      ACTS_VERBOSE("Have converted particle: " << particleSimulated);

      // Find the decay time of the particle, by looking for the first
      // daughter that marks that it's the endpoint of the parent: this
      // daughter's creation time is the decay time of the parent.
      float time = inParticle.getTime() * Acts::UnitConstants::ns;
      ACTS_VERBOSE("Particle has " << inParticle.getDaughters().size()
                                   << " daughters");
      for (const auto& daughter : inParticle.getDaughters()) {
        if (!daughter.vertexIsNotEndpointOfParent()) {
          Acts::Vector4 pos4 = convertPosition4(daughter);
          time = static_cast<float>(pos4[Acts::eFreeTime]);

          break;
        }
      }

      particleSimulated.finalState().setPosition4(
          inParticle.getEndpoint()[0] * Acts::UnitConstants::mm,
          inParticle.getEndpoint()[1] * Acts::UnitConstants::mm,
          inParticle.getEndpoint()[2] * Acts::UnitConstants::mm, time);

      Acts::Vector3 momentumFinal = {inParticle.getMomentumAtEndpoint()[0],
                                     inParticle.getMomentumAtEndpoint()[1],
                                     inParticle.getMomentumAtEndpoint()[2]};
      particleSimulated.finalState().setDirection(momentumFinal.normalized());
      particleSimulated.finalState().setAbsoluteMomentum(momentumFinal.norm());

      ACTS_VERBOSE(
          "- Updated particle initial -> final, position: "
          << particleSimulated.initialState().fourPosition().transpose()
          << " -> "
          << particleSimulated.finalState().fourPosition().transpose());
      ACTS_VERBOSE(
          "                                     momentum: "
          << particleSimulated.initialState().fourMomentum().transpose()
          << " -> "
          << particleSimulated.finalState().fourMomentum().transpose());

      particlesSimulatedUnordered->push_back(particleSimulated);
    }
  }

  ACTS_DEBUG("Converted " << particlesGeneratorUnordered->size()
                          << " generator particles");
  ACTS_DEBUG("Converted " << particlesSimulatedUnordered->size()
                          << " simulated particles");

  ACTS_DEBUG("Skipped particles: " << nSkipped << " (no hits or not accepted)");

  std::ranges::sort(*particlesGeneratorUnordered, detail::CompareParticleId{});
  std::ranges::sort(*particlesSimulatedUnordered, detail::CompareParticleId{});

  SimParticleContainer particlesGenerator{particlesGeneratorUnordered->begin(),
                                          particlesGeneratorUnordered->end()};

  particlesGeneratorUnordered.reset();

  SimParticleContainer particlesSimulated{particlesSimulatedUnordered->begin(),
                                          particlesSimulatedUnordered->end()};

  particlesSimulatedUnordered.reset();

  std::optional simHitsUnordered = std::vector<SimHit>{};
  ACTS_DEBUG("Reading sim hits from " << m_cfg.inputSimHits.size()
                                      << " sim hit collections");
  {
    Acts::ScopedTimer timer("Reading sim hits", logger(), Acts::Logging::DEBUG);

    const auto& vm = m_cfg.dd4hepDetector->dd4hepDetector().volumeManager();

    auto geometryMapper = [&](std::uint64_t cellId) {
      ACTS_VERBOSE("CellID: " << cellId);

      const auto detElement = vm.lookupDetElement(cellId);

      ACTS_VERBOSE(" -> detElement: " << detElement.name());
      ACTS_VERBOSE("   -> id: " << detElement.id());
      ACTS_VERBOSE("   -> key: " << detElement.key());

      Acts::Vector3 position;
      position
          << detElement.nominal().worldTransformation().GetTranslation()[0],
          detElement.nominal().worldTransformation().GetTranslation()[1],
          detElement.nominal().worldTransformation().GetTranslation()[2];
      position *= Acts::UnitConstants::cm;

      ACTS_VERBOSE("   -> detElement position: " << position.transpose());

      auto it = m_surfaceMap.find(detElement.key());
      if (it == m_surfaceMap.end()) {
        ACTS_ERROR("Unable to find surface for detElement "
                   << detElement.name() << " with cellId " << cellId);
        throw std::runtime_error("Unable to find surface for detElement");
      }
      const auto* surface = it->second;
      if (surface == nullptr) {
        ACTS_ERROR("Unable to find surface for detElement "
                   << detElement.name() << " with cellId " << cellId);
        throw std::runtime_error("Unable to find surface for detElement");
      }
      ACTS_VERBOSE("   -> surface: " << surface->geometryId());
      return surface->geometryId();
    };

    auto particleMapper = [&](const auto& inParticle) {
      auto it = edm4hepParticleMap.find(inParticle.getObjectID().index);
      if (it == edm4hepParticleMap.end()) {
        ACTS_ERROR(
            "SimHit has source particle that we did not see "
            "before, particle=\n"
            << inParticle);
        // If we get this here, this is some sort of bug
        throw std::runtime_error(
            "SimHit has source particle that we did not see before, "
            "but expect to be here");
      }

      auto pid = unorderedParticlesInitial.at(it->second.particleIndex);
      return pid;
    };

    // We will (ab)use the sim hit index to store the association with the
    // incoming edm4hep simhit. The reason is that we will not have the final
    // sim hit index in the collection that's sorted by geometry id, so we can't
    // otherwise build an association map. The index will later be overwritten
    // based
    //
    // This index will be across the input collections,
    // so we need to check if the total count is too large
    std::size_t totalSimHitCount = 0;
    std::int32_t simHitIndex = 0;

    if (m_outputSimHitAssociation.isInitialized()) {
      totalSimHitCount = std::accumulate(
          simHitCollections.begin(), simHitCollections.end(), 0u,
          [&](auto sum, const auto* coll) { return sum + coll->size(); });

      if (totalSimHitCount >
          std::numeric_limits<decltype(simHitIndex)>::max()) {
        ACTS_ERROR(
            "Due to the way the conversion uses a 32bit integer to store the "
            "edm4hep sim hit index, the total number of sim hits across all "
            "input collections must be <= "
            << std::numeric_limits<decltype(simHitIndex)>::max() << ", but is "
            << totalSimHitCount);
        throw std::runtime_error("Total sim hit count is too large");
      }
    }

    for (const auto& [inputHits, name] :
         Acts::zip(simHitCollections, m_cfg.inputSimHits)) {
      ACTS_VERBOSE("SimHit collection " << name << " has " << inputHits->size()
                                        << " hits");

      simHitsUnordered->reserve(simHitsUnordered->size() + inputHits->size());

      for (const auto& hit : *inputHits) {
        auto simHit = EDM4hepUtil::readSimHit(hit, particleMapper,
                                              geometryMapper, simHitIndex);

        ACTS_VERBOSE("Converted sim hit for truth particle: "
                     << simHit.particleId() << " at "
                     << simHit.fourPosition().transpose() << " with time "
                     << simHit.time());

        // Increase hit count in generated and simulated particles
        if (auto itSim = particlesSimulated.find(simHit.particleId());
            itSim != particlesSimulated.end()) {
          ACTS_VERBOSE("Found associated simulated particle");
          itSim->finalState().setNumberOfHits(
              itSim->finalState().numberOfHits() + 1);
        } else if (auto itGen = particlesGenerator.find(simHit.particleId());
                   itGen != particlesGenerator.end()) {
          ACTS_VERBOSE("Found associated generator particle");
          itGen->finalState().setNumberOfHits(
              itGen->finalState().numberOfHits() + 1);
        } else {
          const auto& ptcl = ActsPlugins::EDM4hepUtil::getParticle(hit);
          ACTS_ERROR("SimHit (r="
                     << Acts::VectorHelpers::perp(simHit.position())
                     << ", z=" << simHit.position()[Acts::eFreePos2]
                     << ") has source particle that we did not see before:\n"
                     << ptcl);
          double particleR =
              std::hypot(ptcl.getVertex()[0], ptcl.getVertex()[1]) *
              Acts::UnitConstants::mm;
          double particleZ = ptcl.getVertex()[2] * Acts::UnitConstants::mm;

          double particleREnd =
              std::hypot(ptcl.getEndpoint()[0], ptcl.getEndpoint()[1]) *
              Acts::UnitConstants::mm;
          double particleZEnd = ptcl.getEndpoint()[2] * Acts::UnitConstants::mm;

          ACTS_ERROR("Particle loc: " << particleR << ", " << particleZ
                                      << " -> " << particleREnd << ", "
                                      << particleZEnd);
          continue;
        }

        simHitsUnordered->push_back(std::move(simHit));
        simHitIndex += 1;
      }
    }
  }

  ACTS_DEBUG("Read " << simHitsUnordered->size() << " sim hits in total");

  std::ranges::sort(*simHitsUnordered, detail::CompareGeometryId{});

  SimHitContainer simHits{simHitsUnordered->begin(), simHitsUnordered->end()};
  simHitsUnordered.reset();

  // We now know the final indices of the indices in the output simhit
  // collection. In the next step, the indices along the particle path will be
  // rewritten, so we can build an assotiaion map here.

  if (m_outputSimHitAssociation.isInitialized()) {
    ActsPlugins::EDM4hepUtil::SimHitAssociation simHitAssociation;
    simHitAssociation.reserve(simHits.size());

    // @TODO: Make this optional depending on the output key setting
    Acts::ScopedTimer timer("Building sim hit association", logger(),
                            Acts::Logging::DEBUG);
    for (const auto&& [indexInColl, hit] : Acts::enumerate(simHits)) {
      std::size_t index = hit.index();
      // find hit for this index in the input collections
      for (const auto&& [name, coll] :
           Acts::zip(m_cfg.inputSimHits, simHitCollections)) {
        if (index >= coll->size()) {
          index -= coll->size();
          continue;
        }

        ACTS_VERBOSE("Hit assoc int -> ext #" << indexInColl << " -> "
                                              << coll->at(index).id() << " "
                                              << hit.position().transpose());

        simHitAssociation.add(indexInColl, coll->at(index));

        break;
      }
    }

    if (simHitAssociation.size() != simHits.size()) {
      ACTS_ERROR("Sim hit association size " << simHitAssociation.size()
                                             << " does not match sim hit size "
                                             << simHits.size());
      throw std::runtime_error("Sim hit association size mismatch");
    }

    m_outputSimHitAssociation(ctx, std::move(simHitAssociation));
  }

  if (m_cfg.sortSimHitsInTime) {
    Acts::ScopedTimer timer("Sorting sim hits in time", logger(),
                            Acts::Logging::DEBUG);
    std::multimap<ActsFatras::Barcode, std::size_t> hitsByParticle;

    for (std::size_t i = 0; i < simHits.size(); ++i) {
      hitsByParticle.insert({simHits.nth(i)->particleId(), i});
    }

    for (auto it = hitsByParticle.begin(), end = hitsByParticle.end();
         it != end; it = hitsByParticle.upper_bound(it->first)) {
      ACTS_VERBOSE("Particle " << it->first << " has "
                               << hitsByParticle.count(it->first) << " hits");

      std::vector<std::size_t> hitIndices;
      hitIndices.reserve(hitsByParticle.count(it->first));
      for (auto hitIndex : makeRange(hitsByParticle.equal_range(it->first))) {
        hitIndices.push_back(hitIndex.second);
      }

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        ACTS_VERBOSE("Before sorting:");
        for (const auto& hitIdx : hitIndices) {
          ACTS_VERBOSE(" - " << hitIdx << " / " << simHits.nth(hitIdx)->index()
                             << " t=" << simHits.nth(hitIdx)->time());
        }
      }

      std::ranges::sort(hitIndices, {}, [&simHits](std::size_t h) {
        return simHits.nth(h)->time();
      });

      for (std::size_t i = 0; i < hitIndices.size(); ++i) {
        auto& hit = *simHits.nth(hitIndices[i]);
        // SimHit does not have setters, so need to create a new one
        hit = {hit.geometryId(),     hit.particleId(),
               hit.fourPosition(),   hit.momentum4Before(),
               hit.momentum4After(), static_cast<std::int32_t>(i)};
      }

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        ACTS_VERBOSE("After sorting:");
        for (const auto& hitIdx : hitIndices) {
          ACTS_VERBOSE(" - " << hitIdx << " / " << simHits.nth(hitIdx)->index()
                             << " t=" << simHits.nth(hitIdx)->time());
        }
      }
    }
  } else {
    // Reset all indices to -1 to indicate we don't know the ordering along the
    // particle
    for (auto& hit : simHits) {
      // SimHit does not have setters, so need to create a new one
      hit = {hit.geometryId(),      hit.particleId(),     hit.fourPosition(),
             hit.momentum4Before(), hit.momentum4After(), -1};
    }
  }

  std::vector<SimVertex> simVerticesUnordered;

  auto maybeAddVertex = [&](const Acts::Vector4& vtxPos4,

                            SimVertexBarcode vtxId) -> SimVertex& {
    auto getMinDistance = [&]() {
      std::stringstream sstr;
      auto closestIt = std::ranges::min_element(
          simVerticesUnordered, {}, [&vtxPos4](const auto& v) {
            return (v.position4.template head<3>() - vtxPos4.template head<3>())
                .norm();
          });

      if (closestIt != simVerticesUnordered.end()) {
        sstr << (closestIt->position4.head<3>() - vtxPos4.head<3>()).norm();
      } else {
        sstr << "[NONE]";
      }

      return sstr.str();
    };

    auto vertexIt =
        std::ranges::find_if(simVerticesUnordered, [&](const auto& v) {
          return (v.position4.template head<3>() - vtxPos4.template head<3>())
                         .norm() < Acts::UnitConstants::mm * 1e-3 &&
                 v.id == vtxId;
        });

    SimVertex* vertex = nullptr;
    // We don't have a vertex for this position + id yet
    if (vertexIt == simVerticesUnordered.end()) {
      ACTS_VERBOSE("Adding new vertex: position="
                   << vtxPos4.template head<3>().transpose() << " id=" << vtxId
                   << " (closest existing=" << getMinDistance() << ")");
      vertex = &simVerticesUnordered.emplace_back(vtxId, vtxPos4);
    } else {
      vertex = &*vertexIt;
      ACTS_VERBOSE("Reusing existing vertex: position="
                   << vtxPos4.template head<3>().transpose()
                   << " id=" << vertex->id);
    }

    assert(vertex != nullptr);

    return *vertex;
  };

  {
    Acts::ScopedTimer timer("Finding source vertices for particles", logger(),
                            Acts::Logging::DEBUG);

    std::size_t nParticlesWithHits = 0;
    for (const auto& particle : particlesSimulated) {
      // Add current particle to the outgoing particles of the vertex

      if (particle.finalState().numberOfHits() == 0) {
        // Only produce vertices for particles that actually produced any hits
        continue;
      }
      nParticlesWithHits += 1;

      SimVertex& vertex = maybeAddVertex(
          particle.fourPosition(), SimVertexBarcode{particle.particleId()});
      vertex.outgoing.insert(particle.particleId());
    }

    ACTS_DEBUG("Made " << simVerticesUnordered.size() << " vertices from "
                       << nParticlesWithHits
                       << " simulated particles with hits");
  }

  std::ranges::sort(simVerticesUnordered, detail::CompareVertexId{});

  ACTS_DEBUG("Converted number of vertices: " << simVerticesUnordered.size());
  SimVertexContainer simVertices{simVerticesUnordered.begin(),
                                 simVerticesUnordered.end()};

  m_outputParticlesGenerator(ctx, std::move(particlesGenerator));
  m_outputParticlesSimulation(ctx, std::move(particlesSimulated));
  m_outputSimVertices(ctx, std::move(simVertices));

  for (const auto&& [i, hit] : Acts::enumerate(simHits)) {
    ACTS_VERBOSE("- " << i << " " << hit.fourPosition().transpose());
  }

  m_outputSimHits(ctx, std::move(simHits));

  return ProcessCode::SUCCESS;
}

void EDM4hepSimInputConverter::processChildren(
    const edm4hep::MCParticle& inParticle, SimBarcode parentId,
    std::vector<SimBarcode>& particles, ParentRelationship& parentRelationship,
    std::unordered_map<int, detail::ParticleInfo>& particleMap,
    std::size_t& nSecondaryVertices, std::size_t& maxGen,
    const std::function<std::uint8_t(const edm4hep::MCParticle&)>& getNumHits)
    const {
  constexpr auto indent = [&](std::size_t n) {
    std::string result;
    for (std::size_t i = 0; i < n; ++i) {
      result += "  ";
    }
    return result;
  };

  const std::size_t gen = parentId.generation();
  maxGen = std::max(maxGen, gen);

  ACTS_VERBOSE(indent(gen) << "  - processing daughters for input particle "
                           << inParticle.id());
  ACTS_VERBOSE(indent(gen) << "    -> found " << inParticle.daughters_size()
                           << " daughter(s)");

  bool parentDecayed =
      std::any_of(inParticle.daughters_begin(), inParticle.daughters_end(),
                  [](const edm4hep::MCParticle& daughter) {
                    return !daughter.vertexIsNotEndpointOfParent();
                  });
  std::size_t secondaryVertex = 0;

  Acts::Vector3 start{inParticle.getVertex().x, inParticle.getVertex().y,
                      inParticle.getVertex().z};
  Acts::Vector3 end{inParticle.getEndpoint().x, inParticle.getEndpoint().y,
                    inParticle.getEndpoint().z};
  double distance = (end - start).norm() * Acts::UnitConstants::mm;

  constexpr double tolerance = Acts::s_onSurfaceTolerance;
  if (parentDecayed && distance > tolerance) {
    ACTS_VERBOSE(indent(gen) << "    -> parent decays");
    secondaryVertex = ++nSecondaryVertices;
  }

  std::size_t parentIndex = particles.size() - 1;

  std::size_t nParticles = 0;
  for (const auto& daughter : inParticle.getDaughters()) {
    if (auto pIt = particleMap.find(daughter.getObjectID().index);
        pIt != particleMap.end()) {
      ACTS_WARNING("Skipping particle #"
                   << daughter.getObjectID()
                   << " that we've already processed, this should not happen "
                      "at this stage.");
      continue;
    }

    // Check if we want to keep this particle (includes checking if any of the
    // descendants have hits that we'll convert)
    if (!acceptParticle(daughter) &&
        !particleOrDescendantsHaveHits(daughter, getNumHits)) {
      ACTS_VERBOSE(indent(gen + 1) << "  - skipping particle " << daughter.id()
                                   << " (no hits AND not accepted)");
      // Only reject particles if they and their descendants don't have simhits
      continue;
    }

    SimParticle particle = EDM4hepUtil::readParticle(daughter);

    auto pid = parentId.makeDescendant(nParticles);
    nParticles += 1;
    if (daughter.vertexIsNotEndpointOfParent()) {
      // incoming particle survived, interaction via descendant
    } else {
      // incoming particle decayed
      pid = pid.withVertexSecondary(secondaryVertex);
    }
    particle.setParticleId(pid);

    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "+ add particle " << particle << " ("
                 << particle.particleId() << ") from #"
                 << daughter.getObjectID());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - generation: " << particle.particleId().generation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - at " << particle.fourPosition().transpose());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - createdInSim: "
                 << daughter.isCreatedInSimulation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - vertexIsNotEndpointOfParent: "
                 << daughter.vertexIsNotEndpointOfParent());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - isStopped: " << daughter.isStopped());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - endpoint: " << daughter.getEndpoint().x << ", "
                 << daughter.getEndpoint().y << ", "
                 << daughter.getEndpoint().z);

    particles.push_back(particle.particleId());
    particleMap[daughter.getObjectID().index] =
        detail::ParticleInfo{.particleIndex = particles.size() - 1};
    parentRelationship[particles.size() - 1] = parentIndex;

    processChildren(daughter, pid, particles, parentRelationship, particleMap,
                    nSecondaryVertices, maxGen, getNumHits);
  }
}

void EDM4hepSimInputConverter::setSubParticleIds(
    std::span<SimBarcode> particles) {
  std::vector<std::size_t> numByGeneration;
  numByGeneration.reserve(10);

  for (auto& particle : particles) {
    if (particle.generation() >= numByGeneration.size()) {
      numByGeneration.resize(particle.generation() + 1, 0);
    }
    unsigned long nextSubParticle = numByGeneration[particle.generation()];
    numByGeneration[particle.generation()] += 1;

    particle = particle.withSubParticle(nextSubParticle);
  }
}

}  // namespace ActsExamples
