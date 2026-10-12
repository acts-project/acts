// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Podio/PodioReader.hpp"
#include "ActsExamples/Io/Podio/PodioWriter.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "ReconstructionFixture.hpp"
#ifdef ACTS_CALORIMETER_EDM4HEP_JETS
#include "Jets/JetFixture.hpp"
#endif

namespace {
using namespace ActsExamples;
namespace Fixture = CalorimeterFixture;
void require(bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}
void checkOutput(const std::filesystem::path& path) {
  ActsPlugins::PodioUtil::ROOTReader reader;
  reader.openFile(path.string());
  require(reader.getEntries("events") == 3,
          "Unexpected calorimeter event count");
  for (unsigned int event = 0; event < 3; ++event) {
    podio::Frame frame = reader.readEntry("events", event);
    const auto& hits = frame.get<edm4hep::CalorimeterHitCollection>("CaloHits");
    const auto& links =
        frame.get<edm4hep::CaloHitSimCaloHitLinkCollection>("CaloLinks");
    const auto& simHits =
        frame.get<edm4hep::SimCalorimeterHitCollection>("SimCaloHits");
    const auto& clusters =
        frame.get<edm4hep::ClusterCollection>("CaloClusters");
    const auto& times =
        frame.get<podio::UserDataCollection<double>>("CaloClusterTimes");
    const auto& seeds = frame.get<podio::UserDataCollection<std::uint64_t>>(
        "CaloClusterSeedCellIds");
    require(hits.size() == 4 && links.size() == 5 && simHits.size() == 5 &&
                clusters.size() == 3 && times.size() == 3 && seeds.size() == 3,
            "Unexpected persistent collection size");
    require(hits[0].getCellID() == 7 && hits[1].getCellID() == 42 &&
                hits[2].getCellID() == Fixture::highCellId &&
                hits[3].getCellID() == Fixture::highCellId + 1 &&
                std::abs(hits[2].getEnergy() - 0.08) < 1e-7 &&
                std::abs(hits[2].getTime() - 4) < 1e-6 &&
                hits[2].getPosition().x == 1500,
            "Unexpected persistent hit response");
    double weightSum = 0;
    for (const auto& link : links) {
      require(link.getFrom().isAvailable() && link.getTo().isAvailable() &&
                  link.getFrom().getCellID() == link.getTo().getCellID(),
              "Broken truth relation");
      if (link.getFrom().getCellID() == Fixture::highCellId) {
        weightSum += link.getWeight();
        const double expected =
            link.getTo().getObjectID().index == 0 ? 0.75 : 0.25;
        require(std::abs(link.getWeight() - expected) < 1e-6,
                "Unexpected truth weight");
      }
      for (const auto& contribution : link.getTo().getContributions()) {
        require(contribution.isAvailable() &&
                    contribution.getParticle().isAvailable() &&
                    contribution.getParticle().getPDG() == 211,
                "Broken contribution/particle relation");
      }
    }
    require(std::abs(weightSum - 1) < 1e-6, "Truth weights do not sum to one");
    require(std::abs(clusters[0].getEnergy() - 0.06) < 1e-7 &&
                std::abs(clusters[1].getEnergy() - 0.04) < 1e-7 &&
                std::abs(clusters[2].getEnergy() - 0.12) < 1e-7 &&
                clusters[0].getHits(0) == hits[0] &&
                clusters[1].getHits(0) == hits[1] &&
                clusters[2].hits_size() == 2 &&
                clusters[2].getHits(0) == hits[2] &&
                clusters[2].getHits(1) == hits[3],
            "Broken cluster energy/hit relations");
    require(
        seeds[0] == 7 && seeds[1] == 42 && seeds[2] == Fixture::highCellId &&
            std::abs(times[0] - 6) < 1e-6 && std::abs(times[1] - 2) < 1e-6 &&
            std::abs(times[2] - 10. / 3) < 1e-6 &&
            std::abs(clusters[2].getPosition().y - 40. / 3) < 1e-5,
        "Broken cluster sidecars or centre");
    require(
        frame.getParameter<std::string>("generator.tag") == "synthetic" &&
            frame.getParameter<double>("acts.calo.schemaVersion") == 1 &&
            frame.getParameter<double>("acts.calo.response.energyScale") == 2 &&
            frame.getParameter<std::vector<double>>(
                "acts.calo.response.timeWindow") ==
                std::vector<double>({0, 10}) &&
            frame.getParameter<std::string>("acts.calo.output.clusterTimes") ==
                "CaloClusterTimes" &&
            frame.getParameter<std::vector<std::string>>(
                "acts.calo.clustering.neighbours") ==
                std::vector<std::string>{
                    std::to_string(Fixture::highCellId) + ":" +
                    std::to_string(Fixture::highCellId + 1)},
        "Missing or incorrect reconstruction metadata");
#ifdef ACTS_CALORIMETER_EDM4HEP_JETS
    const auto& jets =
        frame.get<edm4hep::ReconstructedParticleCollection>("CaloJets");
    require(jets.size() == 2 && jets[0].clusters_size() == 2 &&
                jets[1].clusters_size() == 1 &&
                jets[0].getClusters(0) == clusters[0] &&
                jets[0].getClusters(1) == clusters[2] &&
                jets[1].getClusters(0) == clusters[1] &&
                std::abs(jets[0].getEnergy() - 0.18) < 1e-7 &&
                std::abs(jets[1].getEnergy() - 0.04) < 1e-7 &&
                jets[0].getMass() > 0,
            "Broken persistent jet relations, energy or mass");
    const Acts::Vector3 expected =
        0.06 * Acts::Vector3(1500, 100, 5).normalized() +
        0.12 * Acts::Vector3(1500, 40. / 3, 5).normalized();
    const auto momentum = jets[0].getMomentum();
    require(
        (Acts::Vector3(momentum.x, momentum.y, momentum.z) - expected).norm() <
                1e-7 &&
            frame.getParameter<double>("acts.calo.jets.radius") == 0.4 &&
            frame.getParameter<std::vector<double>>("acts.calo.jets.origin") ==
                std::vector<double>({0, 0, 0}),
        "Incorrect persistent jet momentum or configuration");
    // Resolve the complete jet -> cluster -> hit -> sim hit -> contribution ->
    // MC chain.
    for (const auto& jet : jets) {
      for (const auto& cluster : jet.getClusters()) {
        for (const auto& hit : cluster.getHits()) {
          bool matched = false;
          for (const auto& link : links) {
            if (link.getFrom() == hit) {
              matched = true;
            }
          }
          require(matched, "Jet constituent lacks a truth association");
        }
      }
    }
#endif
  }
}
}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " OUTPUT_DIRECTORY\n";
    return 1;
  }
  const std::filesystem::path directory(argv[1]);
  std::filesystem::create_directories(directory);
  const auto inputPath = directory / "calorimeter-input.root";
  const auto outputPath = directory / "calorimeter-output.root";
  {
    ActsPlugins::PodioUtil::ROOTWriter writer(inputPath.string());
    for (int event = 0; event < 3; ++event) {
      auto frame = Fixture::makeFrame(true);
      frame.putParameter("generator.tag", std::string("synthetic"));
      writer.writeFrame(frame, "events");
    }
    writer.finish();
  }
  Sequencer::Config sequence;
  sequence.events = 3;
  sequence.numThreads = 1;
  sequence.outputDir = directory.string();
  Sequencer sequencer(sequence);
  PodioReader::Config reader;
  reader.inputPath = inputPath;
  sequencer.addReader(
      std::make_shared<PodioReader>(reader, Acts::Logging::INFO));
  auto metadata = Fixture::metadataConfig();
  sequencer.addAlgorithm(
      std::make_shared<EDM4hepCalorimeterInputConverter>(metadata.input));
  sequencer.addAlgorithm(
      std::make_shared<CalorimeterDigitizationAlgorithm>(metadata.response));
  sequencer.addAlgorithm(
      std::make_shared<CalorimeterClusteringAlgorithm>(metadata.clustering));
  auto hitConverter =
      std::make_shared<EDM4hepCalorimeterOutputConverter>(metadata.hitOutput);
  auto clusterConverter =
      std::make_shared<EDM4hepCalorimeterClusterOutputConverter>(
          metadata.clusterOutput);
  sequencer.addAlgorithm(hitConverter);
  sequencer.addAlgorithm(clusterConverter);
  auto collections = hitConverter->collections();
  const auto clusterCollections = clusterConverter->collections();
  collections.insert(collections.end(), clusterCollections.begin(),
                     clusterCollections.end());
#ifdef ACTS_CALORIMETER_EDM4HEP_JETS
  const auto jets = Fixture::jetConfig();
  sequencer.addAlgorithm(std::make_shared<CalorimeterJetAlgorithm>(jets));
  auto jetConverter = std::make_shared<EDM4hepCalorimeterJetOutputConverter>(
      Fixture::jetOutputConfig());
  sequencer.addAlgorithm(jetConverter);
  metadata.additional = jetConverter->metadata();
  collections.push_back(jetConverter->collections().front());
#endif
  sequencer.addAlgorithm(
      std::make_shared<EDM4hepCalorimeterMetadata>(metadata));
  PodioWriter::Config writer;
  writer.inputFrame = metadata.outputFrame;
  writer.outputPath = outputPath.string();
  writer.category = "events";
  writer.collections = collections;
  sequencer.addWriter(
      std::make_shared<PodioWriter>(writer, Acts::Logging::INFO));
  const auto result = sequencer.run();
  if (result != 0) {
    return result;
  }
  checkOutput(outputPath);
  std::cout
      << "Calorimeter reconstruction ROOT round trip passed for three events\n";
  return 0;
}
