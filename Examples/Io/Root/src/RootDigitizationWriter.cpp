// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootDigitizationWriter.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimIdentifier.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/Intersection.hpp>
#include <Acts/Utilities/Units.hpp>

#include <ios>
#include <optional>
#include <stdexcept>

#include <TFile.h>
#include <TString.h>

ActsExamples::RootDigitizationWriter::RootDigitizationWriter(
    const ActsExamples::RootDigitizationWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputMeasurements, "RootDigitizationWriter", lvl),
      m_cfg(cfg),
      m_outputFile(cfg.rootFile) {
  // inputClusters is already checked by base constructor
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();

  // Analyze the smearers
  if (not m_cfg.smearers.empty()) {
    ACTS_DEBUG("Smearers are present, preparing trees.");
    for (auto& [key, value] : m_cfg.smearers) {
      auto& geoID = key;
      std::visit(
          [&](auto&& smearer) {
            auto dTree = std::make_unique<DigitizationTree>(geoID);
            typename decltype(smearer.first)::ParSet::ParametersVector pv;
            typename decltype(smearer.first)::ParSet pset(std::nullopt, pv);
            dTree->setupBoundRecBranches(pset);
            m_outputTrees.emplace_hint(m_outputTrees.end(), geoID,
                                       std::move(dTree));
          },
          value);
    }
  }
}

ActsExamples::RootDigitizationWriter::~RootDigitizationWriter() {
  /// Close the file if it's yours
  if (m_cfg.rootFile == nullptr) {
    m_outputFile->cd();
    for (auto& [key, value] : m_outputTrees) {
      value->tree->Write();
    }
    m_outputFile->Close();
  }
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::GeometryIdMultimap<
        Acts::FittableMeasurement<ActsExamples::DigitizedHit>>& measurements) {
  // TODO move to retrieved simulated hits container rather than object
  // composition
  // const auto& simHitContainer =
  //    ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  for (auto& [key, value] : measurements) {
    std::visit(
        [&](auto&& m) {
          const auto& surface = m.referenceObject();
          Acts::GeometryIdentifier sIdentifier = surface.geometryId();
          // TODO fix the dedicated search on volume to hiearachical one
          auto dTreeItr = m_outputTrees.find(
              Acts::GeometryIdentifier().setVolume(sIdentifier.volume()));
          if (dTreeItr != m_outputTrees.end()) {
            auto& dTree = *(dTreeItr->second).get();
            // Fill the identification
            dTree.fillIdentification(ctx.eventNumber, sIdentifier);
            // TODO change the contributed sim hits collection to indices
            auto simHits = m.sourceLink().simulatedHits();
            auto tParams = truthParameters(ctx.geoContext, surface, simHits);
            dTree.fillTruthParameters(std::get<0>(tParams),
                                      std::get<1>(tParams),
                                      std::get<2>(tParams));
            dTree.fillBoundMeasurement(m);
            dTree.tree->Fill();
          }
        },
        value);
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  return ActsExamples::ProcessCode::SUCCESS;
}

std::tuple<Acts::Vector2D, Acts::Vector4D, Acts::Vector3D>
ActsExamples::RootDigitizationWriter::truthParameters(
    const Acts::GeometryContext& gCtx, const Acts::Surface& surface,
    const std::vector<ActsExamples::SimHit>& simulatedHits) {
  Acts::Vector3D position(0., 0., 0.);
  Acts::Vector3D direction(0., 0., 0.);
  double ctime = 0.;

  if (simulatedHits.size() == 1 and
      surface.isOnSurface(gCtx, simulatedHits[0].position(),
                          simulatedHits[0].unitDirection())) {
    position = simulatedHits[0].position();
    direction = simulatedHits[0].unitDirection();
    ctime = simulatedHits[0].time();

  } else if (simulatedHits.size() > 1) {
    Acts::Vector4D avePos(0., 0., 0., 0.);
    if (simulatedHits.size() > 1) {
      for (const auto& hit : simulatedHits) {
        avePos += hit.position4();
        direction += hit.unitDirection();
      }
      double denom = 1. / simulatedHits.size();
      avePos *= denom;
      direction *= denom;
      direction = direction.normalized();
    }
    auto sIntersection =
        surface.intersect(gCtx, avePos.segment<3>(0), direction, false);

    position = sIntersection.intersection.position;
    ctime = avePos[Acts::eTime];
  }
  const auto& lposition =
      surface.globalToLocal(gCtx, position, direction).value();

  return {lposition, Acts::VectorHelpers::makeVector4(position, ctime),
          direction};
}
