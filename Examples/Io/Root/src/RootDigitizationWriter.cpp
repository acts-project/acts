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
      m_cfg(cfg) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  // Setup ROOT File
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
  }

  m_outputFile->cd();

  // Analyze the smearers
  std::vector<
      std::pair<Acts::GeometryIdentifier, std::unique_ptr<DigitizationTree>>>
      dTrees;
  if (not m_cfg.smearers.empty()) {
    ACTS_DEBUG("Smearers are present, preparing trees.");
    for (size_t ikv = 0; ikv < m_cfg.smearers.size(); ++ikv) {
      auto geoID = m_cfg.smearers.idAt(ikv);
      auto value = m_cfg.smearers.valueAt(ikv);
      std::visit(
          [&](auto&& smearer) {
            auto dTree = std::make_unique<DigitizationTree>(geoID);
            typename decltype(smearer.first)::ParSet::ParametersVector pv;
            typename decltype(smearer.first)::ParSet pset(std::nullopt, pv);
            dTree->setupBoundRecBranches(pset);
            dTrees.push_back({geoID, std::move(dTree)});
          },
          value);
    }
  }
  m_outputTrees = Acts::GeometryHierarchyMap<std::unique_ptr<DigitizationTree>>(
      std::move(dTrees));
}

ActsExamples::RootDigitizationWriter::~RootDigitizationWriter() {
  /// Close the file if it's yours
  m_outputFile->cd();
  for (auto dTree = m_outputTrees.begin(); dTree != m_outputTrees.end();
       ++dTree) {
    (*dTree)->tree->Write();
  }
  m_outputFile->Close();
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootDigitizationWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::GeometryIdMultimap<
        Acts::FittableMeasurement<ActsExamples::DigitizedHit>>& measurements) {
  const auto& simHitContainer =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  for (auto& [key, value] : measurements) {
    std::visit(
        [&](auto&& m) {
          const auto& surface = m.referenceObject();
          Acts::GeometryIdentifier sIdentifier = surface.geometryId();
          auto dTreeItr = m_outputTrees.find(sIdentifier);
          if (dTreeItr != m_outputTrees.end()) {
            auto& dTree = *dTreeItr;
            // Fill the identification
            dTree->fillIdentification(ctx.eventNumber, sIdentifier);
            auto hitIndices = m.sourceLink().hitIndices();
            std::vector<SimHit> simHits;
            simHits.reserve(hitIndices.size());
            for (auto& hi : hitIndices) {
              auto nthhit = simHitContainer.nth(hi);
              simHits.push_back(*nthhit);
            }

            auto tParams = truthParameters(ctx.geoContext, surface, simHits);
            dTree->fillTruthParameters(std::get<0>(tParams),
                                       std::get<1>(tParams),
                                       std::get<2>(tParams));
            dTree->fillBoundMeasurement(m);
            dTree->tree->Fill();
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
    for (const auto& hit : simulatedHits) {
      avePos += hit.position4();
      direction += hit.unitDirection();
    }
    double denom = 1. / simulatedHits.size();
    avePos *= denom;
    direction *= denom;
    direction = direction.normalized();
    auto sIntersection =
        surface.intersect(gCtx, avePos.segment<3>(0), direction, false);

    position = sIntersection.intersection.position;
    ctime = avePos[Acts::eTime];
  }
  auto lpResult = surface.globalToLocal(gCtx, position, direction);
  Acts::Vector2D lposition(0., 0.);
  if (not lpResult.ok()) {
    ACTS_WARNING(
        "GlobalToLocal did not succeed, will fill (0.,0.) for local position.");
  } else {
    lposition = lpResult.value();
  }

  return {lposition, Acts::VectorHelpers::makeVector4(position, ctime),
          direction};
}
