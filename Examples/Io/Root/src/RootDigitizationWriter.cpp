// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootDigitizationWriter.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"

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
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputHitSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
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
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  const auto& hitSimHitsMap =
      ctx.eventStore.get<IndexMultimap<Index>>(m_cfg.inputHitSimHitsMap);

  std::vector<Index> simHitIndices;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& meas = measurements[hitIdx];

    std::visit(
        [&](const auto& m) {
          const auto& surface = m.referenceObject();
          Acts::GeometryIdentifier geoId = surface.geometryId();

          // find the corresponding output tree
          auto dTreeItr = m_outputTrees.find(geoId);
          if (dTreeItr == m_outputTrees.end()) {
            return;
          }
          auto& dTree = *dTreeItr;

          // Fill the identification
          dTree->fillIdentification(ctx.eventNumber, geoId);

          // Find the contributing simulated hits
          auto simHitsRange = makeRange(hitSimHitsMap.equal_range(hitIdx));
          simHitIndices.clear();
          simHitIndices.reserve(simHitsRange.size());
          for (auto [_, simHitIdx] : simHitsRange) {
            simHitIndices.push_back(simHitIdx);
          }

          auto [local, pos4, dir] = extractAverageTruthParameters(
              ctx.geoContext, surface, simHits, simHitIndices);
          dTree->fillTruthParameters(local, pos4, dir);
          dTree->fillBoundMeasurement(m);
          dTree->tree->Fill();
        },
        meas);
  }

  return ActsExamples::ProcessCode::SUCCESS;
}

std::tuple<Acts::Vector2D, Acts::Vector4D, Acts::Vector3D>
ActsExamples::RootDigitizationWriter::extractAverageTruthParameters(
    const Acts::GeometryContext& gCtx, const Acts::Surface& surface,
    const SimHitContainer& simHits, const std::vector<Index>& simHitIndices) {
  Acts::Vector2D avgLocal = Acts::Vector2D::Zero();
  Acts::Vector4D avgPos4 = Acts::Vector4D::Zero();
  Acts::Vector3D avgDir = Acts::Vector3D::Zero();

  for (Index simHitIdx : simHitIndices) {
    // we assume that the indices are within valid ranges so we do not need to
    // check their validity again.
    const auto& simHit = *simHits.nth(simHitIdx);

    // transforming first to local positions and average that ensures that the
    // averaged position is still on the surface. the averaged global position
    // might not be on the surface anymore.
    auto result =
        surface.globalToLocal(gCtx, simHit.position(), simHit.unitDirection());
    if (result.ok()) {
      avgLocal += result.value();
    } else {
      ACTS_WARNING("Simulated hit "
                   << simHitIdx << " is not on the corresponding surface "
                   << surface.geometryId() << "; use [0,0] as local position");
    }
    // global position should already be at the intersection. no need to perform
    // an additional intersection call.
    avgPos4 += simHit.position4();
    avgDir += simHit.unitDirection();
  }

  // only need to average if there are at least two inputs
  if (2 <= simHitIndices.size()) {
    double scale = 1 / simHitIndices.size();
    avgLocal *= scale;
    avgPos4 *= scale;
    avgDir.normalize();
  }

  return {avgLocal, avgPos4, avgDir};
}
