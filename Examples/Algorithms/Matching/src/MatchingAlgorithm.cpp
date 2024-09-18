// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Matching/MatchingAlgorithm.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <cstddef>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include "TRandom.h"

ActsExamples::MatchingAlgorithm::MatchingAlgorithm(
    ActsExamples::MatchingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("MatchingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputTrackParametersMS.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackParametersMS) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackParametersMS.emplace_back(
        std::make_unique<ReadDataHandle<TrackParametersContainer>>(
            this, "inputTrackParametersMS#" +
                      std::to_string(m_inputTrackParametersMS.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackContainerMS.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackContainerMS) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackContainerMS.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTrackContainerMS#" +
                      std::to_string(m_inputTrackContainerMS.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackParametersVT.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackParametersVT) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackParametersVT.emplace_back(
        std::make_unique<ReadDataHandle<TrackParametersContainer>>(
            this, "inputTrackParametersVT#" +
                      std::to_string(m_inputTrackParametersVT.size())));
    handle->initialize(spName);
  }

  if (m_cfg.inputTrackContainerVT.empty()) {
    throw std::invalid_argument("Missing track parameter collection");
  }
  for (const auto& spName : m_cfg.inputTrackContainerVT) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputTrackContainerVT.emplace_back(
        std::make_unique<ReadDataHandle<ConstTrackContainer>>(
            this, "inputTrackContainerVT#" +
                      std::to_string(m_inputTrackContainerVT.size())));
    handle->initialize(spName);
  }

  // m_inputTrackParametersMS.initialize(m_cfg.inputTrackParametersMS);
  // m_inputTrackParametersVT.initialize(m_cfg.inputTrackParametersVT);
  // m_inputTrackContainerMS.initialize(m_cfg.inputTrackContainerMS);
  // m_inputTrackContainerVT.initialize(m_cfg.inputTrackContainerVT);

  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
  m_outputMatchedTracks.initialize(m_cfg.outputMatchedTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMapVT.initialize(
      m_cfg.inputMeasurementParticlesMapVT);
  m_inputMeasurementParticlesMapMS.initialize(
      m_cfg.inputMeasurementParticlesMapMS);
}

ActsExamples::ProcessCode ActsExamples::MatchingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // const auto& paramsVTinput = m_inputTrackParametersVT(ctx);
  // const auto& paramsMSinput = m_inputTrackParametersMS(ctx);
  // const auto& trackContainterVT = m_inputTrackContainerVT(ctx);
  // const auto& trackContainterMS = m_inputTrackContainerMS(ctx);
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMapMS = m_inputMeasurementParticlesMapMS(ctx);
  const auto& hitParticlesMapVT = m_inputMeasurementParticlesMapVT(ctx);

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{m_cfg.px, m_cfg.py, m_cfg.pz});

  Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator> extrapolator(
      Acts::EigenStepper<>(m_cfg.magneticField),
      Acts::Navigator({m_cfg.trackingGeometry},
                      logger().cloneWithSuffix("Navigator")),
      logger().cloneWithSuffix("Propagator"));

  // VOID PROPAGATOR
  // Set up EigenStepper
  Acts::EigenStepper<> stepper(m_cfg.magneticField);

  // Set up propagator with void navigator
  // auto propagator = std::make_shared<Acts::Propagator<Acts::EigenStepper<>,
  // Acts::VoidNavigator>>(
  //    stepper, Acts::VoidNavigator{}, logger().cloneWithSuffix("Propagator"));

  Acts::PropagatorOptions<Acts::ActionList<Acts::MaterialInteractor>,
                          Acts::AbortList<Acts::EndOfWorldReached>>
      extrapolationOptions(ctx.geoContext, ctx.magFieldContext);
  Acts::PropagatorOptions<> pOptions(ctx.geoContext, ctx.magFieldContext);

  Acts::FullBilloirVertexFitter::Config vertexFitterCfg;
  vertexFitterCfg.extractParameters
      .connect<&Acts::InputTrack::extractParameters>();

  // using TrackProxyType =
  // Acts::TrackContainer<Acts::ConstVectorTrackContainer,
  // Acts::ConstVectorMultiTrajectory, std::shared_ptr>::ConstTrackProxy;

  //////////////
  // MATCHING
  //////////////
  std::vector<std::pair<TrackProxyType, TrackProxyType>> trackPairs;
  int indexMS = 0;
  int indexVT = 0;

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCountsVT;
  std::vector<ParticleHitCount> particleHitCountsMS;

  // auto inputTracksVT = makeInputTracks(paramsVTinput);
  auto it1MS = m_inputTrackParametersMS.begin();
  auto it2MS = m_inputTrackContainerMS.begin();

  for (; it1MS != m_inputTrackParametersMS.end() &&
         it2MS != m_inputTrackContainerMS.end();
       ++it1MS, ++it2MS) {
    const auto& itrkparMS = *it1MS;
    const auto& itrkconMS = *it2MS;
    // for (const auto& itrkparMS : m_inputTrackParametersMS) {
    const auto& paramsMSinput = (*itrkparMS)(ctx);
    const auto& trackContainterMS = (*itrkconMS)(ctx);

    auto inputTracksMS = makeInputTracks(paramsMSinput);
    indexMS = 0;
    for (auto& trackMS : inputTracksMS) {
      indexVT = -1;
      bool isRec = false;
      std::pair<TrackProxyType, TrackProxyType> trackPairTmp(
          trackContainterMS.getTrack(indexMS),
          trackContainterMS.getTrack(indexMS));
      std::pair<int, int> trackPairIndex(indexMS, indexMS);

      // std::cout<<"Muon bef index: "<<indexMS<<std::endl;
      identifyContributingParticles(hitParticlesMapMS,
                                    trackContainterMS.getTrack(indexMS),
                                    particleHitCountsMS);
      // std::cout<<"Muon aft index: "<<indexMS<<std::endl;

      if (particleHitCountsMS.empty()) {
        ACTS_DEBUG(
            "No truth particle associated with this trajectory with tip index "
            "= "
            << trackContainterMS.getTrack(indexMS).tipIndex());
        continue;
      }
      ActsFatras::Barcode majorityParticleIdMS =
          particleHitCountsMS.front().particleId;
      ActsFatras::Barcode majorityParticleIdVT;
      ActsFatras::Barcode majorityParticleIdVTNoLoc;

      indexMS += 1;

      Acts::BoundTrackParameters params1 =
          vertexFitterCfg.extractParameters(trackMS);
      const auto resMS =
          extrapolator.propagateToSurface(params1, *pSurface, pOptions);

      if (!resMS.ok()) {
        continue;
      }
      float chi2 = 1e12;
      float chi2NoLoc = 1e12;
      const auto& endParamsMS = *resMS;

      Acts::BoundVector paramsMS = endParamsMS.parameters();
      const auto covMS = endParamsMS.covariance();

      Acts::ActsScalar covLoc0MS = (*covMS)(0, 0);
      Acts::ActsScalar covLoc1MS = (*covMS)(1, 1);
      Acts::ActsScalar covPhiMS = (*covMS)(2, 2);
      Acts::ActsScalar covThetaMS = (*covMS)(3, 3);
      Acts::ActsScalar covQOvPMS = (*covMS)(4, 4);

      Acts::ActsScalar loc0MS = paramsMS(Acts::BoundIndices::eBoundLoc0);
      Acts::ActsScalar loc1MS = paramsMS(Acts::BoundIndices::eBoundLoc1);
      Acts::ActsScalar phiMS = paramsMS(Acts::BoundIndices::eBoundPhi);
      Acts::ActsScalar thetaMS = paramsMS(Acts::BoundIndices::eBoundTheta);
      Acts::ActsScalar qOvPMS = paramsMS(Acts::BoundIndices::eBoundQOverP);

      Acts::ActsScalar loc0VTmatch = 0;
      Acts::ActsScalar loc1VTmatch = 0;
      Acts::ActsScalar phiVTmatch = 0;
      Acts::ActsScalar thetaVTmatch = 0;
      Acts::ActsScalar qOvPVTmatch = 0;

      Acts::ActsScalar covLoc0VTmatch = 0;
      Acts::ActsScalar covLoc1VTmatch = 0;
      Acts::ActsScalar covPhiVTmatch = 0;
      Acts::ActsScalar covThetaVTmatch = 0;
      Acts::ActsScalar covQOvPVTmatch = 0;

      auto it1VT = m_inputTrackParametersVT.begin();
      auto it2VT = m_inputTrackContainerVT.begin();
      int counter = 0;
      // auto inputTracksVT = makeInputTracks(paramsVTinput);
      for (; it1VT != m_inputTrackParametersVT.end() &&
             it2VT != m_inputTrackContainerVT.end();
           ++it1VT, ++it2VT) {
        const auto& itrkparVT = *it1VT;
        const auto& itrkconVT = *it2VT;
        // for (const auto& itrkparVT : m_inputTrackParametersVT) {
        const auto& paramsVTinput = (*itrkparVT)(ctx);
        const auto& trackContainterVT = (*itrkconVT)(ctx);
        std::cout<<"new container "<<counter<<std::endl;
        counter++;
        indexVT = -1;
        auto inputTracksVT = makeInputTracks(paramsVTinput);
        for (auto& trackVT : inputTracksVT) {
          indexVT++;
          Acts::BoundTrackParameters params2 =
              vertexFitterCfg.extractParameters(trackVT);
          const auto resVT =
              extrapolator.propagateToSurface(params2, *pSurface, pOptions);

          std::cout << "VT bef index: " << indexVT << " tipIdx: "
                    << trackContainterVT.getTrack(indexVT).tipIndex()
                    << std::endl;
          identifyContributingParticles(hitParticlesMapVT,
                                        trackContainterVT.getTrack(indexVT),
                                        particleHitCountsVT);
          // std::cout<<"VT aft index: "<<indexVT<<" tipIdx:
          // "<<trackContainterVT.getTrack(indexVT).tipIndex()<<std::endl;

          if (particleHitCountsVT.empty()) {
            ACTS_DEBUG(
                "No truth particle associated with this trajectory with tip "
                "index = "
                << trackContainterVT.getTrack(indexVT).tipIndex());
            std::cout
                << "No truth particle associated with this trajectory with "
                   "tip index = "
                << trackContainterVT.getTrack(indexVT).tipIndex() << std::endl;

            continue;
          }
          std::cout << "particle: " << particleHitCountsVT.front().particleId
                    << std::endl;

          // Acts::calculateTrackQuantities(trackVT);
          /*
          const auto resVT = Acts::extrapolateTrackToReferenceSurface(
              trackVT, *pSurface, extrapolator, extrapolationOptions,
              m_cfg.extrapolationStrategy, logger());
              */
          if (!resVT.ok()) {
            continue;
          }

          const auto& endParamsVT = *resVT;
          Acts::BoundVector paramsVT = endParamsVT.parameters();
          const auto covVT = endParamsVT.covariance();

          Acts::ActsScalar covLoc0VT = (*covVT)(0, 0);
          Acts::ActsScalar covLoc1VT = (*covVT)(1, 1);
          Acts::ActsScalar covPhiVT = (*covVT)(2, 2);
          Acts::ActsScalar covThetaVT = (*covVT)(3, 3);
          Acts::ActsScalar covQOvPVT = (*covVT)(4, 4);

          Acts::ActsScalar loc0VT = paramsVT(Acts::BoundIndices::eBoundLoc0);
          Acts::ActsScalar loc1VT = paramsVT(Acts::BoundIndices::eBoundLoc1);
          Acts::ActsScalar phiVT = paramsVT(Acts::BoundIndices::eBoundPhi);
          Acts::ActsScalar thetaVT = paramsVT(Acts::BoundIndices::eBoundTheta);
          Acts::ActsScalar qOvPVT = paramsVT(Acts::BoundIndices::eBoundQOverP);

          float chi2tmp =
              (loc0VT - loc0MS) * (loc0VT - loc0MS) / (covLoc0VT + covLoc0MS);
          chi2tmp +=
              (loc1VT - loc1MS) * (loc1VT - loc1MS) / (covLoc1VT + covLoc1MS);
          chi2tmp += (phiVT - phiMS) * (phiVT - phiMS) / (covPhiVT + covPhiMS);
          chi2tmp += (thetaVT - thetaMS) * (thetaVT - thetaMS) /
                     (covThetaVT + covThetaMS);
          chi2tmp +=
              (qOvPVT - qOvPMS) * (qOvPVT - qOvPMS) / (covQOvPVT + covQOvPMS);

          float chi2tmpNoLoc =
              (phiVT - phiMS) * (phiVT - phiMS) / (covPhiVT + covPhiMS);
          chi2tmp += (thetaVT - thetaMS) * (thetaVT - thetaMS) /
                     (covThetaVT + covThetaMS);
          chi2tmp +=
              (qOvPVT - qOvPMS) * (qOvPVT - qOvPMS) / (covQOvPVT + covQOvPMS);

          // write stuff for chi2 comparisons

          if (majorityParticleIdMS == particleHitCountsVT.front().particleId) {
            std::cout << "Match: chi2  = " << chi2tmp << std::endl;
            std::cout << "Match: chi2Noloc  = " << chi2tmpNoLoc << std::endl;
            isRec = true;
          } else {
            std::cout << "Fake: chi2  = " << chi2tmp << std::endl;
            std::cout << "Fake: chi2Noloc  = " << chi2tmpNoLoc << std::endl;
          }

          if (chi2tmpNoLoc < chi2NoLoc) {
            chi2NoLoc = chi2tmpNoLoc;
            majorityParticleIdVTNoLoc = particleHitCountsVT.front().particleId;
          }
          if (chi2tmp < chi2) {
            chi2 = chi2tmp;
            trackPairTmp.second = trackContainterVT.getTrack(indexVT);
            trackPairIndex.second = indexVT;
            majorityParticleIdVT = particleHitCountsVT.front().particleId;
            loc0VTmatch = loc0VT;
            loc1VTmatch = loc1VT;
            phiVTmatch = phiVT;
            thetaVTmatch = thetaVT;
            qOvPVTmatch = qOvPVT;
            covLoc0VTmatch = covLoc0VT;
            covLoc1VTmatch = covLoc1VT;
            covPhiVTmatch = covPhiVT;
            covThetaVTmatch = covThetaVT;
            covQOvPVTmatch = covQOvPVT;
          }
        }
      }
      if (chi2NoLoc != m_cfg.chi2max &&
          trackPairTmp.first != trackPairTmp.second && isRec) {
        if (majorityParticleIdVTNoLoc == majorityParticleIdMS) {
          std::cout << "NOLOCRECMATCHED4!!!!!  p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        } else {
          std::cout << "NOLOCRECFAKE MATCH4!!! p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        }
      }

      if (chi2NoLoc != m_cfg.chi2max &&
          trackPairTmp.first != trackPairTmp.second) {
        if (majorityParticleIdVTNoLoc == majorityParticleIdMS) {
          std::cout << "NOLOCMATCHED3!!!!!  p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        } else {
          std::cout << "NOLOCFAKE MATCH3!!! p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        }
      }

      if (chi2 != m_cfg.chi2max && trackPairTmp.first != trackPairTmp.second &&
          isRec) {
        if (majorityParticleIdVT == majorityParticleIdMS) {
          std::cout << "RECMATCHED2!!!!!  p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        } else {
          std::cout << "RECFAKE MATCH2!!! p = " << 1. / qOvPMS << ", "
                    << thetaMS << ", " << phiMS << std::endl;
        }
      }

      if (chi2 != m_cfg.chi2max && trackPairTmp.first != trackPairTmp.second) {
        trackPairs.push_back(trackPairTmp);
        if (majorityParticleIdVT == majorityParticleIdMS) {
          std::cout << "MATCHED1!!!!!  p = " << 1. / qOvPMS << ", " << thetaMS
                    << ", " << phiMS << std::endl;
          std::cout << "chi2/ndf " << chi2 << std::endl;
          std::cout << "loc0MS: " << loc0MS << " loc0VT: " << loc0VTmatch
                    << std::endl;
          std::cout << "loc1MS: " << loc1MS << " loc1VT: " << loc1VTmatch
                    << std::endl;
          std::cout << "phiMS: " << phiMS << " phiVT: " << phiVTmatch
                    << std::endl;
          std::cout << "thetaMS: " << thetaMS << " thetaVT: " << thetaVTmatch
                    << std::endl;
          std::cout << "qOvPMS: " << qOvPMS << " qOvPVT: " << qOvPVTmatch
                    << std::endl;

          std::cout << "dloc0VT : "
                    << (loc0VTmatch - loc0MS) * (loc0VTmatch - loc0MS) /
                           (covLoc0VTmatch + covLoc0MS)
                    << std::endl;
          std::cout << "dloc1VT : "
                    << (loc1VTmatch - loc1MS) * (loc1VTmatch - loc1MS) /
                           (covLoc1VTmatch + covLoc1MS)
                    << std::endl;
          std::cout << "dphiVT  : "
                    << (phiVTmatch - phiMS) * (phiVTmatch - phiMS) /
                           (covPhiVTmatch + covPhiMS)
                    << std::endl;
          std::cout << "dthetaVT: "
                    << (thetaVTmatch - thetaMS) * (thetaVTmatch - thetaMS) /
                           (covThetaVTmatch + covThetaMS)
                    << std::endl;
          std::cout << "dqOvPVT : "
                    << (qOvPVTmatch - qOvPMS) * (qOvPVTmatch - qOvPMS) /
                           (covQOvPVTmatch + covQOvPMS)
                    << std::endl;

          std::cout << "MS: " << majorityParticleIdMS << std::endl;
          std::cout << "VT: " << majorityParticleIdVT << std::endl;
        } else {
          std::cout << "FAKE MATCH1!!! p = " << 1. / qOvPMS << ", " << thetaMS
                    << ", " << phiMS << std::endl;
          std::cout << "chi2/ndf " << chi2 << std::endl;
          std::cout << "loc0MS: " << loc0MS << " loc0VT: " << loc0VTmatch
                    << std::endl;
          std::cout << "loc1MS: " << loc1MS << " loc1VT: " << loc1VTmatch
                    << std::endl;
          std::cout << "phiMS: " << phiMS << " phiVT: " << phiVTmatch
                    << std::endl;
          std::cout << "thetaMS: " << thetaMS << " thetaVT: " << thetaVTmatch
                    << std::endl;
          std::cout << "qOvPMS: " << qOvPMS << " qOvPVT: " << qOvPVTmatch
                    << std::endl;

          std::cout << "dloc0VT : "
                    << (loc0VTmatch - loc0MS) * (loc0VTmatch - loc0MS) /
                           (covLoc0VTmatch + covLoc0MS)
                    << std::endl;
          std::cout << "dloc1VT : "
                    << (loc1VTmatch - loc1MS) * (loc1VTmatch - loc1MS) /
                           (covLoc1VTmatch + covLoc1MS)
                    << std::endl;
          std::cout << "dphiVT  : "
                    << (phiVTmatch - phiMS) * (phiVTmatch - phiMS) /
                           (covPhiVTmatch + covPhiMS)
                    << std::endl;
          std::cout << "dthetaVT: "
                    << (thetaVTmatch - thetaMS) * (thetaVTmatch - thetaMS) /
                           (covThetaVTmatch + covThetaMS)
                    << std::endl;
          std::cout << "dqOvPVT : "
                    << (qOvPVTmatch - qOvPMS) * (qOvPVTmatch - qOvPMS) /
                           (covQOvPVTmatch + covQOvPMS)
                    << std::endl;

          std::cout << "MS: " << majorityParticleIdMS << std::endl;
          std::cout << "VT: " << majorityParticleIdVT << std::endl;
        }
      }
    }
  }

  //////////////
  // REFIT PSEUDOCODE
  //////////////

  // Perform the fit for each input track
  std::cout << "START REFITTING" << std::endl;
  std::vector<Acts::SourceLink> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;
  RefittingCalibrator calibrator;
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  /*
  std::cout<<"MY REFITTING"<<std::endl;
  for (auto& pair : trackPairs) {
    auto trackMS = pair.first;
    auto trackVT = pair.second;

    trackSourceLinks.clear();
    surfSequence.clear();

    TrackFitterFunction::GeneralFitterOptions options{
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        &trackVT.referenceSurface(), Acts::PropagatorPlainOptions()};

    const Acts::BoundTrackParameters initialParams(
        trackVT.referenceSurface().getSharedPtr(), trackVT.parameters(),
        trackVT.covariance(), trackVT.particleHypothesis());

    for (auto state : trackVT.trackStatesReversed()) {
      surfSequence.push_back(&state.referenceSurface());
      auto sl = RefittingCalibrator::RefittingSourceLink{state};
      trackSourceLinks.push_back(Acts::SourceLink{sl});
    }
    */
  /*
  for (auto state : trackMS.trackStatesReversed()) {
    surfSequence.push_back(&state.referenceSurface());
    auto sl = RefittingCalibrator::RefittingSourceLink{state};
    trackSourceLinks.push_back(Acts::SourceLink{sl});
  }
  */
  /*
   if (surfSequence.empty()) {
     ACTS_WARNING("Empty track found.");
     continue;
   }
   //auto result = (*m_cfg.fit)(trackSourceLinks, initialParams);
   //auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options);
   //auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
 calibrator);
   //auto result = (*m_cfg.fit)(trackSourceLinks, initialParams, options,
 calibrator, surfSequence); auto result = (*m_cfg.fit)(trackSourceLinks,
 initialParams, options, calibrator, surfSequence, tracks);


   if (result.ok()) {
     // Get the fit output object
     const auto& refittedTrack = result.value();
     if (refittedTrack.hasReferenceSurface()) {
       ACTS_VERBOSE("Refitted parameters for track ");
       ACTS_VERBOSE("  " << trackMS.parameters().transpose());
     } else {
       ACTS_DEBUG("No refitted parameters for track ");
     }
   } else {
     ACTS_WARNING("Fit failed for track with error: " << result.error() << ", "
                  << result.error().message());
   }

 }
 //m_outputTrackParameters(ctx, std::move(outputTrackParameter));
 ConstTrackContainer constTracks{
     std::make_shared<Acts::ConstVectorTrackContainer>(
         std::move(*trackContainer)),
     std::make_shared<Acts::ConstVectorMultiTrajectory>(
         std::move(*trackStateContainer))};

 m_outputTracks(ctx, std::move(constTracks));
 */

  m_outputMatchedTracks(ctx, std::move(trackPairs));

  return ActsExamples::ProcessCode::SUCCESS;
}