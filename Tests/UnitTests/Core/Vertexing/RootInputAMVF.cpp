// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/GridDensityVertexFinder.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include "AMVFTestData.ipp"

#include <chrono>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

namespace Acts {
namespace Test {

bool runLocally = true;
TString inputString = "";
TString outputString = "";

int nMaxEvents = 10;

using namespace Acts::UnitLiterals;

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

std::vector<std::pair<Vertex<BoundParameters>, std::vector<BoundParameters>>>
readRootFile();

BOOST_AUTO_TEST_CASE(RootInputAMVF) {
  if (runLocally) {
    outputString = "/Users/bschlag/atlas/acts-core/build/ACTS.root";
  } else {
    outputString =
        "/afs/cern.ch/user/b/baschlag/work/private/acts-core/ACTS.root";
  }

  // OUTPUT
  TFile outFile(outputString, "RECREATE");
  TTree* outtree = new TTree("mytree", "An example of a ROOT tree");

  std::vector<double>* m_vx = 0;
  std::vector<double>* m_vy = 0;
  std::vector<double>* m_vz = 0;
  std::vector<double>* m_vCovx = 0;
  std::vector<double>* m_vCovy = 0;
  std::vector<double>* m_vCovz = 0;
  std::vector<int>* m_nTracks = 0;

  double m_time = 0;

  outtree->Branch("time", &m_time);

  outtree->Branch("vx", &m_vx);
  outtree->Branch("vy", &m_vy);
  outtree->Branch("vz", &m_vz);

  outtree->Branch("vCovx", &m_vCovx);
  outtree->Branch("vCovy", &m_vCovy);
  outtree->Branch("vCovz", &m_vCovz);
  outtree->Branch("nTracks", &m_nTracks);

  // auto events = getAthenaTracks();

  auto events = readRootFile();

  // Set debug mode
  bool debugMode = false;
  // Set up constant B-Field
  ConstantBField bField(Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  // EigenStepper<ConstantBField> stepper(bField);
  EigenStepper<ConstantBField> stepper(bField);

  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(stepper);

  // IP 3D Estimator
  using IPEstimator = ImpactPointEstimator<BoundParameters, Propagator>;

  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  AnnealingUtility::Config annealingConfig(temperatures);
  AnnealingUtility annealingUtility(annealingConfig);

  using Fitter = AdaptiveMultiVertexFitter<BoundParameters, Linearizer>;

  Fitter::Config fitterCfg(ipEstimator);

  fitterCfg.annealingTool = annealingUtility;

  // Linearizer for BoundParameters type test
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Test smoothing
  fitterCfg.doSmoothing = true;

  Fitter fitter(fitterCfg);

  using SeedFinder = GridDensityVertexFinder<2000, 35>;
  SeedFinder::Config seedFinderCfg;
  seedFinderCfg.cacheGridStateForTrackRemoval = true;
  SeedFinder seedFinder(seedFinderCfg);

  using Finder = AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer);

  // TODO: test this as well!
  // finderConfig.useBeamSpotConstraint = false;

  finderConfig.refitAfterBadVertex = false;

  Finder finder(finderConfig);

  Finder::State state;

  VertexingOptions<BoundParameters> vertexingOptions(tgContext, mfContext);

  std::vector<BoundParameters> tracks;
  for (auto evt : events) {
    m_vx->clear();
    m_vy->clear();
    m_vz->clear();
    m_vCovx->clear();
    m_vCovy->clear();
    m_vCovz->clear();
    m_nTracks->clear();

    tracks = evt.second;

    std::vector<const BoundParameters*> tracksPtr;
    for (const auto& trk : tracks) {
      tracksPtr.push_back(&trk);
    }

    Vertex<BoundParameters> constraintVtx = evt.first;

    vertexingOptions.vertexConstraint = constraintVtx;

    auto now11 = std::chrono::system_clock::now();
    auto findResult = finder.find(tracksPtr, vertexingOptions, state);
    auto now12 = std::chrono::system_clock::now();
    auto ms1 =
        std::chrono::duration_cast<std::chrono::milliseconds>(now12 - now11)
            .count();

    std::cout << "time needed: " << ms1 << std::endl;
    m_time = ms1;

    if (!findResult.ok()) {
      std::cout << findResult.error().message() << std::endl;
    }

    std::vector<Vertex<BoundParameters>> allVertices = *findResult;

    std::cout << "reco vertices: " << allVertices.size() << std::endl;

    for (auto vtx : allVertices) {
      auto pos = vtx.position();
      auto cov = vtx.covariance();
      auto nTracks = vtx.tracks().size();

      m_vx->push_back(pos[0]);
      m_vy->push_back(pos[1]);
      m_vz->push_back(pos[2]);

      m_vCovx->push_back(cov(0, 0));
      m_vCovy->push_back(cov(1, 1));
      m_vCovz->push_back(cov(2, 2));

      m_nTracks->push_back(nTracks);
    }

    outtree->Fill();
  }

  outFile.cd();
  outtree->Write();

  outFile.Close();
}

std::vector<std::pair<Vertex<BoundParameters>, std::vector<BoundParameters>>>
readRootFile() {
  std::vector<std::pair<Vertex<BoundParameters>, std::vector<BoundParameters>>>
      events;

  if (runLocally) {
    inputString =
        "/Users/bschlag/atlas/acts-core/Tests/UnitTests/Core/Vertexing/"
        "inputfile_1000events.root";
  } else {
    inputString =
        "/afs/cern.ch/user/b/baschlag/work/private/acts-core/"
        "amvf_athena_output.root";
  }

  TFile fileIn(inputString);

  TTree* theTree = nullptr;
  fileIn.GetObject("VertexAnalysis", theTree);

  int nEvents = theTree->GetEntries();

  int m_nVtx = 0;
  /// The track parameter
  std::vector<double>* m_d0 = 0;
  std::vector<double>* m_z0 = 0;
  std::vector<double>* m_phi = 0;
  std::vector<double>* m_theta = 0;
  std::vector<double>* m_qp = 0;
  std::vector<int>* m_vtxID = 0;

  /// The track covariance matrix
  std::vector<double>* m_cov11 = 0;
  std::vector<double>* m_cov12 = 0;
  std::vector<double>* m_cov13 = 0;
  std::vector<double>* m_cov14 = 0;
  std::vector<double>* m_cov15 = 0;

  std::vector<double>* m_cov21 = 0;
  std::vector<double>* m_cov22 = 0;
  std::vector<double>* m_cov23 = 0;
  std::vector<double>* m_cov24 = 0;
  std::vector<double>* m_cov25 = 0;

  std::vector<double>* m_cov31 = 0;
  std::vector<double>* m_cov32 = 0;
  std::vector<double>* m_cov33 = 0;
  std::vector<double>* m_cov34 = 0;
  std::vector<double>* m_cov35 = 0;

  std::vector<double>* m_cov41 = 0;
  std::vector<double>* m_cov42 = 0;
  std::vector<double>* m_cov43 = 0;
  std::vector<double>* m_cov44 = 0;
  std::vector<double>* m_cov45 = 0;

  std::vector<double>* m_cov51 = 0;
  std::vector<double>* m_cov52 = 0;
  std::vector<double>* m_cov53 = 0;
  std::vector<double>* m_cov54 = 0;
  std::vector<double>* m_cov55 = 0;

  // beamspot
  std::vector<double>* m_bsPosX = 0;
  std::vector<double>* m_bsPosY = 0;
  std::vector<double>* m_bsPosZ = 0;

  std::vector<double>* m_bsCovXX = 0;
  std::vector<double>* m_bsCovYY = 0;
  std::vector<double>* m_bsCovZZ = 0;

  theTree->SetBranchAddress("nVtx", &m_nVtx);

  theTree->SetBranchAddress("d0", &m_d0);
  theTree->SetBranchAddress("z0", &m_z0);
  theTree->SetBranchAddress("phi", &m_phi);
  theTree->SetBranchAddress("theta", &m_theta);
  theTree->SetBranchAddress("qp", &m_qp);
  theTree->SetBranchAddress("vtxID", &m_vtxID);

  theTree->SetBranchAddress("cov11", &m_cov11);
  theTree->SetBranchAddress("cov12", &m_cov12);
  theTree->SetBranchAddress("cov13", &m_cov13);
  theTree->SetBranchAddress("cov14", &m_cov14);
  theTree->SetBranchAddress("cov15", &m_cov15);

  theTree->SetBranchAddress("cov21", &m_cov21);
  theTree->SetBranchAddress("cov22", &m_cov22);
  theTree->SetBranchAddress("cov23", &m_cov23);
  theTree->SetBranchAddress("cov24", &m_cov24);
  theTree->SetBranchAddress("cov25", &m_cov25);

  theTree->SetBranchAddress("cov31", &m_cov31);
  theTree->SetBranchAddress("cov32", &m_cov32);
  theTree->SetBranchAddress("cov33", &m_cov33);
  theTree->SetBranchAddress("cov34", &m_cov34);
  theTree->SetBranchAddress("cov35", &m_cov35);

  theTree->SetBranchAddress("cov41", &m_cov41);
  theTree->SetBranchAddress("cov42", &m_cov42);
  theTree->SetBranchAddress("cov43", &m_cov43);
  theTree->SetBranchAddress("cov44", &m_cov44);
  theTree->SetBranchAddress("cov45", &m_cov45);

  theTree->SetBranchAddress("cov51", &m_cov51);
  theTree->SetBranchAddress("cov52", &m_cov52);
  theTree->SetBranchAddress("cov53", &m_cov53);
  theTree->SetBranchAddress("cov54", &m_cov54);
  theTree->SetBranchAddress("cov55", &m_cov55);

  theTree->SetBranchAddress("bsPosX", &m_bsPosX);
  theTree->SetBranchAddress("bsPosY", &m_bsPosY);
  theTree->SetBranchAddress("bsPosZ", &m_bsPosZ);

  theTree->SetBranchAddress("bsCovX", &m_bsCovXX);
  theTree->SetBranchAddress("bsCovY", &m_bsCovYY);
  theTree->SetBranchAddress("bsCovZ", &m_bsCovZZ);

  if (nEvents > nMaxEvents) {
    nEvents = nMaxEvents;
  }

  std::cout << nEvents << " events to run over ..." << std::endl;
  for (int iEvent = 0; iEvent < nEvents; iEvent++) {
    theTree->GetEvent(iEvent);

    std::cout << iEvent << ". event, nvertices: " << m_nVtx << std::endl;

    int nTracks = m_cov33->size();
    std::vector<BoundParameters> tracks;

    std::shared_ptr<PerigeeSurface> perigeeSurface =
        Surface::makeShared<PerigeeSurface>(
            Vector3D((*m_bsPosX)[0], (*m_bsPosY)[0], (*m_bsPosZ)[0]));

    for (int i = 0; i < nTracks; i++) {
      BoundVector params;
      params << (*m_d0)[i], (*m_z0)[i], (*m_phi)[i], (*m_theta)[i],
          (*m_qp)[i] * 1. / (1_MeV), 0;

      Covariance covMat;
      covMat << (*m_cov11)[i], (*m_cov12)[i], (*m_cov13)[i], (*m_cov14)[i],
          (*m_cov15)[i] * 1. / (1_MeV), 0, (*m_cov21)[i], (*m_cov22)[i],
          (*m_cov23)[i], (*m_cov24)[i], (*m_cov25)[i] * 1. / (1_MeV), 0,
          (*m_cov31)[i], (*m_cov32)[i], (*m_cov33)[i], (*m_cov34)[i],
          (*m_cov35)[i] * 1. / (1_MeV), 0, (*m_cov41)[i], (*m_cov42)[i],
          (*m_cov43)[i], (*m_cov44)[i], (*m_cov45)[i] * 1. / (1_MeV), 0,
          (*m_cov51)[i] * 1. / (1_MeV), (*m_cov52)[i] * 1. / (1_MeV),
          (*m_cov53)[i] * 1. / (1_MeV), (*m_cov54)[i] * 1. / (1_MeV),
          (*m_cov55)[i] * 1. / (1_MeV), 0, 0, 0, 0, 0, 0, 1;

      auto boundParams =
          BoundParameters(tgContext, covMat, params, perigeeSurface);

      tracks.push_back(boundParams);
    }

    ActsSymMatrixD<3> cov = ActsSymMatrixD<3>::Zero();
    cov(0, 0) = (*m_bsCovXX)[0];
    cov(1, 1) = (*m_bsCovYY)[0];
    cov(2, 2) = (*m_bsCovZZ)[0];

    Vector3D pos = Vector3D((*m_bsPosX)[0], (*m_bsPosY)[0], (*m_bsPosZ)[0]);

    Vertex<BoundParameters> constraintVtx(pos);
    constraintVtx.setCovariance(cov);

    events.push_back({constraintVtx, tracks});
  }
  return events;
}

}  // namespace Test
}  // namespace Acts