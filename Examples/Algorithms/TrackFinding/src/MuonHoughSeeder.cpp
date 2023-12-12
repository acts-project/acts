// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MuonHoughSeeder.hpp"
#include "ActsExamples/EventData/MuonSimHit.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <variant>

namespace{
  static const std::map<int, std::string> stationDict{
        {0, "BIL"}, {1, "BIS"}, {7, "BIR"},
        {2, "BML"}, {3, "BMS"}, {8, "BMF"}, {53, "BME"}, {54, "BMG"}, {52, "BIM"},
        {4, "BOL"}, {5, "BOS"}, {9, "BOF"}, {10, "BOG"},
        {6, "BEE"}, {14, "EEL"}, {15, "EES"}, 
        {13, "EIL"}, {49, "EIS"},
        {17, "EML"}, {18, "EMS"}, 
        {20, "EOL"}, {21, "EOS"}
    };

}

ActsExamples::MuonHoughSeeder::MuonHoughSeeder(
    ActsExamples::MuonHoughSeeder::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("MuonHoughSeeder", lvl),
      m_cfg(std::move(cfg)),
      m_logger(Acts::getDefaultLogger("MuonHoughSeeder", lvl)) {
  // require spacepoints or input measurements (or both), but at least one kind
  // of input
  bool foundInput = false;

  if (m_cfg.inDriftCircles.empty()) {
    throw std::invalid_argument(
        "MuonHoughSeeder: Missing drift circle collection");
  }
  if (m_cfg.inSimHits.empty()) {
    throw std::invalid_argument(
        "MuonHoughSeeder: Missing sim hit collection");
  }

  m_inputDriftCircles.initialize(m_cfg.inDriftCircles);
  m_inputSimHits.initialize(m_cfg.inSimHits);

}
using PatternSeed = std::pair<double, double>; // y0, tan theta
ActsExamples::ProcessCode ActsExamples::MuonHoughSeeder::execute(
    const AlgorithmContext& ctx) const {
  // std::cout << "hi!"<<std::endl;
  // ACTS_WARNING("HI!!!");

  auto gotSH = m_inputSimHits(ctx);
  auto gotDC = m_inputDriftCircles(ctx);

  Acts::HoughTransformUtils::HoughPlaneConfig planeCfg;
  planeCfg.fracForWidth = 0.6;
  planeCfg.nBinsX = 1000;
  planeCfg.nBinsY = 1000;
  planeCfg.nLayers = 8;
  planeCfg.threshold = 3;
  planeCfg.xMin = -3;
  planeCfg.xMax = 3;
  planeCfg.yMin = -2000;
  planeCfg.yMax = 2000;
  planeCfg.localMaxWindowSize = 50; 

  auto houghParam_fromDC_left = [](double tanTheta, const DriftCircle& DC){
    return DC.y() - tanTheta * DC.z() - DC.rDrift() / std::cos(std::atan(tanTheta));
  };
  auto houghParam_fromDC_right = [](double tanTheta, const DriftCircle& DC){
    return DC.y() - tanTheta * DC.z() + DC.rDrift() / std::cos(std::atan(tanTheta));
  };
  auto houghWidth_fromDC = [](double tanTheta, const DriftCircle& DC){
    return DC.rDriftError() * 3.;
  };
  auto id_fromDC = [](const DriftCircle& DC)->Acts::HoughTransformUtils::idType {
    muonMdtIdentifierFields idf;
    idf.multilayer = DC.multilayer();
    idf.stationEta = DC.stationEta();
    idf.stationPhi = DC.stationPhi();
    idf.stationName = DC.stationName();
    idf.tubeLayer = DC.tubeLayer();
    idf.tube = DC.tube();
    return compressId(idf);
  };
  auto layer_from_DC = [](const DriftCircle& DC) -> unsigned{
    return 3 * (DC.multilayer()-1) + (DC.tubeLayer()-1); 
  };


  std::vector<PatternSeed> truePatterns; 
  // std::map<
  for (auto & SH : gotSH){
    muonMdtIdentifierFields detailedInfo = ActsExamples::splitId(SH.geometryId().value()); 
    // double theta = std::atan(SH.direction().y() / SH.direction().z());
    truePatterns.emplace_back(SH.fourPosition().y(), SH.direction().y() / SH.direction().z()); 
    // auto ID = SH.
    // std::cout <<" saw a true sim hit in station "<<(int)detailedInfo.stationName<<", eta "<<(int)detailedInfo.stationEta<<", phi "<<(int)detailedInfo.stationPhi<<" @ y0,tantheta = "<<truePatterns.back().first<<", "<<truePatterns.back().second<<std::endl; 

    m_houghHist->Reset();
    Acts::HoughTransformUtils::HoughPlane houghPlane(planeCfg); 

    for (DriftCircle & DC : gotDC){
      if (  DC.stationEta() == detailedInfo.stationEta
          &&DC.stationPhi() == detailedInfo.stationPhi
          &&DC.stationName() == detailedInfo.stationName){
              // std::cout <<"   --> DC station "<<(int)DC.stationName()<<", eta "<<(int)DC.stationEta()<<", phi "<<(int)DC.stationPhi()<<" @ y,z "<<DC.y()<<", "<<DC.z()<<" with r = "<<DC.rDrift()<<std::endl;
              houghPlane.fill<DriftCircle>(DC,houghParam_fromDC_left,houghWidth_fromDC, id_fromDC, layer_from_DC);
              houghPlane.fill<DriftCircle>(DC,houghParam_fromDC_right,houghWidth_fromDC, id_fromDC, layer_from_DC);
              for (int thetaBin = 1; thetaBin < m_houghHist->GetNbinsY()+1; ++thetaBin){
                double tanTheta = m_houghHist->GetYaxis()->GetBinCenter(thetaBin); 
                for (int factor = -1; factor < 2; factor+=2){
                  double y0Main = DC.y() - tanTheta * DC.z() + factor* DC.rDrift() / std::cos(std::atan(tanTheta)); 
                  double nomBin = m_houghHist->GetXaxis()->FindBin(y0Main);
                  // std::cout <<" tanTheta is "<<tanTheta<<" and y0 "<<y0Main<<std::endl;
                  double tubeRadius = std::max(DC.rDriftError(), 1.) * 3 ;  
                  for (int bX =  nomBin; (bX == nomBin) || (std::abs(m_houghHist->GetXaxis()->GetBinCenter(bX) - y0Main) < tubeRadius && bX >0); --bX){
                    m_houghHist->Fill(m_houghHist->GetXaxis()->GetBinCenter(bX),tanTheta);
                  }
                  for (int bX =  nomBin; std::abs(m_houghHist->GetXaxis()->GetBinCenter(bX) - y0Main) < tubeRadius && bX < m_houghHist->GetNbinsX()+1; ++bX){
                    m_houghHist->Fill(m_houghHist->GetXaxis()->GetBinCenter(bX),tanTheta);
                  }
                }
              }
          }
    }
    m_outCanvas->SetTitle(Form("Station %s, Eta %i, Phi %i",stationDict.at(detailedInfo.stationName).c_str(),(int)detailedInfo.stationEta,(int)detailedInfo.stationPhi)); 
    m_houghHist->SetTitle(Form("Station %s, Eta %i, Phi %i",stationDict.at(detailedInfo.stationName).c_str(),(int)detailedInfo.stationEta,(int)detailedInfo.stationPhi)); 
    m_outCanvas->cd();
    m_houghHist->SetContour(m_houghHist->GetMaximum()+1);
    for (size_t k = 0; k < m_houghHist->GetMaximum()+1; ++k){
      m_houghHist->SetContourLevel(k, k - 0.5);
    }
    double m = m_houghHist->GetMaximumBin();
    double tolerance = m_houghHist->GetMaximum()-1; 
    std::vector<std::unique_ptr<TMarker>> markers; 
    m_houghHist->SetContourLevel(m_houghHist->GetMaximum()+1, m_houghHist->GetMaximum()+0.5);
    m_houghHist->Draw("COLZ");
    for (int bx = 0; bx < m_houghHist->GetXaxis()->GetNbins(); ++bx){
      for (int by = 0; by < m_houghHist->GetYaxis()->GetNbins(); ++by){
        if (m_houghHist->GetBinContent(bx,by) >= tolerance && m_houghHist->GetBinContent(bx,by) > 2){
            auto PeakPos = std::make_unique<TMarker>(m_houghHist->GetXaxis()->GetBinCenter(bx),
                                            m_houghHist->GetYaxis()->GetBinCenter(by),kFullSquare);
            PeakPos->SetMarkerSize(0.3);
            PeakPos->SetMarkerColor(kBlue);
            markers.push_back(std::move(PeakPos));
            markers.back()->Draw();
        }
      }
    }
    auto trueMarker = std::make_unique<TMarker>(truePatterns.back().first, truePatterns.back().second,kOpenCrossX);
    trueMarker->SetMarkerSize(3);
    trueMarker->SetMarkerColor(kRed);
    trueMarker->Draw();

    auto maxima = houghPlane.getMaxima();
    std::cout <<" found "<<maxima.size()<<" maxima! "<<std::endl;
    for (auto & max : maxima){
        markers.push_back(std::make_unique<TMarker>(max.x, max.y,kOpenCircle));
        markers.back()->SetMarkerSize(3);
        markers.back()->SetMarkerColor(kOrange-3);
        markers.back()->Draw();
        std::cout <<"  --> at "<<max.x<<" , "<<max.y<<std::endl;
    }
    m_outCanvas->SaveAs("HoughHistograms.pdf");

    TH2D hJahred("jahred","hatAuchEinHough;x;y",planeCfg.nBinsY, planeCfg.yMin, planeCfg.yMax, planeCfg.nBinsX, planeCfg.xMin, planeCfg.xMax);
     for (int bx = 0; bx < hJahred.GetXaxis()->GetNbins(); ++bx){
      for (int by = 0; by < hJahred.GetYaxis()->GetNbins(); ++by){
        hJahred.SetBinContent(by+1, bx+1, houghPlane(bx,by));
      }
     }
     hJahred.Draw("COLZ"); 
    for (auto & max : maxima){
        markers.push_back(std::make_unique<TMarker>(max.y, max.x,kOpenCircle));
        markers.back()->SetMarkerSize(3);
        markers.back()->SetMarkerColor(kOrange-3);
        markers.back()->Draw();
        std::cout <<"  --> at "<<max.x<<" , "<<max.y<<std::endl;
    }
    trueMarker->Draw();
    m_outCanvas->SaveAs("HoughHistograms.pdf");

  }

  ACTS_WARNING("SH: "<<gotSH.size()); 
  ACTS_WARNING("DC: "<<gotDC.size()); 
  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::MuonHoughSeeder::initialize(){
  m_outCanvas = std::make_unique<TCanvas>("canvas","",800,800);
  m_outCanvas->SaveAs("HoughHistograms.pdf[");
  m_outCanvas->SetRightMargin(0.12);
  m_outCanvas->SetLeftMargin(0.12); 
  m_houghHist = std::make_unique<TH2D>("hallo",";y0;tanTheta;density",1000,-2000,2000,1000,-3,3);
  m_houghHist->SetDirectory(0);
  gStyle->SetPalette(kGreyScale);
  gStyle->SetOptStat(0);
  return ProcessCode::SUCCESS;
}
ActsExamples::ProcessCode ActsExamples::MuonHoughSeeder::finalize(){
  m_outCanvas->SaveAs("HoughHistograms.pdf]");
  return ProcessCode::SUCCESS;
}
