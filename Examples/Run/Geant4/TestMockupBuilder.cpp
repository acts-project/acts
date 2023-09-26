// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MuonSpectrometerMockupDetector/MockupSectorBuilder.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <memory>
#include <limits>
#include <chrono>

using namespace Acts;
using namespace Acts::Experimental;
using namespace ActsExamples;

using namespace Acts::UnitLiterals;


struct StrawSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Surface& sf) const {
    return (sf.type() == Surface::Straw);
  }
};

int main() {

  //Measure the cpu time elapsed
  auto begin = std::chrono::high_resolution_clock::now();

  Acts::GeometryContext tContext;
  Acts::MagneticFieldContext mContext;

  auto mockup_config = MockupSectorBuilder::Config();
  mockup_config.gdmlPath =
      " ../../../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
      "MuonChamber.gdml";
  mockup_config.NumberOfSectors = 1;
  mockup_config.robustMode = true;

  MockupSectorBuilder mockup_builder(mockup_config);

  // ***Inner Chamber***

 // configure the multi layers for the inner chamber

  MockupSectorBuilder::MultiLayerConfig multilayer1_config{"MultiLayer1_InnerChamber", {"Inner_Skin_ml0"}, {"xx"}};
  MockupSectorBuilder::MultiLayerConfig multilayer2_config{"MultiLayer2_InnerChamber", {"Inner_Skin_ml1"}, {"xx"}};

  auto ml1_detcomp = mockup_builder.buildMultiLayer(multilayer1_config);
  auto ml2_detcomp = mockup_builder.buildMultiLayer(multilayer2_config);

  std::vector<std::shared_ptr<DetectorVolume>> ml_volumes = {ml1_detcomp->volumes.front(), ml2_detcomp->volumes.front()};

  // configure the inner chamber
  MockupSectorBuilder::ChamberConfig mockup_chamberConfig_inner{"Inner_Chamber", ml_volumes};

  auto inner_chamber_detcomp = mockup_builder.buildChamber(mockup_chamberConfig_inner);

  // ***Middle Chamber***

  // configure the multi layers for the middle chamber

  multilayer1_config.name = "MultiLayer1_MiddleChamber";
  multilayer1_config.SensitiveNames = {"Middle_Skin_ml0"};
 

  multilayer2_config.name = "MultiLayer2_MiddleChamber";
  multilayer2_config.SensitiveNames = {"Middle_Skin_ml1"};


  ml1_detcomp = mockup_builder.buildMultiLayer(multilayer1_config);
  ml2_detcomp = mockup_builder.buildMultiLayer(multilayer2_config);

  ml_volumes.clear();
  ml_volumes = {ml1_detcomp->volumes.front(), ml2_detcomp->volumes.front()};

  // configure the middle chamber
 MockupSectorBuilder::ChamberConfig mockup_chamberConfig_middle{"Middle_Chamber", ml_volumes};

  auto middle_chamber_detcomp = mockup_builder.buildChamber(mockup_chamberConfig_middle);

  //***Outer Chamber***

  // configure the multi layers for the outer chamber

  multilayer1_config.name = "MultiLayer1_OuterChamber";
  multilayer1_config.SensitiveNames = {"Outer_Skin_ml0"};
  

  multilayer2_config.name = "MultiLayer2_OuterChamber";
  multilayer2_config.SensitiveNames = {"Outer_Skin_ml1"};


  ml1_detcomp = mockup_builder.buildMultiLayer(multilayer1_config);
  ml2_detcomp = mockup_builder.buildMultiLayer(multilayer2_config);

  ml_volumes.clear();
  ml_volumes = {ml1_detcomp->volumes.front(), ml2_detcomp->volumes.front()};

  // configure the outer chamber
  MockupSectorBuilder::ChamberConfig mockup_chamberConfig_outer{"Outer_Chamber", ml_volumes};

  auto outer_chamber_detcomp = mockup_builder.buildChamber(mockup_chamberConfig_outer);

  // Construct the detector component sector

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> chamber_volumes = {inner_chamber_detcomp->volumes.front(), 
  middle_chamber_detcomp->volumes.front(), outer_chamber_detcomp->volumes.front()};

  auto detector_sector = mockup_builder.buildSector(chamber_volumes);

  //Navigate through the sector detector

   auto detector = Acts::Experimental::Detector::makeShared(
      "Detector", {detector_sector->volumes.front()}, Acts::Experimental::tryRootVolumes());

     mockup_builder.drawSector(detector_sector->volumes.front(), "sector_volume.obj");

  using StrawsCollector = SurfaceCollector<StrawSelector>;

  using ActionListType = Acts::ActionList<StrawsCollector>;
  using AbortListType = Acts::AbortList<>;

  auto bField = std::make_shared<Acts::ConstantBField>(
      Acts::Vector3(2* Acts::UnitConstants::T, 0, 0));

  Acts::Experimental::DetectorNavigator::Config navCfg;
  navCfg.detector = detector.get();

  auto stepper = Acts::EigenStepper<>(bField);
  auto navigator = Acts::Experimental::DetectorNavigator(
      navCfg, Acts::getDefaultLogger("DetectorNavigator",
                                     Acts::Logging::Level::VERBOSE));
  auto options = Acts::PropagatorOptions<ActionListType, AbortListType>(
      tContext, mContext);
  options.pathLimit = 10000;

  auto propagator = Acts::Propagator<Acts::EigenStepper<>,
                                     Acts::Experimental::DetectorNavigator>(
      stepper, navigator, Acts::getDefaultLogger("Propagator", Acts::Logging::Level::VERBOSE));

  Acts::Vector4 pos(0., 4752, 0., 0.);
  Acts::Vector3 mom(+2, +2, 0);


  Acts::CurvilinearTrackParameters start(pos, mom, +1_e/ mom.norm());
  // propagate to the sector and collect the straw surfaces of the multi layers
  const auto& presult = propagator.propagate(start, options).value();
  auto collector_result = presult.template get<StrawsCollector::result_type>();

  auto end = std::chrono::high_resolution_clock::now();

  auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
  std::cout<<"phi ="<< Acts::VectorHelpers::phi(mom)<<std::endl;
  std::cout<<mom.norm()<<std::endl;

  std::cout<<"collected straw surfaces="<<collector_result.collected.size()<<std::endl;
  std::cout<<"Total CPu time in seconds: "<<elapsed.count()* 1e-9<<std::endl;

  std::cout<<"sector angle ="<<detector_sector->volumes.front()->volumeBounds().values()[3];
  auto sector_angle = detector_sector->volumes.front()->volumeBounds().values()[3];
  auto start_angle = 0.5 * (M_PI - sector_angle);

  std::size_t nTracks = 1000;



}
