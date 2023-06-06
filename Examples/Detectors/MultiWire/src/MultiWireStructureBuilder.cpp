// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MultiWire/MultiWireStructureBuilder.hpp"

#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"
#include "ActsExamples/MultiWire/MultiWireHelper.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

ActsExamples::MultiWireStructureBuilder::MultiWireStructureBuilder(
    const ActsExamples::MultiWireStructureBuilder::Config& config) {
  mCfg = config;

  // check if the surfaces are set
  if (mCfg.strawSurfaces.empty()) {
    std::cout << ("MultiWireStructureBuilder: No surfaces provided- They will "
                  "be set by gdml file")
              << std::endl;
    mCfg.strawSurfaces = ActsExamples::MultiWireHelper::getStrawSurfaces(
        mCfg.sensitiveNames, mCfg.passiveNames);
  }

  multiWireBounds =
      ActsExamples::MultiWireHelper::getMultiWireBounds(mCfg.strawSurfaces);

  // check if the binning is set- if not it will be set using information from
  // the surfaces
  if (mCfg.lbinning.size() < 2u) {
    std::cout << "The binning of the layer structure is not set or it is 1 "
                 "dimensional-it will be set by the standard binning using "
                 "surfaces' information"
              << std::endl;
    mCfg.lbinning = ActsExamples::MultiWireHelper::layerBinning(
        mCfg.strawSurfaces, multiWireBounds);
  }
}

Acts::Experimental::DetectorComponent
ActsExamples::MultiWireStructureBuilder::construct(
    Acts::Experimental::RootDetectorVolumes& roots,
    const Acts::GeometryContext& tContext) {
  auto radius = mCfg.strawSurfaces.front()->bounds().values()[0];

  // Configure the internal structure builder for the internal structure
  Acts::Experimental::LayerStructureBuilder::Config lConfig;
  lConfig.surfaces = [this]() -> std::vector<std::shared_ptr<Acts::Surface>> {
    return mCfg.strawSurfaces;
  };
  lConfig.binnings = mCfg.lbinning;

  // Configure the external structure builder for the external structure
  Acts::ActsScalar hx = mCfg.strawSurfaces.front()->bounds().values()[1];
  Acts::ActsScalar hy =
      abs(multiWireBounds[1].second - multiWireBounds[1].first + 2 * radius);
  Acts::ActsScalar hz =
      abs(multiWireBounds[2].second - multiWireBounds[3].first + 2 * radius);

  Acts::Extent cuboidExtent;
  cuboidExtent.set(Acts::binX, -hx, hx);
  cuboidExtent.set(Acts::binY, multiWireBounds[1].first - radius,
                   multiWireBounds[1].second + radius);
  cuboidExtent.set(Acts::binZ, multiWireBounds[2].first - radius,
                   multiWireBounds[2].second + radius);

  Acts::Experimental::VolumeStructureBuilder::Config vsConfig;
  vsConfig.boundsType = Acts::VolumeBounds::eCuboid;
  vsConfig.extent = cuboidExtent;

  Acts::Experimental::DetectorVolumeBuilder::Config dvConfig;
  dvConfig.auxilliary = "*** Multi Wire Volume ***";
  dvConfig.name = mCfg.name;
  dvConfig.internalsBuilder =
      ActsExamples::MultiWireHelper::internalLayerBuilder(lConfig);
  dvConfig.externalsBuilder =
      ActsExamples::MultiWireHelper::externalVolumeBuilder(vsConfig);
  auto dvBuilder = std::make_shared<Acts::Experimental::DetectorVolumeBuilder>(
      dvConfig,
      Acts::getDefaultLogger("DetectorVolumeBuilder", Acts::Logging::VERBOSE));
  // Acts::Experimental::RootDetectorVolumes roots;
  auto dvComponent = dvBuilder->construct(roots, tContext);

  return dvComponent;
}
