// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"
#include "ACTS/Utilities/Units.hpp"

Acts::DD4hepDetElement::DD4hepDetElement(
    const DD4hep::Geometry::DetElement detElement,
    const std::string&                 axes,
    double                             scalor)
  : Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),
                              detElement.worldTransformation(),
                              detElement.placement().ptr(),
                              axes,
                              scalor)
  , m_detElement(std::move(detElement))

{
  // access the segmentation
  if (m_detElement.volume().isSensitive()) {
    DD4hep::Geometry::SensitiveDetector sensDet(
        m_detElement.volume().sensitiveDetector());
    sensDet.verifyObject();
    m_segmentation = sensDet.readout().segmentation();
  } else
    throw "Detector element is not declared sensitive, can not "
          "access segmentation";
  // @todo when sensitive detector is component
}
