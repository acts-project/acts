// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <DD4hep/CartesianGridXY.h>

Acts::DD4hepDetectorElement::DD4hepDetectorElement(
    const dd4hep::DetElement detElement, const std::string& axes, double scalor,
    bool /*isDisc*/, std::shared_ptr<const Acts::ISurfaceMaterial> material,
    std::shared_ptr<const Acts::DigitizationModule> /*digitizationModule*/)
    : Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),
                                *(detElement.placement().ptr()),
                                detElement.nominal().worldTransformation(),
                                axes, scalor, std::move(material)),
      m_detElement(detElement) {}
