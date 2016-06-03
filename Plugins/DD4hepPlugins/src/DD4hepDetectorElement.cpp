// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"

Acts::DD4hepDetElement::DD4hepDetElement(const DD4hep::Geometry::DetElement detElement, const DD4hep::Geometry::Segmentation segmentation, std::shared_ptr<const Acts::Transform3D> motherTransform) :
Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),detElement.placement().ptr(),motherTransform),
m_detElement(std::move(detElement)),
m_segmentation(std::move(segmentation))
{}

Acts::DD4hepDetElement::~DD4hepDetElement()
{}