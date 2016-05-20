// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DD4hepGeometryHelper.hpp"
// Geometry Module
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
// Root TGeo
#include "TGeoManager.h"

std::shared_ptr<Acts::Transform3D> Acts::DD4hepGeometryHelper::extractTransform(DD4hep::Geometry::DetElement& detElement)
{
    //get the placement and orientation in respect to its mother
    const Double_t* rotation    = (detElement.placement().ptr()->GetMatrix()->GetRotationMatrix());
    const Double_t* translation = (detElement.placement().ptr()->GetMatrix()->GetTranslation());
    auto transform = std::make_shared<Acts::Transform3D>(Acts::Vector3D(rotation[0],rotation[3],rotation[6]),Acts::Vector3D(rotation[1],rotation[4],rotation[7]),Acts::Vector3D(rotation[2],rotation[5],rotation[8]), Acts::Vector3D(translation[0],translation[1], translation[2]));
    return (transform);
}

std::shared_ptr<const Acts::VolumeBounds> Acts::DD4hepGeometryHelper::extractVolumeBounds(DD4hep::Geometry::DetElement& detElement)
{
    TGeoShape* geoShape = detElement.placement().ptr()->GetVolume()->GetShape();
    TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoShape);
    if (!tube) throw "Volume has wrong shape - needs to be TGeoConeSeg!";
    auto cylinderBounds = std::make_shared<const Acts::CylinderVolumeBounds>(tube->GetRmin1(),tube->GetRmax1(),tube->GetDz());
    return cylinderBounds;
}
