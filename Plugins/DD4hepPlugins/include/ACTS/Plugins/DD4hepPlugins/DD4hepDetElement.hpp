// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepDetElement.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
#define ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H 1

//TGeoDetectorElement
#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
//DD4hep
#include "DD4hep/Detector.h"

namespace Acts {
    
    /** @class DD4hepDetElement
     
     DetectorElement plugin for DD4hep detector elements. DD4hep is based on TGeo shapes, therefore the DD4hepDetElement inherits from TGeoDetectorElement. The full geometrical information is provided by the TGeoDetectorElement. Only the TGeoShape and the DD4hepIdentifier need to be provided.
     
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class DD4hepDetElement : public TGeoDetectorElement {
        
    public:
        
        /** Constructor */
        DD4hepDetElement(const DD4hep::Geometry::DetElement detElement, const DD4hep::Geometry::Segmentation segmentation, std::shared_ptr<const Acts::Transform3D> motherTransform=nullptr);
        /** Desctructor */
        virtual ~DD4hepDetElement();
        
    private:
        /** DD4hep detector element */
        DD4hep::Geometry::DetElement            m_detElement;
        /** DD4hep segmentation */
        DD4hep::Geometry::Segmentation          m_segmentation;
        
    };
}


#endif //ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
