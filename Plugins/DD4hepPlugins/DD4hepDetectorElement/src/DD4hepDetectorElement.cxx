//
//  DD4hepDetectorElement.cxx
//  
//
//  Created by Julia Hrdinka on 26/02/16.
//
//

#include "DD4hepDetectorElement/DD4hepDetectorElement.h"

Ats::DD4hepDetectorElement::DD4hepDetectorElement(const DD4hep::Geometry::DetElement& detElement, const DD4hep::Geometry::Segmentation& segmentation, std::shared_ptr<const Ats::Transform3D> motherTransform) :
Ats::TGeoDetectorElement(Identifier(detElement.volumeID()),detElement.placement().ptr()->GetVolume()->GetShape(),motherTransform),
m_detElement(detElement),
m_segmentation(segmentation)
{}

Ats::DD4hepDetectorElement::~DD4hepDetectorElement()
{}