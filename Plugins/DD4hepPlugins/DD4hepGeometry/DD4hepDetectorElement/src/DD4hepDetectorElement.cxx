//
//  DD4hepDetElement.cxx
//  
//
//  Created by Julia Hrdinka on 26/02/16.
//
//

#include "DD4hepDetectorElement/DD4hepDetElement.h"

Add4hep::DD4hepDetElement::DD4hepDetElement(const DD4hep::Geometry::DetElement& detElement, const DD4hep::Geometry::Segmentation& segmentation, std::shared_ptr<const Ats::Transform3D> motherTransform) :
Atgeo::TGeoDetectorElement(Identifier(detElement.volumeID()),detElement.placement().ptr()->GetVolume()->GetShape(),motherTransform),
m_detElement(detElement),
m_segmentation(segmentation)
{}

Add4hep::DD4hepDetElement::~DD4hepDetElement()
{}