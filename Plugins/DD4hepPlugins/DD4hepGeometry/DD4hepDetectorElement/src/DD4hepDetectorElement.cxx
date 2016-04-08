#include "DD4hepDetectorElement/DD4hepDetElement.h"

Acts::DD4hepDetElement::DD4hepDetElement(const DD4hep::Geometry::DetElement& detElement, const DD4hep::Geometry::Segmentation& segmentation, std::shared_ptr<const Acts::Transform3D> motherTransform) :
Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),detElement.placement().ptr(),motherTransform),
m_detElement(detElement),
m_segmentation(segmentation)
{}

Acts::DD4hepDetElement::~DD4hepDetElement()
{}