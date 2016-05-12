#include "ACTS/Plugins/DD4hepPlugins/DD4hepDetElement.hpp"

Acts::DD4hepDetElement::DD4hepDetElement(const DD4hep::Geometry::DetElement detElement, const DD4hep::Geometry::Segmentation segmentation, std::shared_ptr<const Acts::Transform3D> motherTransform) :
Acts::TGeoDetectorElement(Identifier(detElement.volumeID()),detElement.placement().ptr(),motherTransform),
m_detElement(std::move(detElement)),
m_segmentation(std::move(segmentation))
{}

Acts::DD4hepDetElement::~DD4hepDetElement()
{}