#include "DD4hepGeometryUtils/DD4hepGeometryHelper.h"
// Geometry Module
#include "Volumes/CylinderVolumeBounds.h"
// Root TGeo
#include "TGeoManager.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"

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
    if (!tube) throw GaudiException("Volume has wrong shape - needs to be TGeoConeSeg!", "FATAL", StatusCode::FAILURE);
    auto cylinderBounds = std::make_shared<const Acts::CylinderVolumeBounds>(tube->GetRmin1(),tube->GetRmax1(),tube->GetDz());
    return cylinderBounds;
}