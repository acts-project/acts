#include "TGeoDetectorElement/TGeoDetectorElement.h"
// ACTS
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Surfaces/RectangleBounds.h"

#include "Material/Material.h"
#include "Material/MaterialProperties.h"
#include "Material/HomogeneousSurfaceMaterial.h"

#include "TGeoBBox.h"
#include "TGeoTrd2.h"

#include "Algebra/AlgebraDefinitions.h"

#include <iostream>


Acts::TGeoDetectorElement::TGeoDetectorElement(const Identifier& identifier,
                TGeoNode* tGeoDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform) :
  Acts::DetectorElementBase(),
  m_detElement(tGeoDetElement),
  m_identifier(identifier),
  m_thickness(0.)
{
    //get the placement and orientation in respect to its mother
    const Double_t* rotation    = (m_detElement->GetMatrix()->GetRotationMatrix());
    const Double_t* translation = (m_detElement->GetMatrix()->GetTranslation());
    //currently only surfaces with rectangular or trapezoidal shape implemented
    TGeoBBox* box       = dynamic_cast<TGeoBBox*>(m_detElement->GetVolume()->GetShape());
    TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(m_detElement->GetVolume()->GetShape());
    TGeoMaterial* mat = m_detElement->GetVolume()->GetMaterial();
    Material moduleMaterial (mat->GetRadLen(),mat->GetIntLen(),mat->GetA(),mat->GetZ(),mat->GetDensity());
    if (trapezoid) {
        // if the shape is TGeoTrd2 y and z axes needs to be exchanged, since in TGei the description is different
        m_transform = std::make_shared<Acts::Transform3D>(Acts::Vector3D(rotation[0],rotation[3],rotation[6]),Acts::Vector3D(rotation[1],rotation[4],rotation[7]),Acts::Vector3D(rotation[2],rotation[5],rotation[8]), Acts::Vector3D(translation[0],translation[1], translation[2]));
        //now calculate the global transformation
        if (motherTransform) (*m_transform) = (*motherTransform)*(*m_transform);
       // Rotation3D rotation = AngleAxis3D(0.5*M_PI, Vector3D::UnitZ());//*AngleAxis3D(0.5*M_PI, Vector3D::UnitX());
        (*m_transform) *= AngleAxis3D(0.5*M_PI, Vector3D::UnitX());
        //extract the surface bounds
        auto trapezoidBounds = std::make_shared<const Acts::TrapezoidBounds>(trapezoid->GetDx1(),trapezoid->GetDx2(),trapezoid->GetDz());
        m_bounds = trapezoidBounds;
        m_surface = std::make_shared<const Acts::PlaneSurface>(*this,m_identifier);
        MaterialProperties moduleMaterialProperties(moduleMaterial,0.5*(trapezoid->GetDy1()+trapezoid->GetDy1()));
        m_surface->setSurfaceMaterial(std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(moduleMaterialProperties)));
    }
    else {
        m_transform = std::make_shared<Acts::Transform3D>(Acts::Vector3D(rotation[0],rotation[3],rotation[6]),Acts::Vector3D(rotation[1],rotation[4],rotation[7]),Acts::Vector3D(rotation[2],rotation[5],rotation[8]), Acts::Vector3D(translation[0],translation[1], translation[2]));
        //now calculate the global transformation
        if (motherTransform) (*m_transform) = (*motherTransform)*(*m_transform);
        //extract the surface bounds
        auto rectangleBounds = std::make_shared<const Acts::RectangleBounds>(box->GetDX(),box->GetDY());
        m_bounds = rectangleBounds;
        m_surface = std::make_shared<const Acts::PlaneSurface>(*this,m_identifier);
        MaterialProperties moduleMaterialProperties(moduleMaterial,box->GetDZ());
        m_surface->setSurfaceMaterial(std::shared_ptr<const SurfaceMaterial>(new HomogeneousSurfaceMaterial(moduleMaterialProperties)));
    }
}

Acts::TGeoDetectorElement::~TGeoDetectorElement()
{}

const Acts::Vector3D& Acts::TGeoDetectorElement::center() const
{
    if (!m_center) m_center = std::make_shared<const Acts::Vector3D>(m_transform->translation());
    return (*m_center);
}


const Acts::Vector3D& Acts::TGeoDetectorElement::center(const Identifier&) const
{
    if (!m_center) m_center = std::make_shared<const Acts::Vector3D>(m_transform->translation());
    return (*m_center);
}


const Acts::Vector3D& Acts::TGeoDetectorElement::normal() const
{
    if (!m_normal) m_normal = std::make_shared<const Acts::Vector3D>(Acts::Vector3D(m_transform->rotation().col(2)));
    return(*m_normal);
}


const Acts::Vector3D& Acts::TGeoDetectorElement::normal(const Identifier&) const
{
    if (!m_normal) m_normal = std::make_shared<const Acts::Vector3D>(Acts::Vector3D(m_transform->rotation().col(2)));
    return(*m_normal);
}

