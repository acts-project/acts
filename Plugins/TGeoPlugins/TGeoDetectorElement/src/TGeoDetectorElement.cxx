#include "TGeoDetectorElement/TGeoDetectorElement.h"
//ATS
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/TrapezoidBounds.h"
#include "Surfaces/RectangleBounds.h"

#include "TGeoBBox.h"
#include "TGeoTrd2.h"

Acts::TGeoDetectorElement::TGeoDetectorElement(const Identifier& identifier,
                TGeoShape* tGeoDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform) :
Acts::DetectorElementBase(),
m_detElement(tGeoDetElement),
m_identifier(identifier)
{
    //get the placement and orientation in respect to its mother
    const Double_t* rotation    = (m_detElement->GetTransform()->GetRotationMatrix());
    const Double_t* translation = (m_detElement->GetTransform()->GetTranslation());
    auto transform = std::make_shared<Acts::Transform3D>(Acts::Vector3D(rotation[0],rotation[3],rotation[6]),Acts::Vector3D(rotation[1],rotation[4],rotation[7]),Acts::Vector3D(rotation[2],rotation[5],rotation[8]), Acts::Vector3D(translation[0],translation[1], translation[2]));
    //now calculate the global transformation
    if (motherTransform) (*transform) = (*motherTransform)*(*transform);
    m_transform  = transform;
    //currently only surfaces with rectangular or trapezoidal shape implemented
    TGeoBBox* box       = dynamic_cast<TGeoBBox*>(tGeoDetElement);
    TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(tGeoDetElement);
    if (trapezoid) {
        //extract the surface bounds
        auto trapezoidBounds = std::make_shared<const Acts::TrapezoidBounds>(trapezoid->GetDx1(),trapezoid->GetDx2(),trapezoid->GetDz());
        m_bounds = trapezoidBounds;
        m_surface = std::make_shared<const Acts::PlaneSurface>(transform,trapezoidBounds);
    }
    else {
        //extract the surface bounds
        auto rectangleBounds = std::make_shared<const Acts::RectangleBounds>(box->GetDX(),box->GetDY());
        m_bounds = rectangleBounds;
        m_surface = std::make_shared<const Acts::PlaneSurface>(transform,rectangleBounds);
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

