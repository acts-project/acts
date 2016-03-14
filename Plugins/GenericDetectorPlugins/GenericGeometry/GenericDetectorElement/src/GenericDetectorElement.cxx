///////////////////////////////////////////////////////////////////
// GenericDetectorElement.cxx, ATS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

// Ats module
#include "GenericDetectorElement/GenericDetectorElement.h"
// Geometry module
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/DiscBounds.h"
#include "Surfaces/DiscSurface.h"

/** Constructor for single sided detector element - bound to a Plane Suface */
Ats::GenericDetectorElement::GenericDetectorElement(const Identifier& identifier,
                                      std::shared_ptr<Ats::Transform3D> transform, 
                                      std::shared_ptr<const Ats::PlanarBounds> pBounds,
                                      double thickness,
                                      std::shared_ptr<const Ats::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(pBounds.get()),
    m_elementSurface(new Ats::PlaneSurface(*this)),
    m_elementThickness(thickness),
    m_elementCenter(transform->translation()),                    
    m_elementNormal(transform->rotation().col(2)),                                  
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(pBounds),
    m_elementDiscBounds(nullptr)
{
    m_elementSurface->setSurfaceMaterial(material);    
}

/** Constructor for single sided detector element - bound to a Disc Suface */
Ats::GenericDetectorElement::GenericDetectorElement(const Identifier& identifier,
                                   std::shared_ptr<Ats::Transform3D> transform, 
                                   std::shared_ptr<const Ats::DiscBounds> dBounds,
                                   double thickness,
                                   std::shared_ptr<const Ats::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(dBounds.get()),
    m_elementSurface(new Ats::DiscSurface(*this)),
    m_elementThickness(thickness),
    m_elementCenter(transform->translation()),                    
    m_elementNormal(transform->rotation().col(2)),                               
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(nullptr),
    m_elementDiscBounds(dBounds)
{
    m_elementSurface->setSurfaceMaterial(material);
}

/**  Destructor */
Ats::GenericDetectorElement::~GenericDetectorElement()
{}
