///////////////////////////////////////////////////////////////////
// GenericDetectorElement.cxx, ACTS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

// Acts module
#include "GenericDetectorElement/GenericDetectorElement.h"
// Geometry module
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/DiscBounds.h"
#include "Surfaces/DiscSurface.h"

/** Constructor for single sided detector element - bound to a Plane Suface */
Acts::GenericDetectorElement::GenericDetectorElement(const Identifier& identifier,
                                      std::shared_ptr<Acts::Transform3D> transform, 
                                      std::shared_ptr<const Acts::PlanarBounds> pBounds,
                                      double thickness,
                                      std::shared_ptr<const Acts::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(pBounds.get()),
    m_elementSurface(new Acts::PlaneSurface(*this)),
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
Acts::GenericDetectorElement::GenericDetectorElement(const Identifier& identifier,
                                   std::shared_ptr<Acts::Transform3D> transform, 
                                   std::shared_ptr<const Acts::DiscBounds> dBounds,
                                   double thickness,
                                   std::shared_ptr<const Acts::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(dBounds.get()),
    m_elementSurface(new Acts::DiscSurface(*this)),
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
Acts::GenericDetectorElement::~GenericDetectorElement()
{}
