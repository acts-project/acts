///////////////////////////////////////////////////////////////////
// GenericDetectorElement.cxx, ACTS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

// Acts module
#include "ACTS/Examples/GenericDetectorExample/GenericDetectorElement.hpp"
// Geometry module
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/DiscSurface.hpp"

/** Constructor for single sided detector element - bound to a Plane Suface */
Acts::GenericDetectorElement::GenericDetectorElement(const Identifier identifier,
                                      std::shared_ptr<Acts::Transform3D> transform, 
                                      std::shared_ptr<const Acts::PlanarBounds> pBounds,
                                      double thickness,
                                      std::shared_ptr<const Acts::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(std::move(identifier)),
    m_elementTransform(std::move(transform)),
    m_elementBounds(pBounds.get()),
    m_elementSurface(new Acts::PlaneSurface(*this)),
    m_elementThickness(std::move(thickness)),
    m_elementCenter(m_elementTransform->translation()),
    m_elementNormal(m_elementTransform->rotation().col(2)),
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(std::move(pBounds)),
    m_elementDiscBounds(nullptr)
{
    m_elementSurface->setSurfaceMaterial(material);
}

/** Constructor for single sided detector element - bound to a Disc Suface */
Acts::GenericDetectorElement::GenericDetectorElement(const Identifier identifier,
                                   std::shared_ptr<Acts::Transform3D> transform, 
                                   std::shared_ptr<const Acts::DiscBounds> dBounds,
                                   double thickness,
                                   std::shared_ptr<const Acts::SurfaceMaterial> material) :
    DetectorElementBase(),                                  
    m_elementIdentifier(std::move(identifier)),
    m_elementTransform(std::move(transform)),
    m_elementBounds(dBounds.get()),
    m_elementSurface(new Acts::DiscSurface(*this)),
    m_elementThickness(std::move(thickness)),
    m_elementCenter(transform->translation()),                    
    m_elementNormal(transform->rotation().col(2)),                               
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(nullptr),
    m_elementDiscBounds(std::move(dBounds))
{
    m_elementSurface->setSurfaceMaterial(material);
}

/**  Destructor */
Acts::GenericDetectorElement::~GenericDetectorElement()
{}

