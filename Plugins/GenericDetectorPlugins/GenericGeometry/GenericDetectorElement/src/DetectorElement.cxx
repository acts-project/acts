///////////////////////////////////////////////////////////////////
// DetectorElement.cxx, ATS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

// Agd module
#include "GenericDetectorElement/DetectorElement.h"
// Geometry module
#include "Surfaces/PlanarBounds.h"
#include "Surfaces/PlaneSurface.h"
#include "Surfaces/DiscBounds.h"
#include "Surfaces/DiscSurface.h"

/** Constructor for single sided detector element - bound to a Plane Suface */
Agd::DetectorElement::DetectorElement(const Identifier& identifier,
                                      std::shared_ptr<Ats::Transform3D> transform, 
                                      std::shared_ptr<const Ats::PlanarBounds> pBounds,
                                      double thickness) :
    Ats::DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(pBounds.get()),
    m_elementSurface(new Ats::PlaneSurface(*this)),
    m_elementThickness(thickness),
    m_elementCenter(),                    
    m_elementNormal(),                                  
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(pBounds),
    m_elementDiscBounds(nullptr)
{}

/** Constructor for single sided detector element - bound to a Plane Suface */
Agd::DetectorElement::DetectorElement(const Identifier& identifier,
                                   std::shared_ptr<Ats::Transform3D> transform, 
                                   std::shared_ptr<const Ats::DiscBounds> dBounds,
                                   double thickness) :
    Ats::DetectorElementBase(),                                  
    m_elementIdentifier(identifier),
    m_elementTransform(transform),
    m_elementBounds(dBounds.get()),
    m_elementSurface(new Ats::DiscSurface(*this)),
    m_elementThickness(thickness),
    m_elementCenter(),                    
    m_elementNormal(),                                  
    m_elementSurfaces({m_elementSurface}),
    m_elementPlanarBounds(nullptr),
    m_elementDiscBounds(dBounds)
{}

/**  Destructor */
Agd::DetectorElement::~DetectorElement()
{}
