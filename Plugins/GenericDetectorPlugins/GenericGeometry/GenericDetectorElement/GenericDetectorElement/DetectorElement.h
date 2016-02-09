///////////////////////////////////////////////////////////////////
// DetectorElement.h, ATS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#ifndef AGD_GENERICDETECTORELEMENT_DETECTORELEMENT
#define AGD_GENERICDETECTORELEMENT_DETECTORELEMENT 1

// Algebra and Identifier
#include "Identifier/Identifier.h"
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "DetectorElementBase/DetectorElementBase.h"

namespace Ats {
    class Surface;
    class PlanarBounds;
    class DiscBounds;    
}

// namespace for A generic detector
namespace Agd
{
    
    /** @class DetectorElement
     
     This is easiest type of detector element, simple implements the base class.
    
      @TODO implement double sided
     
     @author Andreas.Salzburger@cern.ch
     
     */
    
    class DetectorElement : public Ats::DetectorElementBase {
        
      public:
        
         /** Constructor for single sided detector element - bound to a Plane Surface  */
         DetectorElement(const Identifier& identifier,
                        std::shared_ptr<Ats::Transform3D> transform, 
                        std::shared_ptr<const Ats::PlanarBounds> pBounds,
                        double thickness);

        /** Constructor for single sided detector element - bound to a Plane Surface / Disc Surface */
        DetectorElement(const Identifier& identifier,
                        std::shared_ptr<Ats::Transform3D> transform, 
                        std::shared_ptr<const Ats::DiscBounds> dBounds,
                        double thickness);
                                        
        /**  Destructor */
         ~DetectorElement();
        
        /** Identifier */
         Identifier identify() const override;
        
        /**Return local to global transform*/
         const Ats::Transform3D& transform() const override;
        
        /**Return local to global transform associated with this identifier*/
         const Ats::Transform3D& transform(const Identifier& id) const override;
        
        /**Return surface associated with this detector element*/
         const Ats::Surface& surface () const override;
        
        /**Return surface associated with this identifier, which should come from the */
         const Ats::Surface& surface (const Identifier& identifier) const override;
        
        /** Returns the full list of all detection surfaces associated to this detector element */
         const std::vector< std::shared_ptr<const Ats::Surface> >& surfaces() const override;
        
        /**Return the boundaries of the element*/
         const Ats::SurfaceBounds& bounds() const override;
        
        /**Return the boundaries of the surface associated with this identifier */
         const Ats::SurfaceBounds& bounds(const Identifier& id) const override;
        
        /**Return the center of the element*/
         const Ats::Vector3D& center() const override;
        
        /**Return the center of the surface associated with this identifier
         In the case of silicon it returns the same as center()*/
         const Ats::Vector3D& center(const Identifier& identifier) const override;
        
        /**Return the normal of the element*/
         const Ats::Vector3D& normal() const override;
        
        /**Return the normal of the surface associated with this identifier
         In the case of silicon it returns the same as normal()*/
         const Ats::Vector3D& normal(const Identifier& id) const override;
         
       private: 
         // the element representation  
         Identifier                                          m_elementIdentifier;
         std::shared_ptr<Ats::Transform3D>                   m_elementTransform;
         const Ats::SurfaceBounds*                           m_elementBounds;

         std::shared_ptr<const Ats::Surface>                 m_elementSurface;
         double                                              m_elementThickness;
                                                        
         // the cache                                   
         Ats::Vector3D                                       m_elementCenter;
         Ats::Vector3D                                       m_elementNormal;
         std::vector< std::shared_ptr<const Ats::Surface> >  m_elementSurfaces;
          
         std::shared_ptr<const Ats::PlanarBounds>            m_elementPlanarBounds;
         std::shared_ptr<const Ats::DiscBounds>              m_elementDiscBounds;   
         
           
        
    };
    
    
    inline Identifier DetectorElement::identify() const { return m_elementIdentifier; }
    
    inline const Ats::Transform3D& DetectorElement::transform() const { return (*m_elementTransform.get()); }
   
    inline const Ats::Transform3D& DetectorElement::transform(const Identifier& ) const { return transform(); }
    
    inline const Ats::Surface& DetectorElement::surface () const { return (*m_elementSurface.get()); }
   
    inline const Ats::Surface& DetectorElement::surface (const Identifier& ) const { return surface(); }
   
    inline const std::vector< std::shared_ptr<const Ats::Surface> >& DetectorElement::surfaces() const { return m_elementSurfaces; }
    
    inline const Ats::SurfaceBounds& DetectorElement::bounds() const { return (*m_elementBounds); }
   
    inline const Ats::SurfaceBounds& DetectorElement::bounds(const Identifier& ) const { return bounds(); }
    
    inline const Ats::Vector3D& DetectorElement::center() const { return m_elementCenter; }
   
    inline const Ats::Vector3D& DetectorElement::center(const Identifier& ) const { return center(); }
   
    inline const Ats::Vector3D& DetectorElement::normal() const { return m_elementNormal; }
   
    inline const Ats::Vector3D& DetectorElement::normal(const Identifier&) const { return normal(); }
    
    
}//end of ns

#endif
