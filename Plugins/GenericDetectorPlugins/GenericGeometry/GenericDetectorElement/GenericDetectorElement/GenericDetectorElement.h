///////////////////////////////////////////////////////////////////
// GenericDetectorElement.h, ACTS project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#ifndef AGD_GENERICDETECTORELEMENT_GENERICDETECTORELEMENT
#define AGD_GENERICDETECTORELEMENT_GENERICDETECTORELEMENT 1

// Algebra and Identifier
#include "CoreIdentifier/Identifier.h"
#include "Algebra/AlgebraDefinitions.h"
// Geometry module
#include "DetectorElementBase/DetectorElementBase.h"

namespace Acts {

    class Surface;
    class PlanarBounds;
    class DiscBounds;   
    class SurfaceMaterial; 
    
    /** @class GenericDetectorElement
     
     This is a lightweight type of detector element, 
     it simply implements the base class.
     
     @author Andreas.Salzburger@cern.ch
     
     */
    
    class GenericDetectorElement : public DetectorElementBase {
        
      public:
        
         /** Constructor for single sided detector element - bound to a Plane Surface  */
         GenericDetectorElement(const Identifier& identifier,
                                std::shared_ptr<Transform3D> transform, 
                                std::shared_ptr<const PlanarBounds> pBounds,
                                double thickness,
                                std::shared_ptr<const SurfaceMaterial> material = nullptr);

        /** Constructor for single sided detector element - bound to a Plane Surface / Disc Surface */
        GenericDetectorElement(const Identifier& identifier,
                               std::shared_ptr<Transform3D> transform, 
                               std::shared_ptr<const DiscBounds> dBounds,
                               double thickness,
                               std::shared_ptr<const SurfaceMaterial> material = nullptr);
                                        
        /**  Destructor */
        ~GenericDetectorElement();
        
        /** Identifier */
        Identifier identify() const override;
        
        /** Return local to global transform*/
        const Transform3D& transform() const override;
        
        /** Return local to global transform associated with this identifier*/
        const Transform3D& transform(const Identifier& id) const override;
        
        /** Return surface associated with this detector element*/
        const Surface& surface () const override;
        
        /** Return surface associated with this identifier, which should come from the */
        const Surface& surface (const Identifier& identifier) const override;
        
        /** Returns the full list of all detection surfaces associated to this detector element */
        const std::vector< std::shared_ptr<const Surface> >& surfaces() const override;
        
        /** Return the boundaries of the element*/
        const SurfaceBounds& bounds() const override;
        
        /** Return the boundaries of the surface associated with this identifier */
        const SurfaceBounds& bounds(const Identifier& id) const override;
        
        /** Return the center of the element*/
        const Vector3D& center() const override;
        
        /** Return the center of the surface associated with this identifier
         In the case of silicon it returns the same as center()*/
        const Vector3D& center(const Identifier& identifier) const override;
        
        /** Return the normal of the element*/
        const Vector3D& normal() const override;
        
        /** Return the normal of the surface associated with this identifier
         In the case of silicon it returns the same as normal()*/
        const Vector3D& normal(const Identifier& id) const override;
         
        /** The maximal thickness of the detector element outside the surface dimension */
        double thickness() const override;
         
       private: 
         // the element representation  
         Identifier                                         m_elementIdentifier;
         std::shared_ptr<Transform3D>                       m_elementTransform;
         const SurfaceBounds*                               m_elementBounds;
                                                            
         std::shared_ptr<const Surface>                     m_elementSurface;
         double                                             m_elementThickness;
                                                            
         // the cache                                       
         Vector3D                                           m_elementCenter;
         Vector3D                                           m_elementNormal;
         std::vector< std::shared_ptr<const Surface> >      m_elementSurfaces;
                                                          
         std::shared_ptr<const PlanarBounds>                m_elementPlanarBounds;
         std::shared_ptr<const DiscBounds>                  m_elementDiscBounds;   


        
    };
    
    inline Identifier GenericDetectorElement::identify() const { return m_elementIdentifier; }
    
    inline const Transform3D& GenericDetectorElement::transform() const { return (*m_elementTransform.get()); }
   
    inline const Transform3D& GenericDetectorElement::transform(const Identifier& ) const { return transform(); }
    
    inline const Surface& GenericDetectorElement::surface () const { return (*m_elementSurface.get()); }
   
    inline const Surface& GenericDetectorElement::surface (const Identifier& ) const { return surface(); }
   
    inline const std::vector< std::shared_ptr<const Surface> >& GenericDetectorElement::surfaces() const { return m_elementSurfaces; }
    
    inline const SurfaceBounds& GenericDetectorElement::bounds() const { return (*m_elementBounds); }
   
    inline const SurfaceBounds& GenericDetectorElement::bounds(const Identifier& ) const { return bounds(); }
    
    inline const Vector3D& GenericDetectorElement::center() const { return m_elementCenter; }
   
    inline const Vector3D& GenericDetectorElement::center(const Identifier& ) const { return center(); }
   
    inline const Vector3D& GenericDetectorElement::normal() const { return m_elementNormal; }
   
    inline const Vector3D& GenericDetectorElement::normal(const Identifier&) const { return normal(); }

    inline double GenericDetectorElement::thickness() const { return m_elementThickness; }


    
}//end of ns

#endif
