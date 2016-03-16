///////////////////////////////////////////////////////////////////
// TGeoDetectorElement.h, ACTS project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT
#define ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT 1

// Geometry module
#include "DetectorElementBase/DetectorElementBase.h"
//Root
#include "TGeoManager.h"

namespace Acts {
    
    /** @class TGeoDetectorElement
     
     DetectorElement plugin for ROOT TGeo shapes.
     
     @author julia.hrdinka@cern.ch
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class TGeoDetectorElement : public Acts::DetectorElementBase {
    
    public:
        /** Constructor  */
        TGeoDetectorElement(const Identifier& identifier,
                            TGeoShape* tGeoDetElement, std::shared_ptr<const Acts::Transform3D> motherTransform = nullptr);
        
        /**  Destructor */
        ~TGeoDetectorElement();
        
        /** Identifier */
        Identifier identify() const override;
        
        /**Return local to global transform*/
        const Acts::Transform3D& transform() const override;
        
        /**Return local to global transform associated with this identifier*/
        const Acts::Transform3D& transform(const Identifier& id) const override;
        
        /**Return surface associated with this detector element*/
        const Acts::Surface& surface () const override;
        
        /**Return surface associated with this identifier, which should come from the */
        const Acts::Surface& surface (const Identifier& identifier) const override;
        
        /** Returns the full list of all detection surfaces associated to this detector element */
        const std::vector< std::shared_ptr<const Acts::Surface> >& surfaces() const override;
        
        /**Return the boundaries of the element*/
        const Acts::SurfaceBounds& bounds() const override;
        
        /**Return the boundaries of the surface associated with this identifier */
        const Acts::SurfaceBounds& bounds(const Identifier& id) const override;
        
        /**Return the center of the element*/
        const Acts::Vector3D& center() const override;
        
        /**Return the center of the surface associated with this identifier
         In the case of silicon it returns the same as center()*/
        const Acts::Vector3D& center(const Identifier& identifier) const override;
        
        /**Return the normal of the element*/
        const Acts::Vector3D& normal() const override;
        
        /**Return the normal of the surface associated with this identifier
         In the case of silicon it returns the same as normal()*/
        const Acts::Vector3D& normal(const Identifier& id) const override;

        
    private:
        
        /**DD4hep detector element*/
        TGeoShape*                                          m_detElement;
        /**Transformation of the detector element*/
        std::shared_ptr<Acts::Transform3D>                   m_transform;
        /**Center position of the detector element*/
        mutable std::shared_ptr<const Acts::Vector3D>        m_center;
        /**Normal vector to the detector element*/
        mutable std::shared_ptr<const Acts::Vector3D>        m_normal;
        /**Identifier of the detector element*/
        const Identifier                                    m_identifier;
        /**Boundaries of the detector element*/
        std::shared_ptr<const Acts::SurfaceBounds>           m_bounds;
        /**Corresponding Surface*/
        std::shared_ptr<const Acts::Surface>                 m_surface;
        /**possible contained surfaces*/
        std::vector<std::shared_ptr<const Acts::Surface>>    m_surfaces;
    };
}

inline Identifier Acts::TGeoDetectorElement::identify() const {return m_identifier;}

inline const Acts::Transform3D& Acts::TGeoDetectorElement::transform() const {return (*m_transform);}

inline const Acts::Transform3D& Acts::TGeoDetectorElement::transform(const Identifier&) const {return (*m_transform);}

inline const Acts::Surface& Acts::TGeoDetectorElement::surface() const {return (*m_surface);}

inline const Acts::Surface& Acts::TGeoDetectorElement::surface(const Identifier&) const {return (*m_surface);}

inline const std::vector< std::shared_ptr<const Acts::Surface> >& Acts::TGeoDetectorElement::surfaces() const {return (m_surfaces);}

inline const Acts::SurfaceBounds& Acts::TGeoDetectorElement::bounds() const{return (*m_bounds);}

inline const Acts::SurfaceBounds& Acts::TGeoDetectorElement::bounds(const Identifier&) const {return (*m_bounds);}


#endif
