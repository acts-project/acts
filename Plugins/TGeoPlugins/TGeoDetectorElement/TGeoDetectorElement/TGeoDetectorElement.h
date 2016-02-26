///////////////////////////////////////////////////////////////////
// TGeoDetectorElement.h, ATS project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ATS_TGEODETECTORELEMENT_TGEODETECTORELEMENT
#define ATS_TGEODETECTORELEMENT_TGEODETECTORELEMENT 1

// Geometry module
#include "DetectorElementBase/DetectorElementBase.h"
//Root
#include "TGeoManager.h"

namespace Ats {
    
    /** @class TGeoDetectorElement
     
     DetectorElement plugin for ROOT TGeo shapes.
     
     @author julia.hrdinka@cern.ch
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class TGeoDetectorElement : public Ats::DetectorElementBase {
    
    public:
        /** Constructor  */
        TGeoDetectorElement(const Identifier& identifier,
                            TGeoShape* tGeoDetElement, std::shared_ptr<const Ats::Transform3D> motherTransform = nullptr);
        
        /**  Destructor */
        virtual ~TGeoDetectorElement();
        
        /** Identifier */
        virtual Identifier identify() const override;
        
        /**Return local to global transform*/
        virtual const Ats::Transform3D& transform() const override;
        
        /**Return local to global transform associated with this identifier*/
        virtual const Ats::Transform3D& transform(const Identifier& id) const override;
        
        /**Return surface associated with this detector element*/
        virtual const Ats::Surface& surface () const override;
        
        /**Return surface associated with this identifier, which should come from the */
        virtual const Ats::Surface& surface (const Identifier& identifier) const override;
        
        /** Returns the full list of all detection surfaces associated to this detector element */
        virtual const std::vector< std::shared_ptr<const Ats::Surface> >& surfaces() const override;
        
        /**Return the boundaries of the element*/
        virtual const Ats::SurfaceBounds& bounds() const override;
        
        /**Return the boundaries of the surface associated with this identifier */
        virtual const Ats::SurfaceBounds& bounds(const Identifier& id) const override;
        
        /**Return the center of the element*/
        virtual const Ats::Vector3D& center() const override;
        
        /**Return the center of the surface associated with this identifier
         In the case of silicon it returns the same as center()*/
        virtual const Ats::Vector3D& center(const Identifier& identifier) const override;
        
        /**Return the normal of the element*/
        virtual const Ats::Vector3D& normal() const override;
        
        /**Return the normal of the surface associated with this identifier
         In the case of silicon it returns the same as normal()*/
        virtual const Ats::Vector3D& normal(const Identifier& id) const override;

        
    private:
        
        /**DD4hep detector element*/
        TGeoShape*                                          m_detElement;
        /**Transformation of the detector element*/
        std::shared_ptr<Ats::Transform3D>                   m_transform;
        /**Center position of the detector element*/
        mutable std::shared_ptr<const Ats::Vector3D>        m_center;
        /**Normal vector to the detector element*/
        mutable std::shared_ptr<const Ats::Vector3D>        m_normal;
        /**Identifier of the detector element*/
        const Identifier                                    m_identifier;
        /**Boundaries of the detector element*/
        std::shared_ptr<const Ats::SurfaceBounds>           m_bounds;
        /**Corresponding Surface*/
        std::shared_ptr<const Ats::Surface>                 m_surface;
        /**possible contained surfaces*/
        std::vector<std::shared_ptr<const Ats::Surface>>    m_surfaces;
    };
}

inline Identifier Ats::TGeoDetectorElement::identify() const {return m_identifier;}

inline const Ats::Transform3D& Ats::TGeoDetectorElement::transform() const {return (*m_transform);}

inline const Ats::Transform3D& Ats::TGeoDetectorElement::transform(const Identifier&) const {return (*m_transform);}

inline const Ats::Surface& Ats::TGeoDetectorElement::surface() const {return (*m_surface);}

inline const Ats::Surface& Ats::TGeoDetectorElement::surface(const Identifier&) const {return (*m_surface);}

inline const std::vector< std::shared_ptr<const Ats::Surface> >& Ats::TGeoDetectorElement::surfaces() const {return (m_surfaces);}

inline const Ats::SurfaceBounds& Ats::TGeoDetectorElement::bounds() const{return (*m_bounds);}

inline const Ats::SurfaceBounds& Ats::TGeoDetectorElement::bounds(const Identifier&) const {return (*m_bounds);}


#endif
