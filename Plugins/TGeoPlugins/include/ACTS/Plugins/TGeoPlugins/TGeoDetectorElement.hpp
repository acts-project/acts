// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoDetectorElement.h, ACTS project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT
#define ACTS_TGEODETECTORELEMENT_TGEODETECTORELEMENT 1

// Geometry module
#include "ACTS/Detector/DetectorElementBase.hpp"
//Root
#include "TGeoManager.h"
#include<iostream>

namespace Acts {
    
    /** @class TGeoDetectorElement
     
     DetectorElement plugin for ROOT TGeo shapes.
     
     @TODO what if shape conversion failes? add implementation of more than one surface per module, implementing also for other shapes->Cone,ConeSeg,Tube? what if not used with DD4hep?
     
     */
    
    class TGeoDetectorElement : public DetectorElementBase {
    
    public:
        /** Constructor  */
        TGeoDetectorElement(const Identifier& identifier,
                            TGeoNode* tGeoDetElement, std::shared_ptr<const Transform3D> motherTransform = nullptr);
        
        /**  Destructor */
        virtual ~TGeoDetectorElement();
        
        /** Identifier */
        virtual Identifier identify() const override;
        
        /**Return local to global transform associated with this identifier*/
        virtual const Transform3D& transform(const Identifier& identifier = Identifier()) const override;
        
        /**Return surface associated with this identifier, which should come from the */
        virtual const Surface& surface (const Identifier& identifier = Identifier()) const override;
        
        /** Returns the full list of all detection surfaces associated to this detector element */
        virtual const std::vector< std::shared_ptr<const Surface> >& surfaces() const override;
        
        /**Return the boundaries of the surface associated with this identifier */
        virtual const SurfaceBounds& bounds(const Identifier& identifier = Identifier()) const override;
        
        /**Return the center of the surface associated with this identifier
         In the case of silicon it returns the same as center()*/
        virtual const Vector3D& center(const Identifier& identifier = Identifier()) const override;
        
        /**Return the normal of the surface associated with this identifier
         In the case of silicon it returns the same as normal()*/
        virtual const Vector3D& normal(const Identifier& identifier = Identifier()) const override;
        
        /** Returns the thickness of the module */
        virtual double thickness() const override;
        
    private:
        
        /**DD4hep detector element*/
        TGeoNode*                                      m_detElement;
        /**Transformation of the detector element*/
        std::shared_ptr<Transform3D>                   m_transform;
        /**Center position of the detector element*/
        mutable std::shared_ptr<const Vector3D>        m_center;
        /**Normal vector to the detector element*/
        mutable std::shared_ptr<const Vector3D>        m_normal;
        /**Identifier of the detector element*/
        const Identifier                               m_identifier;
        /**Boundaries of the detector element*/
        std::shared_ptr<const SurfaceBounds>           m_bounds;
        /** Thickness of this detector element*/
        double                                         m_thickness; //@TODO implement thickness from TGeoMode
        /**Corresponding Surface*/
        std::shared_ptr<const Surface>                 m_surface;
        /**possible contained surfaces*/
        std::vector<std::shared_ptr<const Surface>>    m_surfaces;

    };

    inline Identifier TGeoDetectorElement::identify() const {return m_identifier;}
    
    inline const Transform3D& TGeoDetectorElement::transform(const Identifier&) const {return (*m_transform);}
    
    inline const Surface& TGeoDetectorElement::surface(const Identifier&) const {return (*m_surface);}
    
    inline const std::vector< std::shared_ptr<const Surface> >& TGeoDetectorElement::surfaces() const {return (m_surfaces);}
    
    inline const SurfaceBounds& TGeoDetectorElement::bounds(const Identifier&) const {return (*m_bounds);}
    
    inline double TGeoDetectorElement::thickness() const { return m_thickness; }

}

#endif
