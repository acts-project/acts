// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DIGITIZATION_DIGITIZATIONMODULE_H
#define ACTS_DIGITIZATION_DIGITIZATIONMODULE_H

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Digitization/DigitizationCell.hpp"
#include "ACTS/Digitization/Segmentation.hpp"
#include <memory>

namespace Acts {

    class Surface;
    typedef std::shared_ptr<const Surface> SurfacePtr;
    typedef std::vector< SurfacePtr > SurfacePtrVector; 
    
    /// @class DigitizationModule
    ///
    /// Class that holds the surfaces for a planar digitization detector module.
    /// 
    /// It needs a descriptor to design different pixel/strixels/strip setups (with a segmentation class)
    ///
    /// The digitizaiton is done in the local frame of the surface and binning can be done
    /// in local x and local y direction. 
    ///
    /// The lorentz angle is assumed to be only in x-direction and constant for the module, 
    /// it is measured from the local z-direction towards the local x-direction.
    ///
    /// The readout direction defines the charge drift either:
    /// a) towards the surface at -halfThickness if readout is defined at -1 
    /// b) towards the surface at +halfThickness if readout is defined at +1
    ///
    /// Conventions: 
    ///   - 3D positions are within the 3D frame of the module 
    ///   - 2D positions are corrected to parameter surface  at the center of the module (and not the readout surface) 
    ///
    /// The lorenzShift is the correction from the readout surface to the parameter surface
    ///     
    class DigitizationModule { 
        public :
            /// Constructor from a Segmentation descriptor 
            ///
            /// @param moduleSegmentation is the segmentation descriptions
            /// @param halfThickness is the half thickness of the module
            /// @param readoutDirection is the readout drift direction
            /// lorentz angle is the lorentz drift angle 
            DigitizationModule(std::shared_ptr<const Segmentation> moduleSegmentation,
                               double halfThickness,
                               int readoutDirection,
                               double lorentzAngle);

            /// Virtual Destructor 
            virtual ~DigitizationModule(){}
            
            /// Return the internal test segmentation surfaces to test between entry 
            /// and exit given by their cell id's - the boundaries are not given
            ///
            /// @param entryCids are the entry digitisation cell ids
            /// @param exitCids are the exit digitisation cell ids
            /// @return object is a vector of shared surfaces  
            const SurfacePtrVector 
            segmentationSurfaces(const DigitizationCell& entryCids, const DigitizationCell& exitCids) const;
            
            /// Get the digitization cell from a position 
            const DigitizationCell 
            cell(const Vector2D& position) const;
            
            /// module thickness 
            double 
            halfThickness() const;
            
            /// return the readout direction 
            int 
            readoutDirection() const;
            
            /// return the lorentz Angle */
            double 
            lorentzAngle() const;
            
            /// return the segmenation 
            const Segmentation& 
            segmentation() const;
            
            /// return the test surfaces between these points
            ///
            /// @param start is the start position of the step
            /// @param end is the end position of the step
            /// @return stepSurfaces are the surfaces to test
            const SurfacePtrVector 
            stepSurfaces(const Vector3D& start, const Vector3D& end) const;
            
            /// Fill the associated digitsation cell from this start and end position, 
            // correct for lorentz effect if needed 
            const DigitizationStep 
            digitizationStep(const Vector3D& start, const Vector3D& end) const;
            
            /// return the bounding surfaces inlcuding top and bottom 
            const SurfacePtrVector& 
            boundarySurfaces() const;
            
            /// return all surfaces in X - excluding the boundaries 
            const SurfacePtrVector& 
            segmentationSurfacesX() const;
            
            /// return all surfaces in Y - excluding the boundaries 
            const SurfacePtrVector& segmentationSurfacesY() const;
                        
        private:
            double  m_halfThickness;
            int     m_readoutDirection;      ///< defines if the readout is along (+1) / (-1) wrt the z axis
            double  m_lorentzAngle;          ///< the lorentz angle
            double  m_tanLorentzAngle;       ///< and the tangent of it
            
            std::shared_ptr<const Segmentation> m_segmentation;  /// segmentation descriptor            
            SurfacePtrVector  m_boundarySurfaces;      ///< boundary surfaces z, x, y 
            SurfacePtrVector  m_segmentationSurfacesX; ///< segmentation surfaces in X - without boundaries
            SurfacePtrVector  m_segmentationSurfacesY; ///< segmentation surfaces in Y - without boundaries
    
    };

    inline double 
    DigitizationModule::halfThickness() const
        { return m_halfThickness; }
    
    inline int 
    DigitizationModule::readoutDirection() const
        { return m_readoutDirection; }
    
    inline double 
    DigitizationModule::lorentzAngle() const
        { return m_lorentzAngle; }
    
    inline const Segmentation& 
    DigitizationModule::segmentation() const
        { return (*(m_segmentation.get())); }
    
    inline const SurfacePtrVector& 
    DigitizationModule::boundarySurfaces() const
        { return m_boundarySurfaces; }
    
    inline const SurfacePtrVector& 
    DigitizationModule::segmentationSurfacesX() const
        { return m_segmentationSurfacesX; }
    
    inline const SurfacePtrVector& 
    DigitizationModule::segmentationSurfacesY() const
        { return m_segmentationSurfacesY; }
    
    inline const DigitizationStep 
    DigitizationModule::digitizationStep(const Vector3D& start, const Vector3D& end) const
        { return m_segmentation->digitizationStep(start,end,m_halfThickness,m_readoutDirection,m_lorentzAngle); }

}

#endif
