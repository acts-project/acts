// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderLayer.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Layers/CylinderLayer.hpp"

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Detector/GenericApproachDescriptor.hpp"
#include "ACTS/Detector/GenericOverlapDescriptor.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"
#include "ACTS/Volumes/BoundarySurfaceFace.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
// Core module

        
Acts::CylinderLayer::CylinderLayer(std::shared_ptr<Acts::Transform3D> transform,
                                  std::shared_ptr<const Acts::CylinderBounds> cBounds,
                                  std::unique_ptr<SurfaceArray> surfaceArray,
                                  double thickness,
                                  Acts::OverlapDescriptor* olap,
                                  Acts::ApproachDescriptor* ades,
                                  int laytyp) :
  CylinderSurface(transform, cBounds),
  Layer(std::move(surfaceArray), thickness, olap, ades, laytyp)
{
    
    // just create a generic overlap descriptor if none is there
    if (!Layer::m_overlapDescriptor) Layer::m_overlapDescriptor = new GenericOverlapDescriptor();
    // create the representing volume
    CylinderVolumeBounds* cvBounds = new CylinderVolumeBounds(cBounds->r()-0.5*thickness,
                                                              cBounds->r()+0.5*thickness,
                                                              cBounds->halflengthZ());
    Layer::m_representingVolume = new AbstractVolume(transform,VolumeBoundsPtr(cvBounds));
    // associate teh layer 
    CylinderSurface::associateLayer(*this);
    // an approach descriptor is automatically created if there's a surface array
    if (!ades && Layer::m_surfaceArray)
        buildApproachDescriptor();
    // register the layer if the approach descriptor was provided
    if (ades) m_approachDescriptor->registerLayer(*this);
    
}
                        
Acts::CylinderLayer::CylinderLayer(const Acts::CylinderLayer& clay, const Acts::Transform3D& transf):
  CylinderSurface(clay,transf),
  Layer(clay)
{
    if (m_surfaceArray) buildApproachDescriptor();
}
   
const Acts::CylinderSurface& Acts::CylinderLayer::surfaceRepresentation() const
{
  return (*this);
}

/** build approach surfaces */
void Acts::CylinderLayer::buildApproachDescriptor() const {
    // delete it
    delete m_approachDescriptor;
    // delete the surfaces    
    // take the boundary surfaces of the representving volume if they exist
    if (m_representingVolume){
        // get teh boundary surfaces
        const std::vector< std::shared_ptr<const Acts::BoundarySurface<Acts::AbstractVolume> > >& bSurfaces = m_representingVolume->boundarySurfaces();
        // fill in the surfaces into the vector
        std::vector< std::shared_ptr<const Acts::BoundarySurface<Acts::AbstractVolume> > > aSurfaces;
        if (bSurfaces.size() > size_t(tubeOuterCover)) aSurfaces.push_back(bSurfaces.at(tubeInnerCover));
        aSurfaces.push_back(bSurfaces.at(tubeOuterCover));
        // create an ApproachDescriptor with Boundary surfaces
        m_approachDescriptor = new Acts::GenericApproachDescriptor<const BoundarySurface<AbstractVolume> >(aSurfaces);   
    } else {        
        // create the new surfaces
        std::vector< const Acts::Surface* > aSurfaces;              
        aSurfaces.push_back(new Acts::CylinderSurface(m_transform, m_bounds->r()-0.5*thickness(), m_bounds->halflengthZ()));
        aSurfaces.push_back(new Acts::CylinderSurface(m_transform, m_bounds->r()+0.5*thickness(), m_bounds->halflengthZ()));
        // create an ApproachDescriptor with standard surfaces surfaces - these will be deleted by the approach descriptor
        m_approachDescriptor = new Acts::GenericApproachDescriptor<const Acts::Surface>(aSurfaces);   
    }
    for (auto& sIter : (m_approachDescriptor->containedSurfaces())){
        if (sIter)
            sIter->associateLayer(*this);
    }
}
