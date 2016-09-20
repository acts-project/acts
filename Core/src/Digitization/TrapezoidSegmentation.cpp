// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidSegmentation.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Digitization/TrapezoidSegmentation.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
#include "ACTS/Utilities/Helpers.hpp"

Acts::TrapezoidSegmentation::TrapezoidSegmentation(std::shared_ptr<const TrapezoidBounds> mBounds,
                                                   size_t numCellsX, size_t numCellsY) :
   m_activeBounds(mBounds),
   m_binUtility(nullptr),
   m_binsX(numCellsX),
   m_binsY(numCellsY)
{
    // first the x dimension if needed
    if (numCellsX > 1) {
         m_binUtility = std::make_unique<BinUtility>(numCellsX, 
                                                     -0.5*(mBounds->minHalflengthX()+mBounds->maxHalflengthX()),
                                                     0.5*(mBounds->minHalflengthX()+mBounds->maxHalflengthX()), 
                                                     open, 
                                                     binX);
    }
    // use y dimension if needed
    if (numCellsY > 1){
        BinUtility yBinUtility(numCellsY, -mBounds->halflengthY(), mBounds->halflengthY(), open, binY);
        if (m_binUtility)
            (*m_binUtility) += yBinUtility;
        else 
            m_binUtility = std::make_unique<BinUtility>(yBinUtility);
    }           
}       

Acts::TrapezoidSegmentation::~TrapezoidSegmentation()
{}

void 
Acts::TrapezoidSegmentation::createSegmenationSurfaces(std::vector< std::shared_ptr< const Surface> >& boundarySurfaces,
                                                       std::vector< std::shared_ptr< const Surface> >& segmentationSurfacesX,
                                                       std::vector< std::shared_ptr< const Surface> >& segmentationSurfacesY,
                                                       double halfThickness,
                                                       int readoutDirection,
                                                       double) const
{
    // The Lorentz angle is not taken into account for trapezoidal segmentation
    // (A) --- top/bottom surfaces -----------------------------------------------------------
    // let's create the top/botten surfaces first - we call them readout / counter readout
    // there are some things to consider 
    // - they share only the readout surface, then the segmentation surfaces are tilted and cannot be shared on the same module
    std::shared_ptr<const PlanarBounds> moduleBounds(new TrapezoidBounds(m_activeBounds->minHalflengthX(),
                                                                         m_activeBounds->maxHalflengthX(),
                                                                         m_activeBounds->halflengthY()));
    // - they are separated by half a thickness in z
    auto readoutPlaneTransform = std::make_shared<Transform3D>(Transform3D::Identity());
    auto counterPlaneTransform = std::make_shared<Transform3D>(Transform3D::Identity());
    // readout and counter readout bounds, the bounds of the readout plane are like the active ones
    std::shared_ptr<const PlanarBounds> readoutPlaneBounds = moduleBounds;
    std::shared_ptr<const PlanarBounds> counterPlaneBounds(nullptr);
    // the transform of the readout plane is always centric
    (*readoutPlaneTransform).translation()     = Vector3D(0.,0.,readoutDirection*halfThickness);
    // no lorentz angle and everything is straight-forward
    counterPlaneBounds = moduleBounds;
    (*counterPlaneTransform).translation()     = Vector3D(0.,0.,-readoutDirection*halfThickness);
    
    // - build the readout & counter readout surfaces
    boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(new PlaneSurface(readoutPlaneTransform,
                                                                                    readoutPlaneBounds)));
    boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(new PlaneSurface(counterPlaneTransform,
                                                                                    counterPlaneBounds)));
    
    // (B) - bin X -----------------------------------------------------------
    // easy stuff first, constant pitch size and 
    double pitchX =  2.*(m_activeBounds->maxHalflengthX()+m_activeBounds->minHalflengthX())*0.5/m_binsX;

    // now the rotation matrix for the xBins
    RotationMatrix3D xBinRotationMatrix;
    xBinRotationMatrix.col(0) = Vector3D::UnitY();
    xBinRotationMatrix.col(1) = Vector3D::UnitZ();
    xBinRotationMatrix.col(2) = Vector3D::UnitX();
    
    // reserve, it's always (number of bins-1) as the boundaries are within the boundarySurfaces
    segmentationSurfacesX.reserve(m_binsX);
    for (size_t ibinx = 0; ibinx <= m_binsX; ++ibinx){
      // the current step x position
      double cPosX = -(m_activeBounds->minHalflengthX()+m_activeBounds->maxHalflengthX())*0.5+ibinx*pitchX;
      
      // set position & rotation for all (boundaries and segmentations) --> Then you separate between them
      Vector3D xPosition = Vector3D(cPosX, 0.,0.);
      double stereoLocal = asin(sinStereoLocal(Vector2D(cPosX, 0.)));
      const RotationMatrix3D xRotation = xBinRotationMatrix*AngleAxis3D(stereoLocal, Vector3D::UnitY()); 
      // build the rotation from it
      auto binTransform = std::make_shared<Transform3D>(getTransformFromRotTransl(xRotation, xPosition));
      // the correct bounds for this
      auto xBinBounds = std::make_shared<const RectangleBounds>(m_activeBounds->halflengthY()/cos(stereoLocal),halfThickness);
      // these are the boundaries
      if (ibinx==0 || ibinx == m_binsX) // (i) this is the low/high boundary --- ( ibin == 0/m_binsX )
        boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(new PlaneSurface(binTransform,xBinBounds)));
      else // these are the bin boundaries
        segmentationSurfacesX.push_back(std::shared_ptr<const PlaneSurface>(new PlaneSurface(binTransform,xBinBounds)));
    }
    
    // (C) - bin Y surfaces - everything is defined -----------------------------------------------------------
    // now the rotation matrix for the yBins - anticyclic
    RotationMatrix3D yBinRotationMatrix;
    yBinRotationMatrix.col(0) = Vector3D::UnitX();
    yBinRotationMatrix.col(1) = Vector3D::UnitZ();
    yBinRotationMatrix.col(2) = Vector3D(0.,-1.,0.);
    // easy stuff first, constant pitch in Y 
    double pitchY             =  2.*m_activeBounds->halflengthY()/m_binsY;
    // reserve, it's always (number of bins-1) as the boundaries are within the boundarySurfaces
    segmentationSurfacesY.reserve(m_binsY);
    for (size_t ibiny = 0; ibiny <= m_binsY; ++ibiny){
      // the position of the bin surface
      double binPosY = -m_activeBounds->halflengthY()+ibiny*pitchY;
      Vector3D binSurfaceCenter(0.,binPosY,0.);
      double localPitchX = PitchX(Vector2D(0., binPosY));
      auto yBinBounds = std::make_shared<const RectangleBounds>(localPitchX*m_binsX*0.5,halfThickness);
      auto binTransform = std::make_shared<Transform3D>(getTransformFromRotTransl(yBinRotationMatrix,
                                                                                  binSurfaceCenter));
      // these are the boundaries
      if (ibiny == 0 || ibiny == m_binsY)
        boundarySurfaces.push_back(std::shared_ptr<PlaneSurface>(new PlaneSurface(binTransform,yBinBounds)));
      else // these are the bin boundaries
        segmentationSurfacesY.push_back(std::shared_ptr<PlaneSurface>(new PlaneSurface(binTransform,yBinBounds)));
    }  
}


const Acts::Vector2D
Acts::TrapezoidSegmentation::cellPosition(const DigitizationCell& dCell) const
{
    // use the bin utility for this job
    // @TODO fix it
    double bX = 0.; //(m_binsX>1) ? m_binUtility->binPosition(projectLocX(Vector2D(dCell.first, dCell.second)),0.,0) : 0.;
    double bY = 0.;// (m_binsY>1) ? m_binUtility->binPosition(dCell.second,0.,1) : 0.;
    return Vector2D(bX,bY);
}


/** Get the digitization cell from 3D position, it used the projection to the readout surface to estimate the 2D positon */
const Acts::DigitizationStep
Acts::TrapezoidSegmentation::digitizationStep(const Vector3D& startStep,
                                        const Vector3D& endStep,
                                        double halfThickness,
                                        int readoutDirection, 
                                        double lorentzAngle) const
{
   
    Vector3D stepCenter = 0.5*(startStep+endStep);
    // project to parameter surface
    double lorentzDeltaX = -readoutDirection*stepCenter.z()*tan(lorentzAngle);
    // take the full drift length
    double driftInZ = (halfThickness-readoutDirection*stepCenter.z());
    double driftLength  = fabs(driftInZ/cos(lorentzAngle)); 
    // the projected center
    Vector2D stepCenterProjected(stepCenter.x()+lorentzDeltaX,stepCenter.y());
    // the cell & its center
    DigitizationCell dCell = cell(stepCenterProjected);
    Vector2D cellCenter = cellPosition(dCell);
    // we are ready to return what we have
    return DigitizationStep((endStep-startStep).mag(),driftLength,dCell,startStep,endStep,stepCenterProjected,cellCenter);   
}
