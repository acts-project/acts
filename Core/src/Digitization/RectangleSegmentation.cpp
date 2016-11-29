// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangularSegmentation.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Digitization/RectangularSegmentation.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Helpers.hpp"
#include <cmath>

Acts::RectangularSegmentation::RectangularSegmentation(
    std::shared_ptr<const Acts::RectangleBounds> mBounds,
    size_t                                       numCellsX,
    size_t                                       numCellsY)
  : m_activeBounds(mBounds)
  , m_binUtility(nullptr)
  , m_binsX(numCellsX)
  , m_binsY(numCellsY)
{
  // first the x dimension if needed
  if (numCellsX > 1)
    m_binUtility = std::make_unique<BinUtility>(numCellsX,
                                                -mBounds->halflengthX(),
                                                mBounds->halflengthX(),
                                                Acts::open,
                                                Acts::binX);
  // use y dimension if needed
  if (numCellsY > 1) {
    BinUtility yBinUtility(numCellsY,
                           -mBounds->halflengthY(),
                           mBounds->halflengthY(),
                           Acts::open,
                           Acts::binY);
    if (m_binUtility)
      (*m_binUtility) += yBinUtility;
    else
      m_binUtility = std::make_unique<BinUtility>(yBinUtility);
  }
}

Acts::RectangularSegmentation::~RectangularSegmentation()
{
}

void
Acts::RectangularSegmentation::createSegmentationSurfaces(
    SurfacePtrVector& boundarySurfaces,
    SurfacePtrVector& segmentationSurfacesX,
    SurfacePtrVector& segmentationSurfacesY,
    double            halfThickness,
    int               readoutDirection,
    double            lorentzAngle) const
{
  // may be needed throughout
  double lorentzAngleTan    = tan(lorentzAngle);
  double lorentzPlaneShiftX = halfThickness * lorentzAngleTan;

  // (A) --- top/bottom surfaces
  // -----------------------------------------------------------
  // let's create the top/botten surfaces first - we call them readout / counter
  // readout
  // there are some things to consider
  // - they share the RectangleBounds only if the lorentzAngle is 0
  // otherwise only the readout surface has full length bounds like the module
  std::shared_ptr<const PlanarBounds> moduleBounds(new RectangleBounds(
      m_activeBounds->halflengthX(), m_activeBounds->halflengthY()));
  // - they are separated by half a thickness in z
  auto readoutPlaneTransform
      = std::make_shared<Transform3D>(Transform3D::Identity());
  auto counterPlaneTransform
      = std::make_shared<Transform3D>(Transform3D::Identity());
  // readout and counter readout bounds, the bounds of the readout plane are
  // like the active ones
  std::shared_ptr<const PlanarBounds> readoutPlaneBounds = moduleBounds;
  std::shared_ptr<const PlanarBounds> counterPlaneBounds(nullptr);
  // the transform of the readout plane is always centric
  (*readoutPlaneTransform).translation()
      = Vector3D(0., 0., readoutDirection * halfThickness);
  // no lorentz angle and everything is straight-forward
  if (lorentzAngle == 0.) {
    counterPlaneBounds = moduleBounds;
    (*counterPlaneTransform).translation()
        = Vector3D(0., 0., -readoutDirection * halfThickness);
  } else {
    // lorentz reduced Bounds
    double lorentzReducedHalfX
        = m_activeBounds->halflengthX() - std::abs(lorentzPlaneShiftX);
    std::shared_ptr<const PlanarBounds> lorentzReducedBounds(
        new RectangleBounds(lorentzReducedHalfX,
                            m_activeBounds->halflengthY()));
    counterPlaneBounds = lorentzReducedBounds;
    // now we shift the counter plane in position - this depends on lorentz
    // angle
    double counterPlaneShift = -readoutDirection * lorentzPlaneShiftX;
    (*counterPlaneTransform).translation()
        = Vector3D(counterPlaneShift, 0., -readoutDirection * halfThickness);
  }
  // - build the readout & counter readout surfaces
  boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(
      new PlaneSurface(readoutPlaneTransform, readoutPlaneBounds)));
  boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(
      new PlaneSurface(counterPlaneTransform, counterPlaneBounds)));

  // (B) - bin X and lorentz surfaces
  // -----------------------------------------------------------
  // easy stuff first, constant pitch size and
  double pitchX = 2. * m_activeBounds->halflengthX() / m_binsX;

  // now, let's create the shared bounds of all surfaces marking x bins - choice
  // fixes orientation of the matrix
  std::shared_ptr<const PlanarBounds> xBinBounds(
      new RectangleBounds(m_activeBounds->halflengthY(), halfThickness));
  // now, let's create the shared bounds of all surfaces marking lorentz planes
  double lorentzPlaneHalfX = std::abs(halfThickness / cos(lorentzAngle));
  // teh bounds of the lorentz plane
  std::shared_ptr<const PlanarBounds> lorentzPlaneBounds = (lorentzAngle == 0.)
      ? xBinBounds
      : std::shared_ptr<const PlanarBounds>(new RectangleBounds(
            m_activeBounds->halflengthY(), lorentzPlaneHalfX));

  // now the rotation matrix for the xBins
  RotationMatrix3D xBinRotationMatrix;
  xBinRotationMatrix.col(0) = Vector3D::UnitY();
  xBinRotationMatrix.col(1) = Vector3D::UnitZ();
  xBinRotationMatrix.col(2) = Vector3D::UnitX();
  // now the lorentz plane rotation should be the xBin rotation, rotated by the
  // lorentz angle around y
  RotationMatrix3D lorentzPlaneRotationMatrix = (lorentzAngle != 0.)
      ? xBinRotationMatrix * AngleAxis3D(lorentzAngle, Vector3D::UnitX())
      : xBinRotationMatrix;

  // reserve, it's always (number of bins-1) as the boundaries are within the
  // boundarySurfaces
  segmentationSurfacesX.reserve(m_binsX);
  // create and fill them
  for (size_t ibinx = 0; ibinx <= m_binsX; ++ibinx) {
    // the current step x position
    double cPosX = -m_activeBounds->halflengthX() + ibinx * pitchX;
    // (i) this is the low/high boundary --- ( ibin == 0/m_binsX )
    if (!ibinx || ibinx == m_binsX) {
      // check if it a straight boundary or not: always straight for no lorentz
      // angle,
      // and either the first boundary or the last dependening on lorentz &
      // readout
      bool boundaryStraight
          = (lorentzAngle == 0.
             || (!ibinx && readoutDirection * lorentzAngle > 0.)
             || (ibinx == m_binsX && readoutDirection * lorentzAngle < 0));
      // set the low boundary parameters : position & rotation
      Vector3D boundaryXPosition = boundaryStraight
          ? Vector3D(cPosX, 0., 0.)
          : Vector3D(cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      // rotation of the boundary: striaght or lorentz
      const RotationMatrix3D& boundaryXRotation
          = boundaryStraight ? xBinRotationMatrix : lorentzPlaneRotationMatrix;
      // build the rotation from it
      auto boundaryXTransform = std::make_shared<Transform3D>(
          getTransformFromRotTransl(boundaryXRotation, boundaryXPosition));
      // the correct bounds for this
      std::shared_ptr<const PlanarBounds> boundaryXBounds
          = boundaryStraight ? xBinBounds : lorentzPlaneBounds;
      // boundary surfaces
      boundarySurfaces.push_back(std::shared_ptr<const PlaneSurface>(
          new PlaneSurface(boundaryXTransform, boundaryXBounds)));
      // (ii) this is the in between bins  --- ( 1 <= ibin < m_mbnsX )
    } else {
      // shift by the lorentz angle
      Vector3D lorentzPlanePosition(
          cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      auto lorentzPlaneTransform
          = std::make_shared<Transform3D>(getTransformFromRotTransl(
              lorentzPlaneRotationMatrix, lorentzPlanePosition));
      // lorentz plane surfaces
      segmentationSurfacesX.push_back(std::shared_ptr<const PlaneSurface>(
          new PlaneSurface(lorentzPlaneTransform, lorentzPlaneBounds)));
    }
  }

  // (C) - bin Y surfaces - everything is defined
  // -----------------------------------------------------------
  // now the rotation matrix for the yBins - anticyclic
  RotationMatrix3D yBinRotationMatrix;
  yBinRotationMatrix.col(0) = Vector3D::UnitX();
  yBinRotationMatrix.col(1) = Vector3D::UnitZ();
  yBinRotationMatrix.col(2) = Vector3D(0., -1., 0.);
  // easy stuff first, constant pitch in Y
  double pitchY = 2. * m_activeBounds->halflengthY() / m_binsY;
  // let's create the shared bounds of all surfaces marking y bins
  std::shared_ptr<const PlanarBounds> yBinBounds(
      new RectangleBounds(m_activeBounds->halflengthX(), halfThickness));
  // reserve, it's always (number of bins-1) as the boundaries are within the
  // boundarySurfaces
  segmentationSurfacesY.reserve(m_binsY);
  for (size_t ibiny = 0; ibiny <= m_binsY; ++ibiny) {
    // the position of the bin surface
    double   binPosY = -m_activeBounds->halflengthY() + ibiny * pitchY;
    Vector3D binSurfaceCenter(0., binPosY, 0.);
    // the binning transform
    auto binTransform = std::make_shared<Transform3D>(
        getTransformFromRotTransl(yBinRotationMatrix, binSurfaceCenter));
    // these are the boundaries
    if (ibiny == 0 || ibiny == m_binsY)
      boundarySurfaces.push_back(std::shared_ptr<PlaneSurface>(
          new PlaneSurface(binTransform, yBinBounds)));
    else  // these are the bin boundaries
      segmentationSurfacesY.push_back(std::shared_ptr<PlaneSurface>(
          new PlaneSurface(binTransform, yBinBounds)));
  }
}

const Acts::Vector2D
Acts::RectangularSegmentation::cellPosition(const DigitizationCell& dCell) const
{
  // @TODO add protection agains 1D binUtility for Y
  double bX
      = m_binsX > 1 ? m_binUtility->binningData()[0].center(dCell.first) : 0.;
  double bY
      = m_binsY > 1 ? m_binUtility->binningData()[1].center(dCell.second) : 0.;
  return Vector2D(bX, bY);

  return Vector2D(bX, bY);
}

/** Get the digitization cell from 3D position, it used the projection to the
 * readout surface to estimate the 2D positon */
const Acts::DigitizationStep
Acts::RectangularSegmentation::digitizationStep(const Vector3D& startStep,
                                                const Vector3D& endStep,
                                                double          halfThickness,
                                                int    readoutDirection,
                                                double lorentzAngle) const
{
  Vector3D stepCenter = 0.5 * (startStep + endStep);
  // take the full drift length
  // this is the absolute drift in z
  double driftInZ = halfThickness - readoutDirection * stepCenter.z();
  // this is the absolute drift length
  double driftLength = driftInZ / cos(lorentzAngle);
  // project to parameter the readout surface
  double lorentzDeltaX = readoutDirection * driftInZ * tan(lorentzAngle);
  // the projected center, it has the lorentz shift applied
  Vector2D stepCenterProjected(stepCenter.x() + lorentzDeltaX, stepCenter.y());
  // the cell & its center
  Acts::DigitizationCell dCell      = cell(stepCenterProjected);
  Vector2D               cellCenter = cellPosition(dCell);
  // we are ready to return what we have
  return DigitizationStep((endStep - startStep).mag(),
                          driftLength,
                          dCell,
                          startStep,
                          endStep,
                          stepCenterProjected,
                          cellCenter);
}
