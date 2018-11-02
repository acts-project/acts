// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CartesianSegmentation.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include <utility>

#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::CartesianSegmentation::CartesianSegmentation(
    const std::shared_ptr<const PlanarBounds>& mBounds,
    size_t                                     numCellsX,
    size_t                                     numCellsY)
  : m_activeBounds(mBounds), m_binUtility(nullptr)
{
  auto mutableBinUtility
      = std::make_shared<BinUtility>(numCellsX,
                                     -mBounds->boundingBox().halflengthX(),
                                     mBounds->boundingBox().halflengthX(),
                                     Acts::open,
                                     Acts::binX);
  (*mutableBinUtility) += BinUtility(numCellsY,
                                     -mBounds->boundingBox().halflengthY(),
                                     mBounds->boundingBox().halflengthY(),
                                     Acts::open,
                                     Acts::binY);
  m_binUtility = std::const_pointer_cast<const BinUtility>(mutableBinUtility);
}

Acts::CartesianSegmentation::CartesianSegmentation(
    std::shared_ptr<const BinUtility>   bUtility,
    std::shared_ptr<const PlanarBounds> mBounds)
  : m_activeBounds(std::move(mBounds)), m_binUtility(std::move(bUtility))
{
  if (!m_activeBounds) {
    m_activeBounds = std::make_shared<const RectangleBounds>(
        m_binUtility->max(0), m_binUtility->max(1));
  }
}

Acts::CartesianSegmentation::~CartesianSegmentation() = default;

void
Acts::CartesianSegmentation::createSegmentationSurfaces(
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
  std::shared_ptr<const PlanarBounds> moduleBounds(
      new RectangleBounds(m_activeBounds->boundingBox().halflengthX(),
                          m_activeBounds->boundingBox().halflengthY()));
  // - they are separated by half a thickness in z
  auto mutableReadoutPlaneTransform
      = std::make_shared<Transform3D>(Transform3D::Identity());
  auto mutableCounterPlaneTransform
      = std::make_shared<Transform3D>(Transform3D::Identity());
  // readout and counter readout bounds, the bounds of the readout plane are
  // like the active ones
  std::shared_ptr<const PlanarBounds> readoutPlaneBounds = moduleBounds;
  std::shared_ptr<const PlanarBounds> counterPlaneBounds(nullptr);
  // the transform of the readout plane is always centric
  (*mutableReadoutPlaneTransform).translation()
      = Vector3D(0., 0., readoutDirection * halfThickness);
  // no lorentz angle and everything is straight-forward
  if (lorentzAngle == 0.) {
    counterPlaneBounds = moduleBounds;
    (*mutableCounterPlaneTransform).translation()
        = Vector3D(0., 0., -readoutDirection * halfThickness);
  } else {
    // lorentz reduced Bounds
    double lorentzReducedHalfX = m_activeBounds->boundingBox().halflengthX()
        - fabs(lorentzPlaneShiftX);
    std::shared_ptr<const PlanarBounds> lorentzReducedBounds(
        new RectangleBounds(lorentzReducedHalfX,
                            m_activeBounds->boundingBox().halflengthY()));
    counterPlaneBounds = lorentzReducedBounds;
    // now we shift the counter plane in position - this depends on lorentz
    // angle
    double counterPlaneShift = -readoutDirection * lorentzPlaneShiftX;
    (*mutableCounterPlaneTransform).translation()
        = Vector3D(counterPlaneShift, 0., -readoutDirection * halfThickness);
  }
  // - finalize the transforms
  auto readoutPlaneTransform = std::const_pointer_cast<const Transform3D>(
      mutableReadoutPlaneTransform);
  auto counterPlaneTransform = std::const_pointer_cast<const Transform3D>(
      mutableCounterPlaneTransform);
  // - build the readout & counter readout surfaces
  boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
      readoutPlaneTransform, readoutPlaneBounds));
  boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
      counterPlaneTransform, counterPlaneBounds));

  // (B) - bin X and lorentz surfaces
  // -----------------------------------------------------------
  // easy stuff first, constant pitch size and
  double pitchX = 2. * m_activeBounds->boundingBox().halflengthX()
      / m_binUtility->bins(0);

  // now, let's create the shared bounds of all surfaces marking x bins - choice
  // fixes orientation of the matrix
  std::shared_ptr<const PlanarBounds> xBinBounds(new RectangleBounds(
      m_activeBounds->boundingBox().halflengthY(), halfThickness));
  // now, let's create the shared bounds of all surfaces marking lorentz planes
  double lorentzPlaneHalfX = std::abs(halfThickness / cos(lorentzAngle));
  // the bounds of the lorentz plane
  std::shared_ptr<const PlanarBounds> lorentzPlaneBounds = (lorentzAngle == 0.)
      ? xBinBounds
      : std::shared_ptr<const PlanarBounds>(new RectangleBounds(
            m_activeBounds->boundingBox().halflengthY(), lorentzPlaneHalfX));

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
  segmentationSurfacesX.reserve(m_binUtility->bins(0));
  // create and fill them
  for (size_t ibinx = 0; ibinx <= m_binUtility->bins(0); ++ibinx) {
    // the current step x position
    double cPosX
        = -m_activeBounds->boundingBox().halflengthX() + ibinx * pitchX;
    // (i) this is the low/high boundary --- ( ibin == 0/m_binUtility->bins(0) )
    if ((ibinx == 0u) || ibinx == m_binUtility->bins(0)) {
      // check if it a straight boundary or not: always straight for no lorentz
      // angle,
      // and either the first boundary or the last dependening on lorentz &
      // readout
      bool boundaryStraight
          = (lorentzAngle == 0.
             || ((ibinx == 0u) && readoutDirection * lorentzAngle > 0.)
             || (ibinx == m_binUtility->bins(0)
                 && readoutDirection * lorentzAngle < 0));
      // set the low boundary parameters : position & rotation
      Vector3D boundaryXPosition = boundaryStraight
          ? Vector3D(cPosX, 0., 0.)
          : Vector3D(cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      // rotation of the boundary: striaght or lorentz
      const RotationMatrix3D& boundaryXRotation
          = boundaryStraight ? xBinRotationMatrix : lorentzPlaneRotationMatrix;
      // build the rotation from it
      auto boundaryXTransform = std::make_shared<const Transform3D>(
          Translation3D(boundaryXPosition) * boundaryXRotation);
      // the correct bounds for this
      std::shared_ptr<const PlanarBounds> boundaryXBounds
          = boundaryStraight ? xBinBounds : lorentzPlaneBounds;
      // boundary surfaces
      boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
          boundaryXTransform, boundaryXBounds));
      // (ii) this is the in between bins  --- ( 1 <= ibin < m_mbnsX )
    } else {
      // shift by the lorentz angle
      Vector3D lorentzPlanePosition(
          cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      auto lorentzPlaneTransform = std::make_shared<const Transform3D>(
          Translation3D(lorentzPlanePosition) * lorentzPlaneRotationMatrix);
      // lorentz plane surfaces
      segmentationSurfacesX.push_back(Surface::makeShared<PlaneSurface>(
          lorentzPlaneTransform, lorentzPlaneBounds));
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
  double pitchY = 2. * m_activeBounds->boundingBox().halflengthY()
      / m_binUtility->bins(1);
  // let's create the shared bounds of all surfaces marking y bins
  std::shared_ptr<const PlanarBounds> yBinBounds(new RectangleBounds(
      m_activeBounds->boundingBox().halflengthX(), halfThickness));
  // reserve, it's always (number of bins-1) as the boundaries are within the
  // boundarySurfaces
  segmentationSurfacesY.reserve(m_binUtility->bins(1));
  for (size_t ibiny = 0; ibiny <= m_binUtility->bins(1); ++ibiny) {
    // the position of the bin surface
    double binPosY
        = -m_activeBounds->boundingBox().halflengthY() + ibiny * pitchY;
    Vector3D binSurfaceCenter(0., binPosY, 0.);
    // the binning transform
    auto binTransform = std::make_shared<const Transform3D>(
        Translation3D(binSurfaceCenter) * yBinRotationMatrix);
    // these are the boundaries
    if (ibiny == 0 || ibiny == m_binUtility->bins(1)) {
      boundarySurfaces.push_back(
          Surface::makeShared<PlaneSurface>(binTransform, yBinBounds));
    } else {  // these are the bin boundaries
      segmentationSurfacesY.push_back(
          Surface::makeShared<PlaneSurface>(binTransform, yBinBounds));
    }
  }
}

Acts::Vector2D
Acts::CartesianSegmentation::cellPosition(const DigitizationCell& dCell) const
{
  double bX = m_binUtility->bins(0) > 1
      ? m_binUtility->binningData()[0].center(dCell.channel0)
      : 0.;
  double bY = m_binUtility->bins(1) > 1
      ? m_binUtility->binningData()[1].center(dCell.channel1)
      : 0.;
  return Vector2D(bX, bY);
}

/** Get the digitization cell from 3D position, it used the projection to the
 * readout surface to estimate the 2D positon */
Acts::DigitizationStep
Acts::CartesianSegmentation::digitizationStep(const Vector3D& startStep,
                                              const Vector3D& endStep,
                                              double          halfThickness,
                                              int             readoutDirection,
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
  return DigitizationStep((endStep - startStep).norm(),
                          driftLength,
                          dCell,
                          startStep,
                          endStep,
                          stepCenterProjected,
                          cellCenter);
}
