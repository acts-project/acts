// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CartesianSegmentation.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Digitization/CartesianSegmentation.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

Acts::CartesianSegmentation::CartesianSegmentation(
    const std::shared_ptr<const PlanarBounds>& mBounds, std::size_t numCellsX,
    std::size_t numCellsY)
    : m_activeBounds(mBounds), m_binUtility(nullptr) {
  auto mutableBinUtility = std::make_shared<BinUtility>(
      numCellsX, -mBounds->boundingBox().halfLengthX(),
      mBounds->boundingBox().halfLengthX(), Acts::open, Acts::binX);
  (*mutableBinUtility) +=
      BinUtility(numCellsY, -mBounds->boundingBox().halfLengthY(),
                 mBounds->boundingBox().halfLengthY(), Acts::open, Acts::binY);
  m_binUtility = std::const_pointer_cast<const BinUtility>(mutableBinUtility);
}

Acts::CartesianSegmentation::CartesianSegmentation(
    std::shared_ptr<const BinUtility> bUtility,
    std::shared_ptr<const PlanarBounds> mBounds)
    : m_activeBounds(std::move(mBounds)), m_binUtility(std::move(bUtility)) {
  if (!m_activeBounds) {
    m_activeBounds = std::make_shared<const RectangleBounds>(
        m_binUtility->max(0), m_binUtility->max(1));
  }
}

Acts::CartesianSegmentation::~CartesianSegmentation() = default;

void Acts::CartesianSegmentation::createSegmentationSurfaces(
    SurfacePtrVector& boundarySurfaces, SurfacePtrVector& segmentationSurfacesX,
    SurfacePtrVector& segmentationSurfacesY, double halfThickness,
    int readoutDirection, double lorentzAngle) const {
  // may be needed throughout
  double lorentzAngleTan = tan(lorentzAngle);
  double lorentzPlaneShiftX = halfThickness * lorentzAngleTan;

  // (A) --- top/bottom surfaces
  // -----------------------------------------------------------
  // let's create the top/botten surfaces first - we call them readout / counter
  // readout
  // there are some things to consider
  // - they share the RectangleBounds only if the lorentzAngle is 0
  // otherwise only the readout surface has full length bounds like the module
  std::shared_ptr<const PlanarBounds> moduleBounds(
      new RectangleBounds(m_activeBounds->boundingBox()));
  // - they are separated by half a thickness in z
  auto readoutPlaneTransform = Transform3::Identity();
  auto counterPlaneTransform = Transform3::Identity();
  // readout and counter readout bounds, the bounds of the readout plane are
  // like the active ones
  std::shared_ptr<const PlanarBounds> readoutPlaneBounds = moduleBounds;
  std::shared_ptr<const PlanarBounds> counterPlaneBounds(nullptr);
  // the transform of the readout plane is always centric
  readoutPlaneTransform.translation() =
      Vector3(0., 0., readoutDirection * halfThickness);
  // no lorentz angle and everything is straight-forward
  if (lorentzAngle == 0.) {
    counterPlaneBounds = moduleBounds;
    counterPlaneTransform.translation() =
        Vector3(0., 0., -readoutDirection * halfThickness);
  } else {
    // lorentz reduced Bounds
    double lorentzReducedHalfX =
        m_activeBounds->boundingBox().halfLengthX() - fabs(lorentzPlaneShiftX);
    std::shared_ptr<const PlanarBounds> lorentzReducedBounds(
        new RectangleBounds(lorentzReducedHalfX,
                            m_activeBounds->boundingBox().halfLengthY()));
    counterPlaneBounds = lorentzReducedBounds;
    // now we shift the counter plane in position - this depends on lorentz
    // angle
    double counterPlaneShift = -readoutDirection * lorentzPlaneShiftX;
    counterPlaneTransform.translation() =
        Vector3(counterPlaneShift, 0., -readoutDirection * halfThickness);
  }
  // - build the readout & counter readout surfaces
  boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
      readoutPlaneTransform, readoutPlaneBounds));
  boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
      counterPlaneTransform, counterPlaneBounds));

  // (B) - bin X and lorentz surfaces
  // -----------------------------------------------------------
  // easy stuff first, constant pitch size and
  double pitchX =
      2. * m_activeBounds->boundingBox().halfLengthX() / m_binUtility->bins(0);

  // now, let's create the shared bounds of all surfaces marking x bins - choice
  // fixes orientation of the matrix
  std::shared_ptr<const PlanarBounds> xBinBounds(new RectangleBounds(
      m_activeBounds->boundingBox().halfLengthY(), halfThickness));
  // now, let's create the shared bounds of all surfaces marking lorentz planes
  double lorentzPlaneHalfX = std::abs(halfThickness / cos(lorentzAngle));
  // the bounds of the lorentz plane
  std::shared_ptr<const PlanarBounds> lorentzPlaneBounds =
      (lorentzAngle == 0.)
          ? xBinBounds
          : std::shared_ptr<const PlanarBounds>(
                new RectangleBounds(m_activeBounds->boundingBox().halfLengthY(),
                                    lorentzPlaneHalfX));

  // now the rotation matrix for the xBins
  RotationMatrix3 xBinRotationMatrix;
  xBinRotationMatrix.col(0) = Vector3::UnitY();
  xBinRotationMatrix.col(1) = Vector3::UnitZ();
  xBinRotationMatrix.col(2) = Vector3::UnitX();
  // now the lorentz plane rotation should be the xBin rotation, rotated by the
  // lorentz angle around y
  RotationMatrix3 lorentzPlaneRotationMatrix =
      (lorentzAngle != 0.)
          ? xBinRotationMatrix * AngleAxis3(lorentzAngle, Vector3::UnitX())
          : xBinRotationMatrix;

  // reserve, it's always (number of bins-1) as the boundaries are within the
  // boundarySurfaces
  segmentationSurfacesX.reserve(m_binUtility->bins(0));
  // create and fill them
  for (std::size_t ibinx = 0; ibinx <= m_binUtility->bins(0); ++ibinx) {
    // the current step x position
    double cPosX =
        -m_activeBounds->boundingBox().halfLengthX() + ibinx * pitchX;
    // (i) this is the low/high boundary --- ( ibin == 0/m_binUtility->bins(0) )
    if ((ibinx == 0u) || ibinx == m_binUtility->bins(0)) {
      // check if it is a straight boundary or not: always straight for no
      // lorentz angle, and either the first boundary or the last depending on
      // lorentz and readout
      bool boundaryStraight =
          (lorentzAngle == 0. ||
           ((ibinx == 0u) && readoutDirection * lorentzAngle > 0.) ||
           (ibinx == m_binUtility->bins(0) &&
            readoutDirection * lorentzAngle < 0));
      // set the low boundary parameters : position & rotation
      Vector3 boundaryXPosition =
          boundaryStraight
              ? Vector3(cPosX, 0., 0.)
              : Vector3(cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      // rotation of the boundary: straight or lorentz
      const RotationMatrix3& boundaryXRotation =
          boundaryStraight ? xBinRotationMatrix : lorentzPlaneRotationMatrix;
      // build the rotation from it
      auto boundaryXTransform =
          Transform3(Translation3(boundaryXPosition) * boundaryXRotation);
      // the correct bounds for this
      std::shared_ptr<const PlanarBounds> boundaryXBounds =
          boundaryStraight ? xBinBounds : lorentzPlaneBounds;
      // boundary surfaces
      boundarySurfaces.push_back(Surface::makeShared<PlaneSurface>(
          boundaryXTransform, boundaryXBounds));
      // (ii) this is the in between bins  --- ( 1 <= ibin < m_mbnsX )
    } else {
      // shift by the lorentz angle
      Vector3 lorentzPlanePosition(
          cPosX - readoutDirection * lorentzPlaneShiftX, 0., 0.);
      auto lorentzPlaneTransform = Transform3(
          Translation3(lorentzPlanePosition) * lorentzPlaneRotationMatrix);
      // lorentz plane surfaces
      segmentationSurfacesX.push_back(Surface::makeShared<PlaneSurface>(
          lorentzPlaneTransform, lorentzPlaneBounds));
    }
  }

  // (C) - bin Y surfaces - everything is defined
  // -----------------------------------------------------------
  // now the rotation matrix for the yBins - anticyclic
  RotationMatrix3 yBinRotationMatrix;
  yBinRotationMatrix.col(0) = Vector3::UnitX();
  yBinRotationMatrix.col(1) = Vector3::UnitZ();
  yBinRotationMatrix.col(2) = Vector3(0., -1., 0.);
  // easy stuff first, constant pitch in Y
  double pitchY =
      2. * m_activeBounds->boundingBox().halfLengthY() / m_binUtility->bins(1);
  // let's create the shared bounds of all surfaces marking y bins
  std::shared_ptr<const PlanarBounds> yBinBounds(new RectangleBounds(
      m_activeBounds->boundingBox().halfLengthX(), halfThickness));
  // reserve, it's always (number of bins-1) as the boundaries are within the
  // boundarySurfaces
  segmentationSurfacesY.reserve(m_binUtility->bins(1));
  for (std::size_t ibiny = 0; ibiny <= m_binUtility->bins(1); ++ibiny) {
    // the position of the bin surface
    double binPosY =
        -m_activeBounds->boundingBox().halfLengthY() + ibiny * pitchY;
    Vector3 binSurfaceCenter(0., binPosY, 0.);
    // the binning transform
    auto binTransform =
        Transform3(Translation3(binSurfaceCenter) * yBinRotationMatrix);
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

Acts::Vector2 Acts::CartesianSegmentation::cellPosition(
    const DigitizationCell& dCell) const {
  double bX = m_binUtility->bins(0) > 1
                  ? m_binUtility->binningData()[0].center(dCell.channel0)
                  : 0.;
  double bY = m_binUtility->bins(1) > 1
                  ? m_binUtility->binningData()[1].center(dCell.channel1)
                  : 0.;
  return Vector2(bX, bY);
}

/** Get the digitization cell from 3D position, it used the projection to the
 * readout surface to estimate the 2D position */
Acts::DigitizationStep Acts::CartesianSegmentation::digitizationStep(
    const Vector3& startStep, const Vector3& endStep, double halfThickness,
    int readoutDirection, double lorentzAngle) const {
  Vector3 stepCenter = 0.5 * (startStep + endStep);
  // take the full drift length
  // this is the absolute drift in z
  double driftInZ = halfThickness - readoutDirection * stepCenter.z();
  // this is the absolute drift length
  double driftLength = driftInZ / cos(lorentzAngle);
  // project to parameter the readout surface
  double lorentzDeltaX = readoutDirection * driftInZ * tan(lorentzAngle);
  // the projected center, it has the lorentz shift applied
  Vector2 stepCenterProjected(stepCenter.x() + lorentzDeltaX, stepCenter.y());
  // the cell & its center
  Acts::DigitizationCell dCell = cell(stepCenterProjected);
  Vector2 cellCenter = cellPosition(dCell);
  // we are ready to return what we have
  return DigitizationStep((endStep - startStep).norm(), driftLength, dCell,
                          startStep, endStep, stepCenterProjected, cellCenter);
}
