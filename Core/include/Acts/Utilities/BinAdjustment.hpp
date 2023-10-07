// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinAdjustment.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"

#include <stdexcept>

namespace Acts {

/// @brief adjust the BinUtility bu to the dimensions of radial bounds
///
/// @param bu BinUtility at source
/// @param rBounds the Radial bounds to adjust to
/// @param transform Transform for the adjusted @c BinUtility
///
/// @return new updated BinUtiltiy
BinUtility adjustBinUtility(const BinUtility& bu, const RadialBounds& rBounds,
                            const Transform3& transform) {
  // Default constructor
  BinUtility uBinUtil(transform);

  // The parameters from the cylinder bounds
  double minR = rBounds.get(RadialBounds::eMinR);
  double maxR = rBounds.get(RadialBounds::eMaxR);
  double minPhi = rBounds.get(RadialBounds::eAveragePhi) -
                  rBounds.get(RadialBounds::eHalfPhiSector);
  double maxPhi = rBounds.get(RadialBounds::eAveragePhi) +
                  rBounds.get(RadialBounds::eHalfPhiSector);
  // Retrieve the binning data
  const std::vector<BinningData>& bData = bu.binningData();
  // Loop over the binning data and adjust the dimensions
  for (auto& bd : bData) {
    // The binning value
    BinningValue bval = bd.binvalue;
    // Throw exceptions is stuff doesn't make sense:
    // - not the right binning value
    // - not equidistant
    if (bd.type == arbitrary) {
      throw std::invalid_argument("Arbirary binning can not be adjusted.");
    } else if (bval != binR and bval != binPhi) {
      throw std::invalid_argument("Disc binning must be: phi, r");
    }
    float min = 0., max = 0.;
    // Perform the value adjustment
    if (bval == binPhi) {
      min = minPhi;
      max = maxPhi;
    } else {
      min = minR;
      max = maxR;
    }
    // Create the updated BinningData
    BinningData uBinData(bd.option, bval, bd.bins(), min, max);
    uBinUtil += BinUtility(uBinData);
  }
  return uBinUtil;
}

/// @brief adjust the BinUtility bu to the dimensions of cylinder bounds
///
/// @param bu BinUtility at source
/// @param cBounds the Cylinder bounds to adjust to
/// @param transform Transform for the adjusted @c BinUtility
///
/// @return new updated BinUtiltiy
BinUtility adjustBinUtility(const BinUtility& bu, const CylinderBounds& cBounds,
                            const Transform3& transform) {
  // Default constructor
  BinUtility uBinUtil(transform);

  // The parameters from the cylinder bounds
  double cR = cBounds.get(CylinderBounds::eR);
  double cHz = cBounds.get(CylinderBounds::eHalfLengthZ);
  double avgPhi = cBounds.get(CylinderBounds::eAveragePhi);
  double halfPhi = cBounds.get(CylinderBounds::eHalfPhiSector);
  double minPhi = avgPhi - halfPhi;
  double maxPhi = avgPhi + halfPhi;

  // Retrieve the binning data
  const std::vector<BinningData>& bData = bu.binningData();
  // Loop over the binning data and adjust the dimensions
  for (auto& bd : bData) {
    // The binning value
    BinningValue bval = bd.binvalue;
    // Throw exceptions if stuff doesn't make sense:
    // - not the right binning value
    // - not equidistant
    if (bd.type == arbitrary) {
      throw std::invalid_argument("Arbitrary binning can not be adjusted.");
    } else if (bval != binRPhi and bval != binPhi and bval != binZ) {
      throw std::invalid_argument("Cylinder binning must be: rphi, phi, z");
    }
    float min = 0., max = 0.;
    // Perform the value adjustment
    if (bval == binPhi) {
      min = minPhi;
      max = maxPhi;
    } else if (bval == binRPhi) {
      min = cR * minPhi;
      max = cR * maxPhi;
    } else {
      min = -cHz;
      max = cHz;
    }
    // Create the updated BinningData
    BinningData uBinData(bd.option, bval, bd.bins(), min, max);
    uBinUtil += BinUtility(uBinData);
  }
  return uBinUtil;
}

/// @brief adjust the BinUtility bu to the dimensions of plane bounds
///
/// @param bu BinUtility at source
/// @param pBounds the Rectangle bounds to adjust to
/// @param transform Transform for the adjusted @c BinUtility
///
/// @return new updated BinUtiltiy
BinUtility adjustBinUtility(const BinUtility& bu,
                            const RectangleBounds& pBounds,
                            const Transform3& transform) {
  // Default constructor
  BinUtility uBinUtil(transform);

  // The parameters from the cylinder bounds
  double minX = pBounds.get(RectangleBounds::eMinX);
  double minY = pBounds.get(RectangleBounds::eMinY);
  double maxX = pBounds.get(RectangleBounds::eMaxX);
  double maxY = pBounds.get(RectangleBounds::eMaxY);

  // Retrieve the binning data
  const std::vector<BinningData>& bData = bu.binningData();
  // Loop over the binning data and adjust the dimensions
  for (auto& bd : bData) {
    // The binning value
    BinningValue bval = bd.binvalue;
    // Throw exceptions if stuff doesn't make sense:
    // - not the right binning value
    // - not equidistant
    if (bd.type == arbitrary) {
      throw std::invalid_argument("Arbitrary binning can not be adjusted.");
    } else if (bval != binX and bval != binY) {
      throw std::invalid_argument("Rectangle binning must be: x, y. ");
    }
    float min = 0., max = 0.;
    // Perform the value adjustment
    if (bval == binX) {
      min = minX;
      max = maxX;
    } else {
      min = minY;
      max = maxY;
    }
    // Create the updated BinningData
    BinningData uBinData(bd.option, bval, bd.bins(), min, max);
    uBinUtil += BinUtility(uBinData);
  }

  return uBinUtil;
}

/// @brief adjust the BinUtility bu to the dimensions of plane bounds
///
/// @param bu BinUtility at source
/// @param pBounds the Trapezoid bounds to adjust to
/// @param transform Transform for the adjusted @c BinUtility
///
/// @return new updated BinUtiltiy
BinUtility adjustBinUtility(const BinUtility& bu,
                            const TrapezoidBounds& pBounds,
                            const Transform3& transform) {
  // Default constructor
  BinUtility uBinUtil(transform);

  // The parameters from the cylinder bounds

  double halfX = std::max(pBounds.get(Acts::TrapezoidBounds::eHalfLengthXnegY),
                          pBounds.get(Acts::TrapezoidBounds::eHalfLengthXposY));
  double halfY = pBounds.get(Acts::TrapezoidBounds::eHalfLengthY);

  // Retrieve the binning data
  const std::vector<BinningData>& bData = bu.binningData();
  // Loop over the binning data and adjust the dimensions
  for (auto& bd : bData) {
    // The binning value
    BinningValue bval = bd.binvalue;
    // Throw exceptions if stuff doesn't make sense:
    // - not the right binning value
    // - not equidistant
    if (bd.type == arbitrary) {
      throw std::invalid_argument("Arbitrary binning can not be adjusted.");
    } else if (bval != binX and bval != binY) {
      throw std::invalid_argument("Rectangle binning must be: x, y. ");
    }
    float min = 0., max = 0.;
    // Perform the value adjustment
    if (bval == binX) {
      min = -1 * halfX;
      max = halfX;
    } else {
      min = -1 * halfY;
      max = halfY;
    }
    // Create the updated BinningData
    BinningData uBinData(bd.option, bval, bd.bins(), min, max);
    uBinUtil += BinUtility(uBinData);
  }

  return uBinUtil;
}

/// @brief adjust the BinUtility bu to a surface
///
/// @param bu BinUtility at source
/// @param surface Surface to which the adjustment is being done
/// @param gctx Geometry context to get the surfaces transform
///
/// @return new updated BinUtiltiy
BinUtility adjustBinUtility(const BinUtility& bu, const Surface& surface,
                            const GeometryContext& gctx) {
  // The surface type is a cylinder
  if (surface.type() == Surface::Cylinder) {
    // Cast to Cylinder bounds and return
    auto cBounds = dynamic_cast<const CylinderBounds*>(&(surface.bounds()));
    // Return specific adjustment
    return adjustBinUtility(bu, *cBounds, surface.transform(gctx));

  } else if (surface.type() == Surface::Disc) {
    // Cast to Cylinder bounds and return
    auto rBounds = dynamic_cast<const RadialBounds*>(&(surface.bounds()));
    // Return specific adjustment
    return adjustBinUtility(bu, *rBounds, surface.transform(gctx));
  } else if (surface.type() == Surface::Plane) {
    if (surface.bounds().type() == SurfaceBounds::eRectangle) {
      // Cast to Plane bounds and return
      auto pBounds = dynamic_cast<const RectangleBounds*>(&(surface.bounds()));
      // Return specific adjustment
      return adjustBinUtility(bu, *pBounds, surface.transform(gctx));
    } else if (surface.bounds().type() == SurfaceBounds::eTrapezoid) {
      // Cast to Plane bounds and return
      auto pBounds = dynamic_cast<const TrapezoidBounds*>(&(surface.bounds()));
      // Return specific adjustment
      return adjustBinUtility(bu, *pBounds, surface.transform(gctx));
    } else {
      throw std::invalid_argument(
          "Bin adjustment not implemented for this type of plane surface yet!");
    }
  }

  throw std::invalid_argument(
      "Bin adjustment not implemented for this surface yet!");

  return BinUtility();
}

}  // namespace Acts
