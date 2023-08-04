// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlanarModuleStepper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Digitization/PlanarModuleStepper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Digitization/DigitizationCell.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <cstddef>
#include <ostream>

Acts::PlanarModuleStepper::PlanarModuleStepper(
    std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {}

std::vector<Acts::DigitizationStep> Acts::PlanarModuleStepper::cellSteps(
    const GeometryContext& gctx, const DigitizationModule& dmodule,
    const Vector3& startPoint, const Vector3& endPoint) const {
  // create the return vector
  std::vector<DigitizationStep> cSteps;

  // get the test surfaces for bin intersections
  auto stepSurfaces = dmodule.stepSurfaces(startPoint, endPoint);

  // the track direction
  Vector3 trackDirection((endPoint - startPoint).normalized());

  // the intersections through the surfaces, start one is the first valid one
  std::vector<Acts::Intersection3D> stepIntersections;
  stepIntersections.reserve(stepSurfaces.size() + 1);

  // run them - and check for the fast exit
  for (auto& sSurface : stepSurfaces) {
    // try it out by intersecting, but do not force the direction
    auto sIntersection =
        sSurface->intersect(gctx, startPoint, trackDirection, true);
    if (bool(sIntersection)) {
      // now record
      stepIntersections.push_back(sIntersection.intersection);
      ACTS_VERBOSE("Boundary Surface intersected with = "
                   << sIntersection.intersection.position.x() << ", "
                   << sIntersection.intersection.position.y() << ", "
                   << sIntersection.intersection.position.z());
    }
  }
  // Last one is also valid - now sort
  stepIntersections.push_back(
      Intersection3D(endPoint, (startPoint - endPoint).norm(),
                     Intersection3D::Status::reachable));
  std::sort(stepIntersections.begin(), stepIntersections.end());

  Vector3 lastPosition = startPoint;
  // reserve the right amount
  cSteps.reserve(stepIntersections.size());
  for (auto& sIntersection : stepIntersections) {
    // create the new digitization step
    cSteps.push_back(
        dmodule.digitizationStep(lastPosition, sIntersection.position));
    lastPosition = sIntersection.position;
  }
  // return all the steps
  return cSteps;
}

// calculate the steps caused by this track - fast simulation interface
std::vector<Acts::DigitizationStep> Acts::PlanarModuleStepper::cellSteps(
    const GeometryContext& gctx, const Acts::DigitizationModule& dmodule,
    const Vector2& moduleIntersection, const Vector3& trackDirection) const {
  // first, intersect the boundary surfaces
  auto boundarySurfaces = dmodule.boundarySurfaces();
  // intersect them - fast exit for cases where
  // readout and counter readout are hit
  Vector3 intersection3D(moduleIntersection.x(), moduleIntersection.y(), 0.);
  size_t attempts = 0;
  // the collected intersections
  std::vector<Acts::Intersection3D> boundaryIntersections;
  // run them - and check for the fast exit
  for (auto& bSurface : boundarySurfaces) {
    // count as an attempt
    ++attempts;
    // try it out by intersecting, but do not force the direction
    auto bIntersection =
        bSurface->intersect(gctx, intersection3D, trackDirection, true);
    if (bool(bIntersection)) {
      // now record
      boundaryIntersections.push_back(bIntersection.intersection);
      ACTS_VERBOSE("Boundary Surface intersected with = "
                   << bIntersection.intersection.position.x() << ", "
                   << bIntersection.intersection.position.y() << ", "
                   << bIntersection.intersection.position.z());
    }
    // fast break in case of readout/counter surface hit
    // the first two attempts are the module faces, if they are hit,
    // the stepper has run ok.
    if (attempts == 2 && boundaryIntersections.size() == attempts) {
      break;
    }
  }
  // Post-process if we have more than 2 intersections
  // only first or last can be wrong after resorting
  if (boundaryIntersections.size() > 2) {
    ACTS_VERBOSE(
        "More than 2 Boundary Surfaces intersected, this is an edge "
        "case, resolving ... ");
    std::sort(boundaryIntersections.begin(), boundaryIntersections.end());
  }
  // if for some reason the intersection does not work
  if (boundaryIntersections.empty()) {
    return std::vector<Acts::DigitizationStep>();
  }
  // return
  return cellSteps(gctx, dmodule, boundaryIntersections[0].position,
                   boundaryIntersections[1].position);
}
