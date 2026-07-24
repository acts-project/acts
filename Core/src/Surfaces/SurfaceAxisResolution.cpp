// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceAxisResolution.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace {

using namespace Acts;

std::string unsupportedSurfaceMessage(const Surface& surface) {
  std::stringstream ss;
  ss << "Axis resolution is not supported for this surface:\n"
     << surface.toStream(GeometryContext::dangerouslyDefaultConstruct());
  return ss.str();
}

AxisResolution resolveCylinder(const CylinderBounds& cBounds,
                               AxisDirection aDir) {
  double r = cBounds.get(CylinderBounds::eR);
  double hz = cBounds.get(CylinderBounds::eHalfLengthZ);
  double avgPhi = cBounds.get(CylinderBounds::eAveragePhi);
  double halfPhi = cBounds.get(CylinderBounds::eHalfPhiSector);
  AxisBoundaryType phiBoundaryType = cBounds.coversFullAzimuth()
                                         ? AxisBoundaryType::Closed
                                         : AxisBoundaryType::Bound;
  switch (aDir) {
    case AxisDirection::AxisRPhi:
      return {r * (avgPhi - halfPhi), r * (avgPhi + halfPhi), phiBoundaryType};
    case AxisDirection::AxisPhi:
      return {avgPhi - halfPhi, avgPhi + halfPhi, phiBoundaryType};
    case AxisDirection::AxisZ:
      return {-hz, hz, AxisBoundaryType::Bound};
    default:
      throw std::invalid_argument(
          "Cylinder axis resolution must be along rphi, phi or z");
  }
}

AxisResolution resolveDisc(const RadialBounds& rBounds, AxisDirection aDir) {
  double avgPhi = rBounds.get(RadialBounds::eAveragePhi);
  double halfPhi = rBounds.get(RadialBounds::eHalfPhiSector);
  switch (aDir) {
    case AxisDirection::AxisR:
      return {rBounds.get(RadialBounds::eMinR),
              rBounds.get(RadialBounds::eMaxR), AxisBoundaryType::Bound};
    case AxisDirection::AxisPhi:
      return {avgPhi - halfPhi, avgPhi + halfPhi,
              rBounds.coversFullAzimuth() ? AxisBoundaryType::Closed
                                          : AxisBoundaryType::Bound};
    default:
      throw std::invalid_argument(
          "Disc axis resolution must be along r or phi");
  }
}

AxisResolution resolveRectangle(const RectangleBounds& pBounds,
                                AxisDirection aDir) {
  switch (aDir) {
    case AxisDirection::AxisX:
      return {pBounds.get(RectangleBounds::eMinX),
              pBounds.get(RectangleBounds::eMaxX), AxisBoundaryType::Bound};
    case AxisDirection::AxisY:
      return {pBounds.get(RectangleBounds::eMinY),
              pBounds.get(RectangleBounds::eMaxY), AxisBoundaryType::Bound};
    default:
      throw std::invalid_argument(
          "Rectangle axis resolution must be along x or y");
  }
}

AxisResolution resolveTrapezoid(const TrapezoidBounds& pBounds,
                                AxisDirection aDir) {
  switch (aDir) {
    case AxisDirection::AxisX: {
      double halfX = std::max(pBounds.get(TrapezoidBounds::eHalfLengthXnegY),
                              pBounds.get(TrapezoidBounds::eHalfLengthXposY));
      return {-halfX, halfX, AxisBoundaryType::Bound};
    }
    case AxisDirection::AxisY: {
      double halfY = pBounds.get(TrapezoidBounds::eHalfLengthY);
      return {-halfY, halfY, AxisBoundaryType::Bound};
    }
    default:
      throw std::invalid_argument(
          "Trapezoid axis resolution must be along x or y");
  }
}

}  // namespace

Acts::AxisResolution Acts::surfaceAxisResolution(const Surface& surface,
                                                 AxisDirection aDir) {
  if (auto b = dynamic_cast<const CylinderBounds*>(&surface.bounds());
      b != nullptr) {
    return resolveCylinder(*b, aDir);
  }
  if (auto b = dynamic_cast<const RadialBounds*>(&surface.bounds());
      b != nullptr) {
    return resolveDisc(*b, aDir);
  }
  if (surface.type() == Surface::Plane) {
    if (auto b = dynamic_cast<const RectangleBounds*>(&surface.bounds());
        b != nullptr) {
      return resolveRectangle(*b, aDir);
    }
    if (auto b = dynamic_cast<const TrapezoidBounds*>(&surface.bounds());
        b != nullptr) {
      return resolveTrapezoid(*b, aDir);
    }
  }
  throw std::invalid_argument(unsupportedSurfaceMessage(surface));
}

std::vector<std::unique_ptr<Acts::IAxis>> Acts::resolveAxes(
    const MultiAxisFactory& multiAxisFactory, const Surface& surface) {
  std::array<AxisDirection, 2> canonical = surface.localAxes();
  std::size_t n = multiAxisFactory.size();
  if (n > canonical.size()) {
    throw std::invalid_argument(
        "resolveAxes: binning has more dimensions than the surface supports");
  }

  std::size_t nDirected = std::ranges::count_if(
      multiAxisFactory.axisFactories(),
      [](const AxisFactory& af) { return af.direction().has_value(); });
  if (nDirected != 0 && nDirected != n) {
    throw std::invalid_argument(
        "resolveAxes: either all or none of the axes must carry a direction");
  }

  // Determine per output slot the canonical direction and the description
  // that binds to it. Without stored directions the match is positional;
  // with stored directions the descriptions are re-ordered to canonical
  // direction order.
  std::vector<std::pair<AxisDirection, const AxisFactory*>> slots;
  slots.reserve(n);
  if (nDirected == 0) {
    for (std::size_t i = 0; i < n; ++i) {
      slots.emplace_back(canonical[i], &multiAxisFactory.axisFactory(i));
    }
  } else {
    for (AxisDirection aDir : canonical) {
      auto it = std::ranges::find_if(
          multiAxisFactory.axisFactories(),
          [aDir](const AxisFactory& af) { return af.direction() == aDir; });
      if (it != multiAxisFactory.axisFactories().end()) {
        slots.emplace_back(aDir, &(*it));
      }
    }
    if (slots.size() != n) {
      std::stringstream ss;
      ss << "resolveAxes: binning directions do not match the surface axes (";
      for (std::size_t i = 0; i < canonical.size(); ++i) {
        ss << (i > 0 ? ", " : "") << axisDirectionName(canonical[i]);
      }
      ss << ")";
      throw std::invalid_argument(ss.str());
    }
  }

  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(n);
  for (const auto& [aDir, axisFactory] : slots) {
    axes.push_back(
        axisFactory->toAxis(surfaceAxisResolution(surface, aDir), aDir));
  }
  return axes;
}
