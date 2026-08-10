// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MultiAxisFactory.hpp"

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <algorithm>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace Acts {

namespace {

std::string unsupportedSurfaceMessage(const Surface& surface) {
  std::stringstream ss;
  ss << "Axis resolution is not supported for this surface:\n"
     << surface.toStream(GeometryContext::dangerouslyDefaultConstruct());
  return ss.str();
}

AxisFactory::Options resolveCylinder(const CylinderBounds& cBounds,
                                     AxisDirection aDir) {
  using enum AxisDirection;

  double r = cBounds.get(CylinderBounds::eR);
  double hz = cBounds.get(CylinderBounds::eHalfLengthZ);
  double avgPhi = cBounds.get(CylinderBounds::eAveragePhi);
  double halfPhi = cBounds.get(CylinderBounds::eHalfPhiSector);
  AxisBoundaryType phiBoundaryType = cBounds.coversFullAzimuth()
                                         ? AxisBoundaryType::Closed
                                         : AxisBoundaryType::Bound;
  switch (aDir) {
    case AxisRPhi:
      return {.min = r * (avgPhi - halfPhi),
              .max = r * (avgPhi + halfPhi),
              .boundaryType = phiBoundaryType};
    case AxisPhi:
      return {.min = avgPhi - halfPhi,
              .max = avgPhi + halfPhi,
              .boundaryType = phiBoundaryType};
    case AxisZ:
      return {.min = -hz, .max = hz, .boundaryType = AxisBoundaryType::Bound};
    default:
      throw std::invalid_argument(
          "Cylinder axis resolution must be along rphi, phi or z");
  }
}

AxisFactory::Options resolveDisc(const RadialBounds& rBounds,
                                 AxisDirection aDir) {
  using enum AxisBoundaryType;

  double avgPhi = rBounds.get(RadialBounds::eAveragePhi);
  double halfPhi = rBounds.get(RadialBounds::eHalfPhiSector);
  switch (aDir) {
    case AxisDirection::AxisR:
      return {.min = rBounds.get(RadialBounds::eMinR),
              .max = rBounds.get(RadialBounds::eMaxR),
              .boundaryType = Bound};
    case AxisDirection::AxisPhi:
      return {.min = avgPhi - halfPhi,
              .max = avgPhi + halfPhi,
              .boundaryType = rBounds.coversFullAzimuth() ? Closed : Bound};
    default:
      throw std::invalid_argument(
          "Disc axis resolution must be along r or phi");
  }
}

AxisFactory::Options resolveRectangle(const RectangleBounds& pBounds,
                                      AxisDirection aDir) {
  using enum AxisBoundaryType;

  switch (aDir) {
    case AxisDirection::AxisX:
      return {.min = pBounds.get(RectangleBounds::eMinX),
              .max = pBounds.get(RectangleBounds::eMaxX),
              .boundaryType = Bound};
    case AxisDirection::AxisY:
      return {.min = pBounds.get(RectangleBounds::eMinY),
              .max = pBounds.get(RectangleBounds::eMaxY),
              .boundaryType = Bound};
    default:
      throw std::invalid_argument(
          "Rectangle axis resolution must be along x or y");
  }
}

AxisFactory::Options resolveTrapezoid(const TrapezoidBounds& pBounds,
                                      AxisDirection aDir) {
  using enum AxisBoundaryType;

  switch (aDir) {
    case AxisDirection::AxisX: {
      double halfX = std::max(pBounds.get(TrapezoidBounds::eHalfLengthXnegY),
                              pBounds.get(TrapezoidBounds::eHalfLengthXposY));
      return {.min = -halfX, .max = halfX, .boundaryType = Bound};
    }
    case AxisDirection::AxisY: {
      double halfY = pBounds.get(TrapezoidBounds::eHalfLengthY);
      return {.min = -halfY, .max = halfY, .boundaryType = Bound};
    }
    default:
      throw std::invalid_argument(
          "Trapezoid axis resolution must be along x or y");
  }
}

AxisFactory::Options resolveAxisOptions(const Surface& surface,
                                        AxisDirection aDir) {
  const SurfaceBounds& bounds = surface.bounds();
  AxisFactory::Options options = [&]() -> AxisFactory::Options {
    switch (bounds.type()) {
      case SurfaceBounds::eCylinder:
        return resolveCylinder(static_cast<const CylinderBounds&>(bounds),
                               aDir);
      case SurfaceBounds::eDisc:
        return resolveDisc(static_cast<const RadialBounds&>(bounds), aDir);
      case SurfaceBounds::eRectangle:
        return resolveRectangle(static_cast<const RectangleBounds&>(bounds),
                                aDir);
      case SurfaceBounds::eTrapezoid:
        return resolveTrapezoid(static_cast<const TrapezoidBounds&>(bounds),
                                aDir);
      default:
        throw std::invalid_argument(unsupportedSurfaceMessage(surface));
    }
  }();
  options.direction = aDir;
  return options;
}

}  // namespace

MultiAxisFactory::MultiAxisFactory(std::vector<AxisFactory> axisFactories)
    : m_axisFactories(std::move(axisFactories)) {
  if (m_axisFactories.empty()) {
    throw std::invalid_argument(
        "MultiAxisFactory: at least one axis description is required");
  }
}

std::size_t MultiAxisFactory::size() const {
  return m_axisFactories.size();
}

const AxisFactory& MultiAxisFactory::axisFactory(std::size_t i) const {
  return m_axisFactories.at(i);
}

const std::vector<AxisFactory>& MultiAxisFactory::axisFactories() const {
  return m_axisFactories;
}

bool MultiAxisFactory::isDeferred() const {
  return std::ranges::any_of(
      m_axisFactories, [](const AxisFactory& af) { return af.isDeferred(); });
}

std::vector<std::unique_ptr<IAxis>> MultiAxisFactory::buildAxes() const {
  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(size());
  for (const AxisFactory& af : m_axisFactories) {
    axes.push_back(af.buildAxis());
  }
  return axes;
}

std::vector<std::unique_ptr<IAxis>> MultiAxisFactory::buildAxes(
    OptionsView options) const {
  if (options.size() != size()) {
    throw std::invalid_argument(
        "MultiAxisFactory: one option set per axis is required");
  }
  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(size());
  for (std::size_t i = 0; i < size(); ++i) {
    axes.push_back(m_axisFactories[i].buildAxis(options[i]));
  }
  return axes;
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::buildMultiAxis() const {
  return makeMultiAxis(buildAxes());
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::buildMultiAxis(
    OptionsView options) const {
  return makeMultiAxis(buildAxes(options));
}

std::unique_ptr<IMultiAxis> MultiAxisFactory::makeMultiAxis(
    const std::vector<std::unique_ptr<IAxis>>& axes) {
  switch (axes.size()) {
    case 1:
      return IMultiAxis::create(*axes[0]);
    case 2:
      return IMultiAxis::create(*axes[0], *axes[1]);
    case 3:
      return IMultiAxis::create(*axes[0], *axes[1], *axes[2]);
    default:
      throw std::domain_error(
          "MultiAxisFactory: multi-axes support at most 3 dimensions");
  }
}

std::string MultiAxisFactory::toString() const {
  std::stringstream ss;
  ss << "MultiAxisFactory: " << size() << " axes [";
  for (std::size_t i = 0; i < size(); ++i) {
    ss << (i > 0 ? "; " : "") << axisFactory(i);
  }
  ss << "]";
  return ss.str();
}

std::unique_ptr<IMultiAxis2D> resolveMultiAxis(
    const MultiAxisFactory2D& multiAxisFactory, const Surface& surface) {
  const std::vector<AxisFactory>& axisFactories =
      multiAxisFactory.axisFactories();
  std::array<AxisDirection, 2> canonical = surface.localAxes();

  const std::size_t nDirected = std::ranges::count_if(
      axisFactories,
      [](const AxisFactory& af) { return af.direction().has_value(); });
  if (nDirected == 1) {
    throw std::invalid_argument(
        "resolveMultiAxis: either both or neither of the axes must carry a "
        "direction");
  }

  // Bind each canonical direction to its description, positionally without
  // stored directions and by direction otherwise
  std::array<const AxisFactory*, 2> slots{};
  for (std::size_t i = 0; i < canonical.size(); ++i) {
    if (nDirected == 0) {
      slots[i] = &axisFactories[i];
      continue;
    }
    AxisDirection aDir = canonical[i];
    auto it = std::ranges::find_if(
        axisFactories,
        [aDir](const AxisFactory& af) { return af.direction() == aDir; });
    if (it == axisFactories.end()) {
      std::stringstream ss;
      ss << "resolveMultiAxis: binning directions do not match the surface "
            "axes ("
         << axisDirectionName(canonical[0]) << ", "
         << axisDirectionName(canonical[1]) << ")";
      throw std::invalid_argument(ss.str());
    }
    slots[i] = std::to_address(it);
  }

  return IMultiAxis::create(
      *slots[0]->buildAxis(resolveAxisOptions(surface, canonical[0])),
      *slots[1]->buildAxis(resolveAxisOptions(surface, canonical[1])));
}

}  // namespace Acts
