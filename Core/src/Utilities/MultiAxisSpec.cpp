// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/MultiAxisSpec.hpp"

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

AxisSpec::Options resolveCylinder(const CylinderBounds& cBounds,
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

AxisSpec::Options resolveDisc(const RadialBounds& rBounds, AxisDirection aDir) {
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

AxisSpec::Options resolveRectangle(const RectangleBounds& pBounds,
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

AxisSpec::Options resolveTrapezoid(const TrapezoidBounds& pBounds,
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

AxisSpec::Options resolveAxisOptions(const Surface& surface,
                                     AxisDirection aDir) {
  const SurfaceBounds& bounds = surface.bounds();
  AxisSpec::Options options = [&]() -> AxisSpec::Options {
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

MultiAxisSpec::MultiAxisSpec(std::vector<AxisSpec> axisSpecs)
    : m_axisSpecs(std::move(axisSpecs)) {
  if (m_axisSpecs.empty()) {
    throw std::invalid_argument(
        "MultiAxisSpec: at least one axis spec is required");
  }
}

std::size_t MultiAxisSpec::size() const {
  return m_axisSpecs.size();
}

const AxisSpec& MultiAxisSpec::axisSpec(std::size_t i) const {
  return m_axisSpecs.at(i);
}

std::span<const AxisSpec> MultiAxisSpec::axisSpecs() const {
  return m_axisSpecs;
}

bool MultiAxisSpec::isDeferred() const {
  return std::ranges::any_of(
      m_axisSpecs, [](const AxisSpec& af) { return af.isDeferred(); });
}

std::unique_ptr<IMultiAxis> MultiAxisSpec::buildMultiAxis(
    const Options& options) const {
  return buildMultiAxisImpl(options);
}

std::unique_ptr<IMultiAxis> MultiAxisSpec::buildMultiAxisImpl(
    std::span<const AxisSpec::Options> options) const {
  if (!options.empty() && options.size() != size()) {
    throw std::invalid_argument(
        "MultiAxisSpec: either one option set per axis or none at all");
  }

  std::vector<std::unique_ptr<IAxis>> axes;
  axes.reserve(size());
  for (std::size_t i = 0; i < size(); ++i) {
    axes.push_back(m_axisSpecs[i].buildAxis(
        options.empty() ? AxisSpec::Options{} : options[i]));
  }

  switch (axes.size()) {
    case 1:
      return IMultiAxis::create(*axes[0]);
    case 2:
      return IMultiAxis::create(*axes[0], *axes[1]);
    case 3:
      return IMultiAxis::create(*axes[0], *axes[1], *axes[2]);
    default:
      throw std::domain_error(
          "MultiAxisSpec: multi-axes support at most 3 dimensions");
  }
}

std::string MultiAxisSpec::toString() const {
  std::stringstream ss;
  ss << "MultiAxisSpec: " << size() << " axes [";
  for (std::size_t i = 0; i < size(); ++i) {
    ss << (i > 0 ? "; " : "") << axisSpec(i);
  }
  ss << "]";
  return ss.str();
}

std::unique_ptr<IMultiAxis2D> resolveMultiAxis(
    const MultiAxisSpec2D& multiAxisSpec, const Surface& surface) {
  std::span<const AxisSpec> axisSpecs = multiAxisSpec.axisSpecs();
  std::array<AxisDirection, 2> canonical = surface.localAxes();

  const std::size_t nDirected = std::ranges::count_if(
      axisSpecs, [](const AxisSpec& af) { return af.direction().has_value(); });
  if (nDirected == 1) {
    throw std::invalid_argument(
        "resolveMultiAxis: either both or neither of the axes must carry a "
        "direction");
  }

  // Bind a canonical direction to its spec, positionally without stored
  // directions and by direction otherwise
  auto slot = [&](std::size_t i) -> const AxisSpec& {
    if (nDirected == 0) {
      return axisSpecs[i];
    }
    AxisDirection aDir = canonical[i];
    auto it = std::ranges::find_if(axisSpecs, [aDir](const AxisSpec& af) {
      return af.direction() == aDir;
    });
    if (it == axisSpecs.end()) {
      std::stringstream ss;
      ss << "resolveMultiAxis: binning directions do not match the surface "
            "axes ("
         << axisDirectionName(canonical[0]) << ", "
         << axisDirectionName(canonical[1]) << ")";
      throw std::invalid_argument(ss.str());
    }
    return *it;
  };

  return IMultiAxis::create(
      *slot(0).buildAxis(resolveAxisOptions(surface, canonical[0])),
      *slot(1).buildAxis(resolveAxisOptions(surface, canonical[1])));
}

}  // namespace Acts
