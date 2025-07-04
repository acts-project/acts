// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/ProtoAxis.hpp"

#include <sstream>

Acts::ProtoAxis::ProtoAxis(Acts::AxisBoundaryType abType,
                           const std::vector<double>& edges)
    : m_axis(IAxis::createVariable(abType, edges)) {}

Acts::ProtoAxis::ProtoAxis(AxisBoundaryType abType, double minE, double maxE,
                           std::size_t nbins)
    : m_axis(IAxis::createEquidistant(abType, minE, maxE, nbins)) {}

Acts::ProtoAxis::ProtoAxis(AxisBoundaryType abType, std::size_t nbins)
    : m_axis(IAxis::createEquidistant(abType, 0., 1., nbins)),
      m_autorange(true) {}

Acts::ProtoAxis::ProtoAxis(const ProtoAxis& other)
    : m_autorange(other.m_autorange) {
  const auto& axis = other.getAxis();
  if (!m_autorange) {
    const auto& edges = axis.getBinEdges();
    if (axis.getType() == AxisType::Variable) {
      m_axis = IAxis::createVariable(axis.getBoundaryType(), edges);
    } else {
      m_axis = IAxis::createEquidistant(axis.getBoundaryType(), edges.front(),
                                        edges.back(), axis.getNBins());
    }
  } else {
    m_axis = IAxis::createEquidistant(axis.getBoundaryType(), 0., 1.,
                                      axis.getNBins());
  }
}

Acts::ProtoAxis& Acts::ProtoAxis::operator=(const ProtoAxis& other) {
  if (this != &other) {
    m_autorange = other.m_autorange;
    const auto& axis = other.getAxis();
    if (!m_autorange) {
      const auto& edges = axis.getBinEdges();
      if (axis.getType() == AxisType::Variable) {
        m_axis = IAxis::createVariable(axis.getBoundaryType(), edges);
      } else {
        m_axis = IAxis::createEquidistant(axis.getBoundaryType(), edges.front(),
                                          edges.back(), axis.getNBins());
      }
    } else {
      m_axis = IAxis::createEquidistant(axis.getBoundaryType(), 0., 1.,
                                        axis.getNBins());
    }
  }
  return *this;
}

const Acts::IAxis& Acts::ProtoAxis::getAxis() const {
  return *m_axis;
}

void Acts::ProtoAxis::setRange(double minE, double maxE) {
  if (minE > maxE) {
    throw std::invalid_argument(
        "ProtoAxis::setRange: minE > maxE is not allowed.");
  }

  if (m_axis->getType() == AxisType::Equidistant) {
    m_axis = IAxis::createEquidistant(m_axis->getBoundaryType(), minE, maxE,
                                      m_axis->getNBins());
  } else {
    std::vector<double> edges = m_axis->getBinEdges();
    // Clip it to min/max
    std::erase_if(edges,
                  [minE, maxE](double e) { return (e < minE || e > maxE); });
    // Add the min and max
    edges.emplace_back(minE);
    edges.emplace_back(maxE);
    std::ranges::sort(edges);
    m_axis = IAxis::createVariable(m_axis->getBoundaryType(), edges);
  }

  // Force autorange to be false
  m_autorange = false;
}

bool Acts::ProtoAxis::isAutorange() const {
  return m_autorange;
}

void Acts::ProtoAxis::toStream(std::ostream& os) const {
  os << toString();
}

std::string Acts::ProtoAxis::toString() const {
  std::stringstream ss;
  const auto& axis = getAxis();
  ss << "ProtoAxis: " << axis.getNBins() << " bins";
  ss << (axis.getType() == AxisType::Variable ? ", variable "
                                              : ", equidistant ");
  if (!m_autorange) {
    const auto& edges = axis.getBinEdges();
    ss << "within [" << edges.front() << ", " << edges.back() << "]";
  } else {
    ss << "within automatic range";
  }
  return ss.str();
}

// Ostream operator implementation
std::ostream& Acts::operator<<(std::ostream& os,
                               const std::vector<ProtoAxis>& as) {
  for (const auto& a : as) {
    os << a.toString() << '\n';
  }
  return os;
}

Acts::DirectedProtoAxis::DirectedProtoAxis(AxisDirection axisDir,
                                           AxisBoundaryType abType,
                                           const std::vector<double>& edges)
    : ProtoAxis(abType, edges), m_direction(axisDir) {}

Acts::DirectedProtoAxis::DirectedProtoAxis(AxisDirection axisDir,
                                           AxisBoundaryType abType, double minE,
                                           double maxE, std::size_t nbins)
    : ProtoAxis(abType, minE, maxE, nbins), m_direction(axisDir) {}

Acts::DirectedProtoAxis::DirectedProtoAxis(AxisDirection axisDir,
                                           AxisBoundaryType abType,
                                           std::size_t nbins)
    : ProtoAxis(abType, nbins), m_direction(axisDir) {}

Acts::AxisDirection Acts::DirectedProtoAxis::getAxisDirection() const {
  return m_direction;
}

void Acts::DirectedProtoAxis::toStream(std::ostream& os) const {
  os << toString();
}

std::string Acts::DirectedProtoAxis::toString() const {
  std::stringstream ss;
  const auto& axis = getAxis();
  ss << "DirectedProtoAxis: " << axis.getNBins() << " bins in "
     << axisDirectionName(m_direction);
  ss << (axis.getType() == AxisType::Variable ? ", variable "
                                              : ", equidistant ");
  if (!m_autorange) {
    const auto& edges = axis.getBinEdges();
    ss << "within [" << edges.front() << ", " << edges.back() << "]";
  } else {
    ss << "within automatic range";
  }
  return ss.str();
}

// Ostream operator implementation vector of directed proto axes
std::ostream& Acts::operator<<(std::ostream& os,
                               const std::vector<DirectedProtoAxis>& as) {
  for (const auto& a : as) {
    os << a << '\n';
  }
  return os;
}
