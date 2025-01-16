// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/ProtoAxis.hpp"

#include <algorithm>

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, Acts::AxisBoundaryType abType,
                           const std::vector<double>& edges)
    : m_axisDir(aDir),
      m_axisType(AxisType::Variable),
      m_axisBoundaryType(abType),
      m_edges(edges),
      m_autorange(false) {

  if (m_edges.size() < 2) {
    throw std::invalid_argument(
        "ProtoBinning: Invalid binning, at least two bin edges are needed.");
  }

  if (!std::is_sorted(m_edges.begin(), m_edges.end())) {
    throw std::invalid_argument(
        "ProtoBinning: Invalid binning, bin edges are not sorted.");
  }

  if (abType == AxisBoundaryType::Closed && aDir != AxisDirection::AxisPhi &&
      aDir != AxisDirection::AxisRPhi) {
    std::string msg =
        "ProtoBinning: Invalid axis boundary type 'Closed' for direction '";
    msg += axisDirectionName(aDir) + "'.";
    throw std::invalid_argument(msg);
  }
}

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, AxisBoundaryType abType,
                           double minE, double maxE, std::size_t nbins)
    : m_axisDir(aDir), m_axisBoundaryType(abType), m_autorange(false) {
  if (minE >= maxE) {
    std::string msg = "ProtoBinning: Invalid axis range for direction '";
    msg += axisDirectionName(m_axisDir);
    msg += "', min edge (" + std::to_string(minE) + ") ";
    msg += " needs to be smaller than max edge (";
    msg += std::to_string(maxE) + ").";
    throw std::invalid_argument(msg);
  }
  if (nbins < 1u) {
    throw std::invalid_argument(
        "ProtoBinning: Invalid binning, at least one bin is needed.");
  }

  if (abType == AxisBoundaryType::Closed && aDir != AxisDirection::AxisPhi &&
      aDir != AxisDirection::AxisRPhi) {
    std::string msg =
        "ProtoBinning: Invalid axis boundary type 'Closed' for direction '";
    msg += axisDirectionName(aDir) + "'.";
    throw std::invalid_argument(msg);
  }

  double stepE = (maxE - minE) / nbins;
  m_edges.reserve(nbins + 1);
  for (std::size_t i = 0; i <= nbins; i++) {
    m_edges.push_back(minE + i * stepE);
  }
}

Acts::ProtoAxis::ProtoAxis(AxisDirection aDir, AxisBoundaryType abType,
                           std::size_t nbins)
    : m_axisDir(aDir),
      m_axisBoundaryType(abType),
      m_edges(nbins + 1, 0.),
      m_autorange(true) {
  if (abType == AxisBoundaryType::Closed && aDir != AxisDirection::AxisPhi &&
      aDir != AxisDirection::AxisRPhi) {
    std::string msg =
        "ProtoBinning: Invalid axis boundary type 'Closed' for direction '";
    msg += axisDirectionName(aDir) + "'.";
    throw std::invalid_argument(msg);
  }
}

bool Acts::ProtoAxis::isEquidistant() const {
  return m_axisType == AxisType::Equidistant;
}

bool Acts::ProtoAxis::isVariable() const {
  return m_axisType == AxisType::Variable;
}

bool Acts::ProtoAxis::isAutorange() const {
  return m_autorange;
}

Acts::AxisDirection Acts::ProtoAxis::getAxisDirection() const {
  return m_axisDir;
}

Acts::AxisType Acts::ProtoAxis::getType() const {
  return m_axisType;
}

Acts::AxisBoundaryType Acts::ProtoAxis::getBoundaryType() const {
  return m_axisBoundaryType;
}

std::vector<double> Acts::ProtoAxis::getBinEdges() const {
  return m_edges;
}

double Acts::ProtoAxis::getMin() const {
  return m_edges.front();
}

double Acts::ProtoAxis::getMax() const {
  return m_edges.back();
}

std::size_t Acts::ProtoAxis::getNBins() const {
  return m_edges.size() - 1;
};

void Acts::ProtoAxis::toStream(std::ostream& os) const {
  os << toString();
}

std::string Acts::ProtoAxis::toString() const {
  std::stringstream ss;
  ss << "ProtoAxis: " << getNBins() << " bins in "
     << axisDirectionName(m_axisDir);
  ss << (m_axisType == AxisType::Variable ? ", variable " : ", equidistant ");
  if (!m_autorange) {
    ss << "within [" << m_edges.front() << ", " << m_edges.back() << "] ";
  } else {
    ss << "within automatic range";
  }
  return ss.str();
}
