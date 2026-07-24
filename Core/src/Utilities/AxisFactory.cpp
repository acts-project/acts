// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/AxisFactory.hpp"

#include <algorithm>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace Acts {

namespace {

void checkStrictlyIncreasing(const std::vector<double>& edges,
                             const std::string& context) {
  if (edges.size() < 2) {
    throw std::invalid_argument(context + ": at least two edges are required");
  }
  if (!std::ranges::is_sorted(edges, std::less_equal<double>())) {
    throw std::invalid_argument(context +
                                ": edges must be strictly increasing");
  }
}

}  // namespace

AxisFactory::AxisFactory(Variant variant,
                         std::optional<AxisDirection> direction)
    : m_variant(std::move(variant)), m_direction(direction) {}

AxisFactory AxisFactory::Equidistant(AxisBoundaryType boundaryType, double min,
                                     double max, std::size_t nBins,
                                     std::optional<AxisDirection> direction) {
  if (min >= max) {
    throw std::invalid_argument("AxisFactory::Equidistant: min must be < max");
  }
  if (nBins == 0) {
    throw std::invalid_argument(
        "AxisFactory::Equidistant: at least one bin is required");
  }
  return AxisFactory(EquidistantParams{boundaryType, min, max, nBins},
                     direction);
}

AxisFactory AxisFactory::Variable(AxisBoundaryType boundaryType,
                                  std::vector<double> edges,
                                  std::optional<AxisDirection> direction) {
  checkStrictlyIncreasing(edges, "AxisFactory::Variable");
  return AxisFactory(VariableParams{boundaryType, std::move(edges)}, direction);
}

AxisFactory AxisFactory::DeferredEquidistant(
    std::size_t nBins, std::optional<AxisDirection> direction) {
  if (nBins == 0) {
    throw std::invalid_argument(
        "AxisFactory::DeferredEquidistant: at least one bin is required");
  }
  return AxisFactory(DeferredEquidistantParams{nBins}, direction);
}

AxisFactory AxisFactory::DeferredVariable(
    std::vector<double> normalizedEdges,
    std::optional<AxisDirection> direction) {
  checkStrictlyIncreasing(normalizedEdges, "AxisFactory::DeferredVariable");
  if (normalizedEdges.front() != 0. || normalizedEdges.back() != 1.) {
    throw std::invalid_argument(
        "AxisFactory::DeferredVariable: edges must be normalized to [0, 1]");
  }
  return AxisFactory(DeferredVariableParams{std::move(normalizedEdges)},
                     direction);
}

AxisFactory AxisFactory::FromAxis(const IAxis& axis) {
  if (axis.getType() == AxisType::Equidistant) {
    return Equidistant(axis.getBoundaryType(), axis.getMin(), axis.getMax(),
                       axis.getNBins(), axis.getDirection());
  }
  return Variable(axis.getBoundaryType(), axis.getBinEdges(),
                  axis.getDirection());
}

AxisFactory AxisFactory::withDirection(AxisDirection direction) const {
  return AxisFactory(m_variant, direction);
}

AxisFactory AxisFactory::toDeferred() const {
  return std::visit(
      [this]<typename T>(const T& params) {
        if constexpr (std::is_same_v<T, EquidistantParams>) {
          return AxisFactory(DeferredEquidistantParams{params.nBins},
                             m_direction);
        } else if constexpr (std::is_same_v<T, VariableParams>) {
          std::vector<double> normalizedEdges = params.edges;
          double min = normalizedEdges.front();
          double max = normalizedEdges.back();
          for (double& edge : normalizedEdges) {
            edge = (edge - min) / (max - min);
          }
          // Force exact endpoints against floating point round-off
          normalizedEdges.front() = 0.;
          normalizedEdges.back() = 1.;
          return AxisFactory(DeferredVariableParams{std::move(normalizedEdges)},
                             m_direction);
        } else {
          return *this;
        }
      },
      m_variant);
}

bool AxisFactory::isDeferred() const {
  return std::holds_alternative<DeferredEquidistantParams>(m_variant) ||
         std::holds_alternative<DeferredVariableParams>(m_variant);
}

bool AxisFactory::isEquidistant() const {
  return std::holds_alternative<EquidistantParams>(m_variant) ||
         std::holds_alternative<DeferredEquidistantParams>(m_variant);
}

bool AxisFactory::isVariable() const {
  return !isEquidistant();
}

std::optional<AxisDirection> AxisFactory::direction() const {
  return m_direction;
}

std::optional<AxisBoundaryType> AxisFactory::boundaryType() const {
  return std::visit(
      []<typename T>(const T& params) -> std::optional<AxisBoundaryType> {
        if constexpr (std::is_same_v<T, EquidistantParams> ||
                      std::is_same_v<T, VariableParams>) {
          return params.boundaryType;
        } else {
          return std::nullopt;
        }
      },
      m_variant);
}

std::size_t AxisFactory::nBins() const {
  return std::visit(
      []<typename T>(const T& params) -> std::size_t {
        if constexpr (std::is_same_v<T, EquidistantParams> ||
                      std::is_same_v<T, DeferredEquidistantParams>) {
          return params.nBins;
        } else if constexpr (std::is_same_v<T, VariableParams>) {
          return params.edges.size() - 1;
        } else {
          return params.normalizedEdges.size() - 1;
        }
      },
      m_variant);
}

const AxisFactory::EquidistantParams& AxisFactory::asEquidistant() const {
  return std::get<EquidistantParams>(m_variant);
}

const AxisFactory::VariableParams& AxisFactory::asVariable() const {
  return std::get<VariableParams>(m_variant);
}

const AxisFactory::DeferredEquidistantParams&
AxisFactory::asDeferredEquidistant() const {
  return std::get<DeferredEquidistantParams>(m_variant);
}

const AxisFactory::DeferredVariableParams& AxisFactory::asDeferredVariable()
    const {
  return std::get<DeferredVariableParams>(m_variant);
}

std::optional<AxisDirection> AxisFactory::resolveDirection(
    std::optional<AxisDirection> direction) const {
  if (m_direction.has_value() && direction.has_value() &&
      m_direction != direction) {
    throw std::invalid_argument(
        "AxisFactory: axis direction " + axisDirectionName(*m_direction) +
        " does not match the direction " + axisDirectionName(*direction) +
        " expected by the caller");
  }
  return direction.has_value() ? direction : m_direction;
}

std::unique_ptr<IAxis> AxisFactory::toAxis(
    std::optional<AxisDirection> direction) const {
  std::optional<AxisDirection> effectiveDirection = resolveDirection(direction);
  return std::visit(
      [&]<typename T>(const T& params) -> std::unique_ptr<IAxis> {
        if constexpr (std::is_same_v<T, EquidistantParams>) {
          return IAxis::createEquidistant(params.boundaryType, params.min,
                                          params.max, params.nBins,
                                          effectiveDirection);
        } else if constexpr (std::is_same_v<T, VariableParams>) {
          return IAxis::createVariable(params.boundaryType, params.edges,
                                       effectiveDirection);
        } else {
          throw std::domain_error(
              "AxisFactory: a deferred axis description requires a range and "
              "boundary type to produce an axis");
        }
      },
      m_variant);
}

std::unique_ptr<IAxis> AxisFactory::toAxis(
    const AxisResolution& resolution,
    std::optional<AxisDirection> direction) const {
  if (resolution.min >= resolution.max) {
    throw std::invalid_argument("AxisFactory::toAxis: min must be < max");
  }
  std::optional<AxisDirection> effectiveDirection = resolveDirection(direction);
  return std::visit(
      [&]<typename T>(const T& params) -> std::unique_ptr<IAxis> {
        if constexpr (std::is_same_v<T, DeferredEquidistantParams>) {
          return IAxis::createEquidistant(resolution.boundaryType,
                                          resolution.min, resolution.max,
                                          params.nBins, effectiveDirection);
        } else if constexpr (std::is_same_v<T, DeferredVariableParams>) {
          std::vector<double> edges = params.normalizedEdges;
          for (double& edge : edges) {
            edge = resolution.min + edge * (resolution.max - resolution.min);
          }
          return IAxis::createVariable(resolution.boundaryType, edges,
                                       effectiveDirection);
        } else {
          throw std::domain_error(
              "AxisFactory: a fully specified axis description cannot be "
              "resolved with a range");
        }
      },
      m_variant);
}

std::string AxisFactory::toString() const {
  std::stringstream ss;
  ss << *this;
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, const AxisFactory& axisFactory) {
  os << "AxisFactory: " << axisFactory.nBins() << " bins";
  if (axisFactory.m_direction.has_value()) {
    os << " in " << axisDirectionName(*axisFactory.m_direction);
  }
  os << (axisFactory.isEquidistant() ? ", equidistant" : ", variable");
  std::visit(
      [&os]<typename T>(const T& params) {
        if constexpr (std::is_same_v<T, AxisFactory::EquidistantParams>) {
          os << " within [" << params.min << ", " << params.max << "], "
             << params.boundaryType;
        } else if constexpr (std::is_same_v<T, AxisFactory::VariableParams>) {
          os << " within [" << params.edges.front() << ", "
             << params.edges.back() << "], " << params.boundaryType;
        } else {
          os << " within deferred range";
        }
      },
      axisFactory.m_variant);
  return os;
}

}  // namespace Acts
