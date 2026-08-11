// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/AxisSpec.hpp"

#include <algorithm>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace Acts {

namespace {

/// Take the value from whichever side gives it, requiring both to agree if
/// both do
template <typename T>
std::optional<T> mergeProperty(const std::optional<T>& specified,
                               const std::optional<T>& supplied,
                               const std::string& what) {
  if (specified.has_value() && supplied.has_value() &&
      *specified != *supplied) {
    throw std::invalid_argument("AxisSpec: the specified axis " + what +
                                " does not match the one supplied");
  }
  return specified.has_value() ? specified : supplied;
}

/// Unwrap a property that has to be known by now
template <typename T>
const T& requireProperty(const std::optional<T>& value,
                         const std::string& what) {
  if (!value.has_value()) {
    throw std::domain_error("AxisSpec: the axis " + what +
                            " is neither specified nor supplied");
  }
  return *value;
}

/// Range of a spec, unset where it is left to the consumer
std::optional<double> minOf(const AxisSpec::EquidistantParams& params) {
  return params.min;
}
std::optional<double> maxOf(const AxisSpec::EquidistantParams& params) {
  return params.max;
}
std::optional<double> minOf(const AxisSpec::VariableParams& params) {
  return params.edges.front();
}
std::optional<double> maxOf(const AxisSpec::VariableParams& params) {
  return params.edges.back();
}
std::optional<double> minOf(
    const AxisSpec::DeferredVariableParams& /*params*/) {
  return std::nullopt;
}
std::optional<double> maxOf(
    const AxisSpec::DeferredVariableParams& /*params*/) {
  return std::nullopt;
}

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

AxisSpec::AxisSpec(Variant variant,
                   std::optional<AxisBoundaryType> boundaryType,
                   std::optional<AxisDirection> direction)
    : m_variant(std::move(variant)),
      m_boundaryType(boundaryType),
      m_direction(direction) {}

AxisSpec AxisSpec::Equidistant(std::size_t nBins, std::optional<double> min,
                               std::optional<double> max,
                               std::optional<AxisBoundaryType> boundaryType,
                               std::optional<AxisDirection> direction) {
  if (min.has_value() && max.has_value() && *min >= *max) {
    throw std::invalid_argument("AxisSpec::Equidistant: min must be < max");
  }
  if (nBins == 0) {
    throw std::invalid_argument(
        "AxisSpec::Equidistant: at least one bin is required");
  }
  return AxisSpec(EquidistantParams{nBins, min, max}, boundaryType, direction);
}

AxisSpec AxisSpec::DeferredEquidistant(std::size_t nBins,
                                       std::optional<AxisDirection> direction) {
  return Equidistant(nBins, std::nullopt, std::nullopt, std::nullopt,
                     direction);
}

AxisSpec AxisSpec::Variable(std::vector<double> edges,
                            std::optional<AxisBoundaryType> boundaryType,
                            std::optional<AxisDirection> direction) {
  checkStrictlyIncreasing(edges, "AxisSpec::Variable");
  return AxisSpec(VariableParams{std::move(edges)}, boundaryType, direction);
}

AxisSpec AxisSpec::DeferredVariable(
    std::vector<double> normalizedEdges,
    std::optional<AxisBoundaryType> boundaryType,
    std::optional<AxisDirection> direction) {
  checkStrictlyIncreasing(normalizedEdges, "AxisSpec::DeferredVariable");
  if (normalizedEdges.front() != 0. || normalizedEdges.back() != 1.) {
    throw std::invalid_argument(
        "AxisSpec::DeferredVariable: edges must be normalized to [0, 1]");
  }
  return AxisSpec(DeferredVariableParams{std::move(normalizedEdges)},
                  boundaryType, direction);
}

AxisSpec AxisSpec::FromAxis(const IAxis& axis) {
  if (axis.getType() == AxisType::Equidistant) {
    return Equidistant(axis.getNBins(), axis.getMin(), axis.getMax(),
                       axis.getBoundaryType(), axis.getDirection());
  }
  return Variable(axis.getBinEdges(), axis.getBoundaryType(),
                  axis.getDirection());
}

AxisSpec AxisSpec::withDirection(AxisDirection direction) const {
  return AxisSpec(m_variant, m_boundaryType, direction);
}

AxisSpec AxisSpec::toDeferred() const {
  return std::visit(
      [this]<typename T>(const T& params) -> AxisSpec {
        if constexpr (std::is_same_v<T, EquidistantParams>) {
          return AxisSpec(EquidistantParams{params.nBins}, std::nullopt,
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
          return AxisSpec(DeferredVariableParams{std::move(normalizedEdges)},
                          std::nullopt, m_direction);
        } else {
          return AxisSpec(DeferredVariableParams{params.normalizedEdges},
                          std::nullopt, m_direction);
        }
      },
      m_variant);
}

bool AxisSpec::isDeferred() const {
  if (!m_boundaryType.has_value()) {
    return true;
  }
  if (std::holds_alternative<DeferredVariableParams>(m_variant)) {
    return true;
  }
  if (const auto* eq = std::get_if<EquidistantParams>(&m_variant);
      eq != nullptr) {
    return !eq->min.has_value() || !eq->max.has_value();
  }
  return false;
}

bool AxisSpec::isEquidistant() const {
  return std::holds_alternative<EquidistantParams>(m_variant);
}

bool AxisSpec::isVariable() const {
  return !isEquidistant();
}

bool AxisSpec::isDeferredVariable() const {
  return std::holds_alternative<DeferredVariableParams>(m_variant);
}

std::optional<AxisDirection> AxisSpec::direction() const {
  return m_direction;
}

std::optional<AxisBoundaryType> AxisSpec::boundaryType() const {
  return m_boundaryType;
}

std::size_t AxisSpec::nBins() const {
  return std::visit(
      []<typename T>(const T& params) -> std::size_t {
        if constexpr (std::is_same_v<T, EquidistantParams>) {
          return params.nBins;
        } else if constexpr (std::is_same_v<T, VariableParams>) {
          return params.edges.size() - 1;
        } else {
          return params.normalizedEdges.size() - 1;
        }
      },
      m_variant);
}

const AxisSpec::EquidistantParams& AxisSpec::asEquidistant() const {
  return std::get<EquidistantParams>(m_variant);
}

const AxisSpec::VariableParams& AxisSpec::asVariable() const {
  return std::get<VariableParams>(m_variant);
}

const AxisSpec::DeferredVariableParams& AxisSpec::asDeferredVariable() const {
  return std::get<DeferredVariableParams>(m_variant);
}

std::unique_ptr<IAxis> AxisSpec::buildAxis(const Options& options) const {
  AxisBoundaryType boundaryType = requireProperty(
      mergeProperty(m_boundaryType, options.boundaryType, "boundary type"),
      "boundary type");
  std::optional<AxisDirection> direction =
      mergeProperty(m_direction, options.direction, "direction");

  return std::visit(
      [&]<typename T>(const T& params) -> std::unique_ptr<IAxis> {
        // For absolute edges the range is implied by them, so merging only
        // validates a supplied one
        std::optional<double> min =
            mergeProperty(minOf(params), options.min, "minimum");
        std::optional<double> max =
            mergeProperty(maxOf(params), options.max, "maximum");

        if constexpr (std::is_same_v<T, VariableParams>) {
          return IAxis::createVariable(boundaryType, params.edges, direction);
        } else {
          double minValue = requireProperty(min, "minimum");
          double maxValue = requireProperty(max, "maximum");
          if (minValue >= maxValue) {
            throw std::invalid_argument(
                "AxisSpec::buildAxis: min must be < max");
          }

          if constexpr (std::is_same_v<T, EquidistantParams>) {
            return IAxis::createEquidistant(boundaryType, minValue, maxValue,
                                            params.nBins, direction);
          } else {
            std::vector<double> edges = params.normalizedEdges;
            for (double& edge : edges) {
              edge = minValue + edge * (maxValue - minValue);
            }
            return IAxis::createVariable(boundaryType, edges, direction);
          }
        }
      },
      m_variant);
}

std::string AxisSpec::toString() const {
  std::stringstream ss;
  ss << "AxisSpec: " << nBins() << " bins";
  if (m_direction.has_value()) {
    ss << " in " << axisDirectionName(*m_direction);
  }
  ss << (isEquidistant() ? ", equidistant" : ", variable");
  std::visit(
      [&ss](const auto& params) {
        std::optional<double> min = minOf(params);
        std::optional<double> max = maxOf(params);
        if (min.has_value() && max.has_value()) {
          ss << " within [" << *min << ", " << *max << "]";
        } else {
          ss << " within deferred range";
        }
      },
      m_variant);
  if (m_boundaryType.has_value()) {
    ss << ", " << *m_boundaryType;
  } else {
    ss << ", deferred boundary type";
  }
  return ss.str();
}

}  // namespace Acts
