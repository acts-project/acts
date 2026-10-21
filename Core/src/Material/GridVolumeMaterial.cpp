// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/GridVolumeMaterial.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>

namespace Acts {

GridVolumeMaterial::GridVolumeMaterial(MultiAxisSpec3D binning, Storage storage,
                                       const Transform3& transform)
    : m_binning(std::move(binning)),
      m_multiAxis(m_binning.buildMultiAxis()),
      m_storage(std::move(storage)),
      m_transform(transform) {
  const auto bins = m_multiAxis->getNTotalBins(true);
  std::visit(
      [bins](const auto& values) {
        using T = std::decay_t<decltype(values)>;
        if constexpr (std::is_same_v<T, Direct>) {
          if (values.size() != bins) {
            throw std::invalid_argument(
                "GridVolumeMaterial: invalid storage size");
          }
        } else {
          if (values.indices.size() != bins) {
            throw std::invalid_argument(
                "GridVolumeMaterial: invalid index count");
          }
          std::size_t size = 0;
          if constexpr (std::is_same_v<T, GloballyIndexed>) {
            if (!values.material) {
              throw std::invalid_argument(
                  "GridVolumeMaterial: null material palette");
            }
            size = values.material->size();
          } else {
            size = values.material.size();
          }
          for (auto index : values.indices) {
            if (index >= size) {
              throw std::invalid_argument(
                  "GridVolumeMaterial: invalid material index");
            }
          }
        }
      },
      m_storage);

  std::array<AxisDirection, 3> directions{
      AxisDirection::AxisX, AxisDirection::AxisY, AxisDirection::AxisZ};
  const bool directed = m_binning.axisSpec(0).direction().has_value();
  Vector3 periods = Vector3::Zero();
  Vector3 origins = Vector3::Zero();
  for (std::size_t i = 0; i < 3; ++i) {
    const auto dir = m_binning.axisSpec(i).direction();
    if (dir.has_value() != directed) {
      throw std::invalid_argument(
          "GridVolumeMaterial: specify all axis directions or none");
    }
    if (dir) {
      directions[i] = *dir;
      if (static_cast<int>(*dir) < 0 ||
          static_cast<std::size_t>(*dir) >= numAxisDirections()) {
        throw std::invalid_argument(
            "GridVolumeMaterial: invalid axis direction");
      }
    }
    const auto& axis = m_multiAxis->getAxis(i);
    if (axis.getBoundaryType() == AxisBoundaryType::Closed) {
      periods[i] = axis.getMax() - axis.getMin();
      origins[i] = axis.getMin();
    }
  }
  if (directions[0] == directions[1] || directions[0] == directions[2] ||
      directions[1] == directions[2]) {
    throw std::invalid_argument(
        "GridVolumeMaterial: duplicate axis directions");
  }
  m_globalToGrid = std::make_shared<const VolumeMaterialCell::Transform>(
      [inverse = Transform3(transform.inverse()), directions, periods,
       origins](const Vector3& position) -> Vector3 {
        const Vector3 local = inverse * position;
        Vector3 result;
        for (std::size_t i = 0; i < 3; ++i) {
          double value = VectorHelpers::cast(local, directions[i]);
          if (periods[i] > 0.) {
            value -= periods[i] * std::floor((value - origins[i]) / periods[i]);
          }
          result[i] = value;
        }
        return result;
      });
}

const Material& GridVolumeMaterial::materialAt(std::size_t bin) const {
  return std::visit(
      overloaded{[bin](const Direct& values) -> const Material& {
                   return values.at(bin);
                 },
                 [bin](const Indexed& values) -> const Material& {
                   return values.material.at(values.indices.at(bin));
                 },
                 [bin](const GloballyIndexed& values) -> const Material& {
                   return values.material->at(values.indices.at(bin));
                 }},
      m_storage);
}

const Material GridVolumeMaterial::material(const Vector3& position) const {
  const Vector3 point = (*m_globalToGrid)(position);
  if (!point.allFinite()) {
    throw std::invalid_argument(
        "GridVolumeMaterial: non-finite lookup position");
  }
  return materialAt(
      m_multiAxis->getGlobalBinFromPoint({point[0], point[1], point[2]}));
}

VolumeMaterialCell GridVolumeMaterial::materialCell(
    const Vector3& position) const {
  const Vector3 point = (*m_globalToGrid)(position);
  if (!point.allFinite()) {
    throw std::invalid_argument(
        "GridVolumeMaterial: non-finite lookup position");
  }
  const auto localBins =
      m_multiAxis->getLocalBinsFromPoint({point[0], point[1], point[2]});
  Vector3 lower, upper;
  constexpr double infinity = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < 3; ++i) {
    const auto& axis = m_multiAxis->getAxis(i);
    const auto bin = localBins[i];
    const auto last = axis.getNBins();
    const bool bound = axis.getBoundaryType() == AxisBoundaryType::Bound;
    lower[i] = (bin == 0 || (bound && bin == 1))
                   ? -infinity
                   : (bin > last ? axis.getMax() : axis.getBinLowerBound(bin));
    upper[i] = (bin > last || (bound && bin == last))
                   ? infinity
                   : (bin == 0 ? axis.getMin() : axis.getBinUpperBound(bin));
  }
  return VolumeMaterialCell(
      materialAt(m_multiAxis->getGlobalBinFromLocalBins(localBins)), lower,
      upper, m_globalToGrid);
}

std::ostream& GridVolumeMaterial::toStream(std::ostream& stream) const {
  return stream << "Acts::GridVolumeMaterial : " << m_binning;
}

}  // namespace Acts
