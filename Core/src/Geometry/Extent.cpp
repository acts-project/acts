// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Extent.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <numbers>

Acts::Extent::Extent(const ExtentEnvelope& envelope)
    : m_constrains(0), m_envelope(envelope) {
  m_range[toUnderlying(AxisDirection::AxisR)] =
      Range1D<double>(0., std::numeric_limits<double>::max());
  m_range[toUnderlying(AxisDirection::AxisPhi)] = Range1D<double>(
      -std::numbers::pi_v<double>, std::numbers::pi_v<double>);
  m_range[toUnderlying(AxisDirection::AxisRPhi)] =
      Range1D<double>(0., std::numeric_limits<double>::max());
  m_range[toUnderlying(AxisDirection::AxisMag)] =
      Range1D<double>(0., std::numeric_limits<double>::max());
}

void Acts::Extent::extend(const Vector3& vtx,
                          const std::vector<AxisDirection>& bValues,
                          bool applyEnv, bool fillHistograms) {
  for (auto bValue : bValues) {
    // Get the casted value given the binning value description
    double cValue = VectorHelpers::cast(vtx, bValue);
    if (fillHistograms) {
      m_valueHistograms[toUnderlying(bValue)].push_back(cValue);
    }
    // Apply envelope as suggested
    double lEnv = applyEnv ? m_envelope[bValue][0] : 0.;
    double hEnv = applyEnv ? m_envelope[bValue][1] : 0.;
    double mValue = cValue - lEnv;
    // Special protection for radial value
    if (bValue == AxisDirection::AxisR && mValue < 0.) {
      mValue = std::max(mValue, 0.);
    }
    if (constrains(bValue)) {
      m_range[toUnderlying(bValue)].expand(mValue, cValue + hEnv);
    } else {
      m_range[toUnderlying(bValue)].shrink(mValue, cValue + hEnv);
    }
    m_constrains.set(toUnderlying(bValue));
  }
}

void Acts::Extent::extend(const Extent& rhs,
                          const std::vector<AxisDirection>& bValues,
                          bool applyEnv) {
  for (auto bValue : bValues) {
    // The value is constraint, envelope can be optional
    if (rhs.constrains(bValue)) {
      double lEnv = applyEnv ? m_envelope[bValue][0] : 0.;
      double hEnv = applyEnv ? m_envelope[bValue][1] : 0.;
      if (constrains(bValue)) {
        m_range[toUnderlying(bValue)].expand(
            rhs.range()[toUnderlying(bValue)].min() - lEnv,
            rhs.range()[toUnderlying(bValue)].max() + hEnv);
      } else {
        m_range[toUnderlying(bValue)].shrink(
            rhs.range()[toUnderlying(bValue)].min() - lEnv,
            rhs.range()[toUnderlying(bValue)].max() + hEnv);
      }
      m_constrains.set(toUnderlying(bValue));
    } else if (rhs.envelope()[bValue] != zeroEnvelope) {
      // Only an envelope given, but value is not constraint -> apply envelope
      m_range[toUnderlying(bValue)].expand(
          m_range[toUnderlying(bValue)].min() - rhs.envelope()[bValue][0],
          m_range[toUnderlying(bValue)].max() + rhs.envelope()[bValue][1]);
      m_constrains.set(toUnderlying(bValue));
    }
  }
}

void Acts::Extent::addConstrain(const Acts::Extent& rhs,
                                const ExtentEnvelope& envelope) {
  for (const auto& bValue : allAxisDirections()) {
    if (rhs.constrains(bValue) && !constrains(bValue)) {
      const auto& cRange = rhs.range(bValue);
      m_range[toUnderlying(bValue)].setMin(cRange.min() - envelope[bValue][0u]);
      m_range[toUnderlying(bValue)].setMax(cRange.max() + envelope[bValue][1u]);
      m_constrains.set(toUnderlying(bValue));
    }
  }
}

void Acts::Extent::set(AxisDirection bValue, double min, double max) {
  double minval = min;
  if (bValue == AxisDirection::AxisR && minval < 0.) {
    minval = 0.;
  }
  m_range[toUnderlying(bValue)] = Range1D<double>{minval, max};
  m_constrains.set(toUnderlying(bValue));
}

void Acts::Extent::setMin(AxisDirection bValue, double min) {
  double minval = min;
  if (bValue == AxisDirection::AxisR && minval < 0.) {
    minval = 0.;
  }
  m_range[toUnderlying(bValue)].setMin(0u, minval);
  m_constrains.set(toUnderlying(bValue));
}

void Acts::Extent::setMax(AxisDirection bValue, double max) {
  m_range[toUnderlying(bValue)].setMax(0u, max);
  m_constrains.set(toUnderlying(bValue));
}

void Acts::Extent::setEnvelope(const ExtentEnvelope& envelope) {
  m_envelope = envelope;
}

bool Acts::Extent::contains(const Vector3& vtx) const {
  Extent checkExtent;
  for (const auto& bv : allAxisDirections()) {
    if (constrains(bv)) {
      double vtxVal = VectorHelpers::cast(vtx, bv);
      checkExtent.set(bv, vtxVal, vtxVal);
    }
  }
  return contains(checkExtent);
}

bool Acts::Extent::contains(const Extent& rhs,
                            std::optional<AxisDirection> bValue) const {
  // Helper to check including a constraint bit set check
  auto checkContainment = [&](AxisDirection bvc) -> bool {
    if (!constrains(bvc)) {
      return true;
    }
    return (rhs.range()[toUnderlying(bvc)] <= m_range[toUnderlying(bvc)]);
  };

  // Check all
  if (!bValue.has_value()) {
    return std::ranges::all_of(allAxisDirections(), checkContainment);
  }
  // Check specific
  return checkContainment(bValue.value());
}

bool Acts::Extent::intersects(const Extent& rhs,
                              std::optional<AxisDirection> bValue) const {
  // Helper to check including a constraint bit set check
  auto checkIntersect = [&](AxisDirection bvc) -> bool {
    if (!constrains(bvc) || !rhs.constrains(bvc)) {
      return false;
    }
    return (m_range[toUnderlying(bvc)] && rhs.range()[toUnderlying(bvc)]);
  };

  // Check all
  if (!bValue.has_value()) {
    return std::ranges::any_of(allAxisDirections(), checkIntersect);
  }
  // Check specific
  return checkIntersect(bValue.value());
}

bool Acts::Extent::constrains(AxisDirection bValue) const {
  return m_constrains.test(static_cast<std::size_t>(bValue));
}

bool Acts::Extent::constrains() const {
  return m_constrains.count() > 0;
}

bool Acts::Extent::operator==(const Extent& e) const {
  if (m_constrains != e.m_constrains) {
    return false;
  }
  if (m_envelope != e.m_envelope) {
    return false;
  }
  if (!(m_range == e.m_range)) {
    return false;
  }
  if (m_valueHistograms != e.m_valueHistograms) {
    return false;
  }
  return true;
}

std::string Acts::Extent::toString(const std::string& indent) const {
  std::stringstream sl;
  sl << indent << "Extent in space :" << std::endl;
  for (const auto& bv : allAxisDirections()) {
    if (constrains(bv)) {
      sl << indent << "  - value :" << std::setw(10) << bv << " | range = ["
         << m_range[toUnderlying(bv)].min() << ", "
         << m_range[toUnderlying(bv)].max() << "]" << std::endl;
    }
  }
  return sl.str();
}

// Overload of << operator for std::ostream for debug output
std::ostream& Acts::operator<<(std::ostream& sl, const Extent& rhs) {
  sl << rhs.toString();
  return sl;
}
