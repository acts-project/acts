// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/GeometricExtent.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include <iomanip>
#include <ostream>

Acts::GeometricExtent::GeometricExtent(
    const std::array<std::array<ActsScalar, 2>, binValues>& envelope)
    : m_constrains(0), m_envelope(envelope) {
  m_range[binR] =
      Range1D<ActsScalar>(0., std::numeric_limits<ActsScalar>::max());
  m_range[binPhi] = Range1D<ActsScalar>(-M_PI, M_PI);
  m_range[binRPhi] =
      Range1D<ActsScalar>(0., std::numeric_limits<ActsScalar>::max());
  m_range[binMag] =
      Range1D<ActsScalar>(0., std::numeric_limits<ActsScalar>::max());
}

void Acts::GeometricExtent::extend(const Vector3& vtx,
                                   const std::vector<BinningValue>& bValues,
                                   bool applyEnv) {
  auto adapt = [&](BinningValue bValue) {
    ActsScalar cValue = VectorHelpers::cast(vtx, bValue);
    ActsScalar lEnv = applyEnv ? m_envelope[bValue][0] : 0.;
    ActsScalar hEnv = applyEnv ? m_envelope[bValue][1] : 0.;
    ActsScalar mValue = cValue - lEnv;
    // Special protection for radial value 
    if (bValue == binR and mValue < 0.) {
      mValue = std::max(mValue, 0.);
    }
    if (constrains(bValue)) {
      m_range[bValue].expand(mValue, cValue + hEnv);
    } else {
      m_range[bValue].shrink(mValue, cValue + hEnv);
    }
    m_constrains.set(bValue);
  };

  for (auto bValue : bValues) {
    adapt(bValue);
  }
}

void Acts::GeometricExtent::extend(const GeometricExtent& rhs,
                                   const std::vector<BinningValue>& bValues,
                                   bool applyEnv) {
  for (auto bValue : bValues) {
    // The value is constraint, envelope can be optional 
    if (rhs.constrains(bValue)) {
      ActsScalar lEnv = applyEnv ? m_envelope[bValue][0] : 0.;
      ActsScalar hEnv = applyEnv ? m_envelope[bValue][1] : 0.;
      if (constrains(bValue)) {
        m_range[bValue].expand(rhs.range()[bValue].min() - lEnv,
                               rhs.range()[bValue].max() + hEnv);
      } else {
        m_range[bValue].shrink(rhs.range()[bValue].min() - lEnv,
                               rhs.range()[bValue].max() + hEnv);
      }
      m_constrains.set(bValue);
    } else if (rhs.envelope()[bValue] != zeroEnvelope) {
      // Only an envelope given, but value is not contraint -> apply envelope
      m_range[bValue].expand(m_range[bValue].min() - rhs.envelope()[bValue][0],
                             m_range[bValue].max() + rhs.envelope()[bValue][1]);
      m_constrains.set(bValue);
    }
  }
}

void Acts::GeometricExtent::set(BinningValue bValue, ActsScalar min,
                                ActsScalar max) {
  ActsScalar minval = min;
  if (bValue == binR and minval < 0.) {
    minval = 0.;
  }
  m_range[bValue] = Range1D{minval, max};
  m_constrains.set(bValue);
}

void Acts::GeometricExtent::setEnvelope(const ExtentEnvelope& envelope) {
  m_envelope = envelope;
}

const Acts::Range1D<Acts::ActsScalar>& Acts::GeometricExtent::range(
    BinningValue bValue) const {
  return m_range[bValue];
}

bool Acts::GeometricExtent::contains(const GeometricExtent& rhs,
                                     BinningValue bValue) const {
  // Helper to check including a constraint bit set check
  auto checkContainment = [&](BinningValue bvc) -> bool {
    if (not constrains(bvc)) {
      return true;
    }
    return (rhs.range()[bvc] <= m_range[bvc]);
  };

  // Check all
  if (bValue == binValues) {
    for (int ibv = 0; ibv < (int)binValues; ++ibv) {
      if (not checkContainment((BinningValue)ibv)) {
        return false;
      }
    }
    return true;
  }
  // Check specific
  return checkContainment(bValue);
}

bool Acts::GeometricExtent::intersects(const GeometricExtent& rhs,
                                       BinningValue bValue) const {
  // Helper to check including a constraint bit set check
  auto checkIntersect = [&](BinningValue bvc) -> bool {
    if (not constrains(bvc) or not rhs.constrains(bvc)) {
      return false;
    }
    return (m_range[bvc] && rhs.range()[bvc]);
  };

  // Check all
  if (bValue == binValues) {
    for (int ibv = 0; ibv < (int)binValues; ++ibv) {
      if (checkIntersect((BinningValue)ibv)) {
        return true;
      }
    }
    return false;
  }
  // Check specific
  return checkIntersect(bValue);
}

bool Acts::GeometricExtent::constrains(BinningValue bValue) const {
  if (bValue == binValues) {
    return (m_constrains.count() > 0);
  }
  return m_constrains.test(size_t(bValue));
}

std::ostream& Acts::GeometricExtent::toStream(std::ostream& sl) const {
  sl << "GeometricExtent in space : " << std::endl;
  for (size_t ib = 0; ib < static_cast<size_t>(binValues); ++ib) {
    if (constrains((BinningValue)ib)) {
      sl << "  - value :" << std::setw(10) << binningValueNames()[ib]
         << " | range = [" << m_range[ib].min() << ", " << m_range[ib].max()
         << "]" << std::endl;
    }
  }
  return sl;
}

// Overload of << operator for std::ostream for debug output
std::ostream& Acts::operator<<(std::ostream& sl, const GeometricExtent& rhs) {
  return rhs.toStream(sl);
}
