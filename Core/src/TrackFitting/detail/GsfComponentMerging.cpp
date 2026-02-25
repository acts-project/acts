// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfComponentMerging.hpp"

#include <iostream>

namespace Acts {

std::tuple<BoundVector, BoundMatrix> detail::Gsf::mergeGaussianMixture(
    std::span<const GsfComponent> mixture, const Surface &surface,
    ComponentMergeMethod method) {
  return mergeGaussianMixture(
      mixture,
      [](const GsfComponent &c) {
        return std::tie(c.weight, c.boundPars, c.boundCov);
      },
      surface, method);
}

GsfComponent detail::Gsf::mergeTwoComponents(const GsfComponent &a,
                                             const GsfComponent &b,
                                             const Surface &surface) {
  assert(a.weight >= 0.0 && b.weight >= 0.0 && "non-positive weight");

  std::array components = {&a, &b};
  const auto proj = [](const GsfComponent *c) {
    return std::tie(c->weight, c->boundPars, c->boundCov);
  };
  auto [mergedPars, mergedCov] =
      angleDescriptionSwitch(surface, [&](const auto &desc) {
        return mergeGaussianMixtureMeanCov(components, proj, desc);
      });

  GsfComponent ret = a;
  ret.boundPars = mergedPars;
  ret.boundCov = mergedCov;
  ret.weight = a.weight + b.weight;
  return ret;
}

namespace detail::Gsf {

double SymmetricKLDistanceMatrix::computeSymmetricKlDivergence(
    const GsfComponent &a, const GsfComponent &b) {
  const double parsA = a.boundPars[eBoundQOverP];
  const double parsB = b.boundPars[eBoundQOverP];
  const double covA = a.boundCov(eBoundQOverP, eBoundQOverP);
  const double covB = b.boundCov(eBoundQOverP, eBoundQOverP);

  assert(covA != 0.0);
  assert(std::isfinite(covA));
  assert(covB != 0.0);
  assert(std::isfinite(covB));

  const double kl = covA * (1 / covB) + covB * (1 / covA) +
                    (parsA - parsB) * (1 / covA + 1 / covB) * (parsA - parsB);

  assert(kl >= 0.0 && "kl-divergence must be non-negative");

  return kl;
}

SymmetricKLDistanceMatrix::SymmetricKLDistanceMatrix(
    std::span<const GsfComponent> cmps)
    : m_distances(Array::Zero(cmps.size() * (cmps.size() - 1) / 2)),
      m_mask(Mask::Ones(cmps.size() * (cmps.size() - 1) / 2)),
      m_mapToPair(m_distances.size()),
      m_numberComponents(cmps.size()) {
  for (std::size_t i = 1; i < m_numberComponents; ++i) {
    const std::size_t indexConst = (i - 1) * i / 2;
    for (std::size_t j = 0; j < i; ++j) {
      m_mapToPair.at(indexConst + j) = {i, j};
      m_distances[indexConst + j] =
          computeSymmetricKlDivergence(cmps[i], cmps[j]);
    }
  }
}

double SymmetricKLDistanceMatrix::at(std::size_t i, std::size_t j) const {
  return m_distances[i * (i - 1) / 2 + j];
}

void SymmetricKLDistanceMatrix::recomputeAssociatedDistances(
    std::size_t n, std::span<const GsfComponent> cmps) {
  assert(cmps.size() == m_numberComponents && "size mismatch");

  setAssociated(n, m_distances, [&](std::size_t i, std::size_t j) {
    return computeSymmetricKlDivergence(cmps[i], cmps[j]);
  });
}

void SymmetricKLDistanceMatrix::maskAssociatedDistances(std::size_t n) {
  setAssociated(n, m_mask, [&](std::size_t, std::size_t) { return false; });
}

std::pair<std::size_t, std::size_t> SymmetricKLDistanceMatrix::minDistancePair()
    const {
  double min = std::numeric_limits<double>::max();
  std::size_t idx = 0;

  for (std::size_t i = 0; i < static_cast<std::size_t>(m_distances.size());
       ++i) {
    if (double new_min = std::min(min, m_distances[i]);
        m_mask[i] && new_min < min) {
      min = new_min;
      idx = i;
    }
  }

  return m_mapToPair.at(idx);
}

std::ostream &SymmetricKLDistanceMatrix::toStream(std::ostream &os) const {
  const auto prev_precision = os.precision();
  const int width = 8;
  const int prec = 2;

  os << "\n";
  os << std::string(width, ' ') << " | ";
  for (std::size_t j = 0ul; j < m_numberComponents - 1; ++j) {
    os << std::setw(width) << j << "  ";
  }
  os << "\n";
  os << std::string((width + 3) + (width + 2) * (m_numberComponents - 1), '-');
  os << "\n";

  for (std::size_t i = 1ul; i < m_numberComponents; ++i) {
    const std::size_t indexConst = (i - 1) * i / 2;
    os << std::setw(width) << i << " | ";
    for (std::size_t j = 0ul; j < i; ++j) {
      os << std::setw(width) << std::setprecision(prec)
         << m_distances[indexConst + j] << "  ";
    }
    os << "\n";
  }
  os << std::setprecision(prev_precision);
  return os;
}

}  // namespace detail::Gsf

}  // namespace Acts
