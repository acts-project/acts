// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/VectorMultiTrajectory.hpp"

#include <boost/histogram.hpp>
#include <boost/histogram/make_histogram.hpp>

namespace Acts {

VectorMultiTrajectory::DynamicColumnBase::~DynamicColumnBase() = default;

auto VectorMultiTrajectory::addTrackState_impl(TrackStatePropMask mask,
                                               IndexType iprevious)
    -> IndexType {
  using PropMask = TrackStatePropMask;

  m_index.emplace_back();
  IndexData& p = m_index.back();
  size_t index = m_index.size() - 1;

  if (iprevious != kInvalid) {
    p.iprevious = iprevious;
  }

  // always set, but can be null
  m_referenceSurfaces.emplace_back(nullptr);

  assert(m_params.size() == m_cov.size());

  if (ACTS_CHECK_BIT(mask, PropMask::Predicted)) {
    m_params.emplace_back();
    m_cov.emplace_back();
    p.ipredicted = m_params.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Filtered)) {
    m_params.emplace_back();
    m_cov.emplace_back();
    p.ifiltered = m_params.size() - 1;
  }

  if (ACTS_CHECK_BIT(mask, PropMask::Smoothed)) {
    m_params.emplace_back();
    m_cov.emplace_back();
    p.ismoothed = m_params.size() - 1;
  }

  assert(m_params.size() == m_cov.size());

  if (ACTS_CHECK_BIT(mask, PropMask::Jacobian)) {
    m_jac.emplace_back();
    p.ijacobian = m_jac.size() - 1;
  }

  m_sourceLinks.push_back(nullptr);
  p.iuncalibrated = m_sourceLinks.size() - 1;

  if (ACTS_CHECK_BIT(mask, PropMask::Calibrated)) {
    m_meas.emplace_back();
    m_measCov.emplace_back();
    p.icalibrated = m_meas.size() - 1;

    m_sourceLinks.push_back(nullptr);
    p.icalibratedsourcelink = m_sourceLinks.size() - 1;

    m_projectors.emplace_back();
    p.iprojector = m_projectors.size() - 1;
  }

  // dynamic columns
  for (auto& [key, vec] : m_dynamic) {
    vec->add();
  }

  return index;
}

void VectorMultiTrajectory::shareFrom_impl(IndexType iself, IndexType iother,
                                           TrackStatePropMask shareSource,
                                           TrackStatePropMask shareTarget) {
  // auto other = getTrackState(iother);
  IndexData& self = m_index[iself];
  IndexData& other = m_index[iother];

  assert(ACTS_CHECK_BIT(getTrackState(iother).getMask(), shareSource) &&
         "Source has incompatible allocation");

  using PM = TrackStatePropMask;

  IndexType sourceIndex{kInvalid};
  switch (shareSource) {
    case PM::Predicted:
      sourceIndex = other.ipredicted;
      break;
    case PM::Filtered:
      sourceIndex = other.ifiltered;
      break;
    case PM::Smoothed:
      sourceIndex = other.ismoothed;
      break;
    case PM::Jacobian:
      sourceIndex = other.ijacobian;
      break;
    default:
      throw std::domain_error{"Unable to share this component"};
  }

  assert(sourceIndex != kInvalid);

  switch (shareTarget) {
    case PM::Predicted:
      assert(shareSource != PM::Jacobian);
      self.ipredicted = sourceIndex;
      break;
    case PM::Filtered:
      assert(shareSource != PM::Jacobian);
      self.ifiltered = sourceIndex;
      break;
    case PM::Smoothed:
      assert(shareSource != PM::Jacobian);
      self.ismoothed = sourceIndex;
      break;
    case PM::Jacobian:
      assert(shareSource == PM::Jacobian);
      self.ijacobian = sourceIndex;
      break;
    default:
      throw std::domain_error{"Unable to share this component"};
  }
}

void VectorMultiTrajectory::unset_impl(TrackStatePropMask target,
                                       IndexType istate) {
  using PM = TrackStatePropMask;

  switch (target) {
    case PM::Predicted:
      m_index[istate].ipredicted = kInvalid;
      break;
    case PM::Filtered:
      m_index[istate].ifiltered = kInvalid;
      break;
    case PM::Smoothed:
      m_index[istate].ismoothed = kInvalid;
      break;
    case PM::Jacobian:
      m_index[istate].ijacobian = kInvalid;
      break;
    case PM::Calibrated:
      m_index[istate].icalibrated = kInvalid;
      break;
    default:
      throw std::domain_error{"Unable to unset this component"};
  }
}

void VectorMultiTrajectory::clear_impl() {
  m_index.clear();
  m_params.clear();
  m_cov.clear();
  m_meas.clear();
  m_measCov.clear();
  m_jac.clear();
  m_sourceLinks.clear();
  m_projectors.clear();
  m_referenceSurfaces.clear();
  for (auto& [key, vec] : m_dynamic) {
    vec->clear();
  }
}

// VectorMultiTrajectory::~VectorMultiTrajectory() {
// if (m_index.empty()) {
// return;
// }

// size_t total = 0;
// std::cout << "VMT Size:" << std::endl;
// #define PRINT_SIZE(x)                                                   \
  // do {                                                                  \
    // constexpr size_t esize = sizeof(decltype(x)::value_type);           \
    // total += esize * x.size();                                          \
    // std::cout << #x << " size: " << esize << "b * " << x.size() << "->" \
              // << (esize * x.size()) / 1024 / 1024 << "M" << std::endl;  \
  // } while (0)

// PRINT_SIZE(m_index);
// PRINT_SIZE(m_previous);
// PRINT_SIZE(m_params);
// PRINT_SIZE(m_cov);
// PRINT_SIZE(m_meas);
// PRINT_SIZE(m_measCov);
// PRINT_SIZE(m_jac);
// PRINT_SIZE(m_sourceLinks);
// PRINT_SIZE(m_projectors);

// #undef PRINT_SIZE

// std::cout << "total: " << total / 1024 / 1024 << "M" << std::endl;
// std::cout << "---" << std::endl;
// }

auto VectorMultiTrajectory::statistics() const -> Statistics {
  using namespace boost::histogram;
  using cat = axis::category<std::string>;

  Statistics::axes_t axes;
  axes.emplace_back(cat({"count", "index", "params", "cov", "meas", "measCov",
                         "jac", "sourceLinks", "projectors"}));

  Statistics::hist_t h = make_histogram(std::move(axes));
  std::cout << "rank: " << h.rank() << std::endl;

#define FILL_SIZE(x)                                              \
  do {                                                            \
    constexpr size_t esize = sizeof(decltype(m_##x)::value_type); \
    h(#x, weight(esize* m_##x.size()));                           \
  } while (0)

  FILL_SIZE(index);
  FILL_SIZE(params);
  FILL_SIZE(cov);
  FILL_SIZE(meas);
  FILL_SIZE(measCov);
  FILL_SIZE(jac);
  FILL_SIZE(sourceLinks);
  FILL_SIZE(projectors);

#undef FILL_SIZE

  h("count", weight(m_index.size()));

  return Statistics{h};
}

std::ostream& operator<<(std::ostream& os,
                         const VectorMultiTrajectory::Statistics& stats) {
  using namespace boost::histogram;
  using cat = axis::category<std::string>;

  auto& h = stats.hist;

  auto column_axis = axis::get<cat>(h.axis(0));

  double total = 0;

  auto p = [&](const auto& key, const double v) {
    os << std::setw(20) << key << ": ";
    os << std::fixed << std::setw(8) << std::setprecision(2)
       << (v / 1024 / 1024) << "M" << std::endl;
  };

  for (auto&& x : indexed(h)) {
    p(column_axis.bin(x.index(0)), *x);
    total += *x;
  }

  p("total", total);

  return os;
}

}  // namespace Acts
