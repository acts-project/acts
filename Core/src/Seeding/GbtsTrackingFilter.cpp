// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsTrackingFilter.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"

#include <cmath>
#include <memory>
#include <utility>

namespace Acts::Experimental {

namespace {

ParabolaSpacePointFilter makeFilter(const GbtsTrackingFilter::Config& cfg) {
  ParabolaSpacePointFilter::Config filterCfg;
  filterCfg.initialVarianceTransverse = cfg.initialVarianceX;
  filterCfg.initialVarianceLongitudinal = cfg.initialVarianceY;

  ParabolaPropagation propagation;
  propagation.cfg.sigmaMS = cfg.sigmaMS;
  propagation.cfg.radLen = cfg.radLen;
  propagation.cfg.referenceCurvature = cfg.maxCurvature;

  ConfiguredResolution resolution;
  resolution.cfg.sigmaX = cfg.sigmaX;
  resolution.cfg.sigmaY = cfg.sigmaY;
  resolution.cfg.sigmaXStrip = cfg.sigmaXStrip;

  return ParabolaSpacePointFilter(filterCfg, propagation, resolution);
}

}  // namespace

GbtsTrackingFilter::GbtsTrackingFilter(
    const Config& config, const std::shared_ptr<const GbtsGeometry>& geometry,
    std::unique_ptr<const Logger> logger)
    : m_cfg(config),
      m_filter(makeFilter(config)),
      m_geometry(geometry),
      m_logger(std::move(logger)) {}

SpacePointFilterInput GbtsTrackingFilter::input(
    const detail::GbtsNodeView& nodeView, const SpacePointIndex index) const {
  const detail::GbtsNodeProxy node = nodeView[index];
  SpacePointFilterInput out;
  out.x = node.x();
  out.y = node.y();
  out.z = node.z();
  out.r = node.r();
  out.measuresR = getLayerType(node.layer()) == GbtsLayerType::Endcap;
  if (const auto* strip = nodeView.strip(index); strip != nullptr) {
    out.stripHalfVector = strip->outerHalfVector;
  }
  return out;
}

detail::GbtsEdgeState GbtsTrackingFilter::initialize(
    const detail::GbtsNodeView& nodeView, const detail::GbtsEdge& pS) const {
  // n2->n1
  return {.j = 0,
          .filter = m_filter.initialize(input(nodeView, pS.n2),
                                        input(nodeView, pS.n1))};
}

bool GbtsTrackingFilter::update(const detail::GbtsNodeView& nodeView,
                                const detail::GbtsEdge& pS,
                                detail::GbtsEdgeState& ts) const {
  const DecoupledSpacePointFilterState& fs = ts.filter;
  if (fs.cx[5] < 0 || fs.cx[3] < 0 || fs.cx[0] < 0) {
    ACTS_DEBUG("Negative cov_x");
  }

  if (fs.cy[2] < 0 || fs.cy[0] < 0) {
    ACTS_DEBUG("Negative cov_y");
  }

  const DecoupledSpacePointFilterPrediction prediction = m_filter.predict(
      ts.filter, input(nodeView, pS.n2), input(nodeView, pS.n1));

  const float dchi2_x = prediction.chi2Transverse();
  const float dchi2_y = prediction.chi2Longitudinal();

  if (dchi2_x > m_cfg.maxDChi2X || dchi2_y > m_cfg.maxDChi2Y) {
    return false;
  }

  // The path search in GraphBasedTrackSeeder::extractBestPathPerEdge skips
  // paths that cannot beat the best one found, assuming a hit adds at most
  // addHit here. A score that can add more must change or remove that bound.
  ts.j += m_cfg.addHit - dchi2_x * m_cfg.weightX - dchi2_y * m_cfg.weightY;

  m_filter.apply(ts.filter, prediction);

  if (std::abs(ts.filter.x[2]) > m_cfg.maxCurvature) {
    return false;
  }

  const float z0 = ts.filter.y[0] - ts.filter.refL * ts.filter.y[1];

  const GbtsZ0Range& z0Range = m_geometry->z0Range();
  if (z0 < z0Range.min || z0 > z0Range.max) {
    return false;
  }

  return true;
}

GbtsLayerType GbtsTrackingFilter::getLayerType(
    const GbtsLayerIndex layerIndex) const {
  return m_geometry->layerDescription(layerIndex).type;
}

}  // namespace Acts::Experimental
