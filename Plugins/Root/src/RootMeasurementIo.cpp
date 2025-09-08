// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMeasurementIo.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TTree.h>

Acts::RootMeasurementIo::RootMeasurementIo(const Config& config)
    : m_cfg(config) {
  clear();
}

void Acts::RootMeasurementIo::connectForWrite(TTree& measurementTree) {
  measurementTree.Branch("event_nr", &m_payload.eventNr);
  measurementTree.Branch("volume_id", &m_payload.volumeID);
  measurementTree.Branch("layer_id", &m_payload.layerID);
  measurementTree.Branch("surface_id", &m_payload.surfaceID);
  measurementTree.Branch("extra_id", &m_payload.extraID);

  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("rec_" + bNames[ib]).c_str(),
                           &m_payload.recBound[ib]);
  }
  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("var_" + bNames[ib]).c_str(),
                           &m_payload.varBound[ib]);
  }

  measurementTree.Branch("rec_gx", &m_payload.recGx);
  measurementTree.Branch("rec_gy", &m_payload.recGy);
  measurementTree.Branch("rec_gz", &m_payload.recGz);

  measurementTree.Branch("clus_size", &m_payload.nch);
  measurementTree.Branch("channel_value", &m_payload.chValue);
  // Both are allocated, but only relevant ones are set
  for (auto ib : m_cfg.clusterIndices) {
    if (static_cast<unsigned int>(ib) < 2) {
      measurementTree.Branch(("channel_" + bNames[ib]).c_str(),
                             &m_payload.chId[ib]);
      measurementTree.Branch(("clus_size_" + bNames[ib]).c_str(),
                             &m_payload.cSize[ib]);
    }
  }

  for (unsigned int ib = 0; ib < eBoundSize; ++ib) {
    measurementTree.Branch(("true_" + bNames[ib]).c_str(),
                           &m_payload.trueBound[ib]);
  }
  measurementTree.Branch("true_x", &m_payload.trueGx);
  measurementTree.Branch("true_y", &m_payload.trueGy);
  measurementTree.Branch("true_z", &m_payload.trueGz);
  measurementTree.Branch("true_incident_phi", &m_payload.incidentPhi);
  measurementTree.Branch("true_incident_theta", &m_payload.incidentTheta);

  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("residual_" + bNames[ib]).c_str(),
                           &m_payload.residual[ib]);
  }
  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("pull_" + bNames[ib]).c_str(), &m_payload.pull[ib]);
  }
  clear();
}

void Acts::RootMeasurementIo::fillIdentification(
    int evnt, const GeometryIdentifier& geoId) {
  m_payload.eventNr = evnt;
  m_payload.volumeID = geoId.volume();
  m_payload.layerID = geoId.layer();
  m_payload.surfaceID = geoId.sensitive();
  m_payload.extraID = geoId.extra();
}

void Acts::RootMeasurementIo::fillTruthParameters(
    const Vector2& lp, const Vector4& xt, const Vector3& dir,
    const std::pair<double, double> angles) {
  m_payload.trueBound[eBoundLoc0] = lp[eBoundLoc0];
  m_payload.trueBound[eBoundLoc1] = lp[eBoundLoc1];
  m_payload.trueBound[eBoundPhi] = VectorHelpers::phi(dir);
  m_payload.trueBound[eBoundTheta] = VectorHelpers::theta(dir);
  m_payload.trueBound[eBoundTime] = xt[eTime];

  m_payload.trueGx = xt[ePos0];
  m_payload.trueGy = xt[ePos1];
  m_payload.trueGz = xt[ePos2];

  m_payload.incidentPhi = angles.first;
  m_payload.incidentTheta = angles.second;
}

void Acts::RootMeasurementIo::fillBoundMeasurement(
    const std::vector<double>& measurements,
    const std::vector<double>& variances,
    const std::vector<unsigned int>& subspaceIndex) {
  if (measurements.size() != subspaceIndex.size() ||
      variances.size() != subspaceIndex.size()) {
    throw std::invalid_argument(
        "Measurement, covariance and access index size mismatch");
  }

  for (auto [im, m] : enumerate(measurements)) {
    auto ib = subspaceIndex[im];

    m_payload.recBound[ib] = m;
    m_payload.varBound[ib] = variances[im];
    m_payload.residual[ib] = m_payload.recBound[ib] - m_payload.trueBound[ib];
    m_payload.pull[ib] =
        m_payload.residual[ib] / std::sqrt(m_payload.varBound[ib]);
  }
}

void Acts::RootMeasurementIo::fillCluster(
    const Vector3& pos,
    const std::vector<std::tuple<int, int, float>>& channels) {
  m_payload.nch = static_cast<int>(channels.size());
  for (auto [ch0, ch1, chv] : channels) {
    m_payload.chId[0].push_back(static_cast<int>(ch0));
    m_payload.chId[1].push_back(static_cast<int>(ch1));
    m_payload.chValue.push_back(static_cast<float>(chv));
  }
  // Calculate cluster size in 0 and 1 direction
  auto [min0, max0] = std::ranges::minmax_element(m_payload.chId[0]);
  auto [min1, max1] = std::ranges::minmax_element(m_payload.chId[1]);
  m_payload.cSize[0] = (*max0 - *min0 + 1);
  m_payload.cSize[1] = (*max1 - *min1 + 1);
  m_payload.recGx = pos.x();
  m_payload.recGy = pos.y();
  m_payload.recGz = pos.z();
}

void Acts::RootMeasurementIo::clear() {
  for (unsigned int ib = 0; ib < eBoundSize; ++ib) {
    m_payload.trueBound[ib] = std::numeric_limits<float>::quiet_NaN();
    m_payload.recBound[ib] = std::numeric_limits<float>::quiet_NaN();
    m_payload.varBound[ib] = std::numeric_limits<float>::quiet_NaN();
    m_payload.residual[ib] = std::numeric_limits<float>::quiet_NaN();
    m_payload.pull[ib] = std::numeric_limits<float>::quiet_NaN();
  }
  m_payload.recGx = std::numeric_limits<float>::quiet_NaN();
  m_payload.recGy = std::numeric_limits<float>::quiet_NaN();
  m_payload.recGz = std::numeric_limits<float>::quiet_NaN();
  m_payload.trueGx = std::numeric_limits<float>::quiet_NaN();
  m_payload.trueGy = std::numeric_limits<float>::quiet_NaN();
  m_payload.trueGz = std::numeric_limits<float>::quiet_NaN();
  m_payload.incidentPhi = std::numeric_limits<float>::quiet_NaN();
  m_payload.incidentTheta = std::numeric_limits<float>::quiet_NaN();
  m_payload.nch = 0;
  m_payload.cSize[0] = 0;
  m_payload.cSize[1] = 0;
  m_payload.chId[0].clear();
  m_payload.chId[1].clear();
  m_payload.chValue.clear();
}
