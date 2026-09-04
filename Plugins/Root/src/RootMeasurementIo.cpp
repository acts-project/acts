// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootMeasurementIo.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <TTree.h>

using namespace Acts;

ActsPlugins::RootMeasurementIo::RootMeasurementIo(const Config& config)
    : m_cfg(config) {
  clear();
}

void ActsPlugins::RootMeasurementIo::connectForWrite(TTree& measurementTree) {
  measurementTree.Branch("event_nr", &m_measurementPayload.eventNr);
  measurementTree.Branch("measurement_id", &m_measurementPayload.measurementID);
  measurementTree.Branch("volume_id", &m_measurementPayload.volumeID);
  measurementTree.Branch("layer_id", &m_measurementPayload.layerID);
  measurementTree.Branch("surface_id", &m_measurementPayload.surfaceID);
  measurementTree.Branch("extra_id", &m_measurementPayload.extraID);

  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("rec_" + bNames[ib]).c_str(),
                           &m_measurementPayload.recBound[ib]);
  }
  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("var_" + bNames[ib]).c_str(),
                           &m_measurementPayload.varBound[ib]);
  }

  for (unsigned int ib = 0; ib < eBoundSize; ++ib) {
    measurementTree.Branch(("true_" + bNames[ib]).c_str(),
                           &m_measurementPayload.trueBound[ib]);
  }
  measurementTree.Branch("true_x", &m_measurementPayload.trueGx);
  measurementTree.Branch("true_y", &m_measurementPayload.trueGy);
  measurementTree.Branch("true_z", &m_measurementPayload.trueGz);
  measurementTree.Branch("true_incident_phi",
                         &m_measurementPayload.incidentPhi);
  measurementTree.Branch("true_incident_theta",
                         &m_measurementPayload.incidentTheta);

  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("residual_" + bNames[ib]).c_str(),
                           &m_measurementPayload.residual[ib]);
  }
  for (auto ib : m_cfg.recoIndices) {
    measurementTree.Branch(("pull_" + bNames[ib]).c_str(),
                           &m_measurementPayload.pull[ib]);
  }
  clear();
}

void ActsPlugins::RootMeasurementIo::fillIdentification(
    int evnt, std::uint64_t measurementId, const GeometryIdentifier& geoId) {
  m_measurementPayload.eventNr = evnt;
  m_measurementPayload.measurementID = measurementId;
  m_measurementPayload.volumeID = static_cast<int>(geoId.volume());
  m_measurementPayload.layerID = static_cast<int>(geoId.layer());
  m_measurementPayload.surfaceID = static_cast<int>(geoId.sensitive());
  m_measurementPayload.extraID = static_cast<int>(geoId.extra());
}

void ActsPlugins::RootMeasurementIo::fillTruthParameters(
    const Vector2& lp, const Vector4& xt, const Vector3& dir,
    const std::pair<double, double> angles) {
  m_measurementPayload.trueBound[eBoundLoc0] = lp[eBoundLoc0];
  m_measurementPayload.trueBound[eBoundLoc1] = lp[eBoundLoc1];
  m_measurementPayload.trueBound[eBoundPhi] = VectorHelpers::phi(dir);
  m_measurementPayload.trueBound[eBoundTheta] = VectorHelpers::theta(dir);
  m_measurementPayload.trueBound[eBoundTime] = xt[eTime];

  m_measurementPayload.trueGx = xt[ePos0];
  m_measurementPayload.trueGy = xt[ePos1];
  m_measurementPayload.trueGz = xt[ePos2];

  m_measurementPayload.incidentPhi = static_cast<float>(angles.first);
  m_measurementPayload.incidentTheta = static_cast<float>(angles.second);
}

void ActsPlugins::RootMeasurementIo::fillBoundMeasurement(
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

    m_measurementPayload.recBound[ib] = static_cast<float>(m);
    m_measurementPayload.varBound[ib] = static_cast<float>(variances[im]);
    m_measurementPayload.residual[ib] =
        m_measurementPayload.recBound[ib] - m_measurementPayload.trueBound[ib];
    m_measurementPayload.pull[ib] =
        m_measurementPayload.residual[ib] /
        std::sqrt(m_measurementPayload.varBound[ib]);
  }
}

void ActsPlugins::RootMeasurementIo::clear() {
  for (unsigned int ib = 0; ib < eBoundSize; ++ib) {
    m_measurementPayload.trueBound[ib] =
        std::numeric_limits<float>::quiet_NaN();
    m_measurementPayload.recBound[ib] = std::numeric_limits<float>::quiet_NaN();
    m_measurementPayload.varBound[ib] = std::numeric_limits<float>::quiet_NaN();
    m_measurementPayload.residual[ib] = std::numeric_limits<float>::quiet_NaN();
    m_measurementPayload.pull[ib] = std::numeric_limits<float>::quiet_NaN();
  }
  m_measurementPayload.trueGx = std::numeric_limits<float>::quiet_NaN();
  m_measurementPayload.trueGy = std::numeric_limits<float>::quiet_NaN();
  m_measurementPayload.trueGz = std::numeric_limits<float>::quiet_NaN();
  m_measurementPayload.incidentPhi = std::numeric_limits<float>::quiet_NaN();
  m_measurementPayload.incidentTheta = std::numeric_limits<float>::quiet_NaN();
}
