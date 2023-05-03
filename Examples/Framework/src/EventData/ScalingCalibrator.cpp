// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Definitions/TrackParametrization.hpp>
#include <ActsExamples/EventData/ScalingCalibrator.hpp>

ActsExamples::ScalingCalibrator::ScalingCalibrator(const char* path) {
  readMap(path);
}

void ActsExamples::ScalingCalibrator::calibrate(
    const MeasurementContainer& measurements,
    const ClusterContainer* /*clusters*/, const Acts::GeometryContext& /*gctx*/,
    Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
        trackState) const {
  const IndexSourceLink& sourceLink =
      trackState.getUncalibratedSourceLink().get<IndexSourceLink>();
  assert((sourceLink.index() < measurements.size()) and
         "Source link index is outside the container bounds");

  auto [xoff, xscale, yoff, yscale] =
      getOffsetAndScale(sourceLink.geometryId());

  std::visit(
      [&](const auto& meas) {
        auto E = meas.expander();
        auto P = meas.projector();

        Acts::ActsVector<Acts::eBoundSize> fpar = E * meas.parameters();

        Acts::ActsSymMatrix<Acts::eBoundSize> fcov =
            E * meas.covariance() * E.transpose();

        fpar[Acts::eBoundLoc0] += xoff;
        fpar[Acts::eBoundLoc1] += yoff;
        fcov(Acts::eBoundLoc0, Acts::eBoundLoc0) *= xscale;
        fcov(Acts::eBoundLoc1, Acts::eBoundLoc1) *= yscale;

        constexpr size_t kSize = meas.size();
        std::array<Acts::BoundIndices, kSize> indices = meas.indices();

        Acts::ActsVector<kSize> cpar = P * fpar;
        Acts::ActsSymMatrix<kSize> ccov = P * fcov * P.transpose();

        Acts::SourceLink sl{sourceLink.geometryId(), sourceLink};

        Acts::Measurement<Acts::BoundIndices, kSize> cmeas(std::move(sl),
                                                           indices, cpar, ccov);

        trackState.allocateCalibrated(cmeas.size());
        trackState.setCalibrated(cmeas);
      },
      (measurements)[sourceLink.index()]);
}

ActsExamples::ScalingCalibrator::ConstantTuple
ActsExamples::ScalingCalibrator::getOffsetAndScale(
    Acts::GeometryIdentifier geoId) const {
  // FIXME
  // return m_offsetScaleMap.at(geoId.getBits(m_geoId_mask));
  return m_offsetScaleMap.at(0);
}

template <typename T>
std::vector<T>* ActsExamples::ScalingCalibrator::readVec(
    TFile& tf, const char* key, size_t size_min, size_t size_max) const {
  std::vector<T>* vec;
  tf.GetObject(key, vec);

  if (!vec) {
    throw std::runtime_error(std::string("Vector \"") + std::string(key) +
                             std::string("\" not found!"));
  }
  if (vec->size() < size_min || vec->size() > size_max) {
    throw std::runtime_error(std::string("Vector \"") + std::string(key) +
                             std::string("\" has invalid size!"));
  }
  return vec;
}

void ActsExamples::ScalingCalibrator::readMap(const char* path) {
  TFile tf(path, "READ");
  if (tf.IsZombie()) {
    throw std::runtime_error(std::string("Unable to open TFile: ") + path);
  }

  std::vector<GeoId>* v_mask = readVec<GeoId>(tf, "v_mask", 1, 1);
  m_geoId_mask = v_mask->at(0);

  std::vector<GeoId>* v_geoId =
      readVec<GeoId>(tf, "v_geoId", 1, std::numeric_limits<size_t>::max());

  std::vector<Constant>* v_x_offset =
      readVec<Constant>(tf, "v_x_offset", v_geoId->size(), v_geoId->size());

  std::vector<Constant>* v_x_scale =
      readVec<Constant>(tf, "v_x_scale", v_geoId->size(), v_geoId->size());

  std::vector<Constant>* v_y_offset =
      readVec<Constant>(tf, "v_y_offset", v_geoId->size(), v_geoId->size());

  std::vector<Constant>* v_y_scale =
      readVec<Constant>(tf, "v_y_scale", v_geoId->size(), v_geoId->size());

  for (size_t i = 0; i < v_geoId->size(); i++) {
    m_offsetScaleMap.insert({v_geoId->at(i),
                             {v_x_offset->at(i), v_x_scale->at(i),
                              v_y_offset->at(i), v_y_scale->at(i)}});
  }

  // Is this needed?
  delete v_mask;
  delete v_geoId;
  delete v_x_offset;
  delete v_x_scale;
  delete v_y_offset;
  delete v_y_scale;
}
