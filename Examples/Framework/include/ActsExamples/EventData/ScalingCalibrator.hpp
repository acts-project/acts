// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// TODO move this somewhere else?

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <ActsExamples/EventData/MeasurementCalibration.hpp>

#include <TFile.h>

namespace ActsExamples {

template <typename indices_t, size_t kSize>
constexpr size_t measurementSize(
    const Acts::Measurement<indices_t, kSize>& meas) {
  return kSize;
}

class ScalingCalibrator : public MeasurementCalibrator {
 public:
  ScalingCalibrator(const char* path);

  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const;

  bool needsClusters() { return true; }

 private:
  using GeoId = Acts::GeometryIdentifier::Value;
  using Constant = Double_t;

  struct ConstantTuple {
    Constant x_offset;
    Constant x_scale;
    Constant y_offset;
    Constant y_scale;
  };

  GeoId m_geoId_mask{0};
  std::unordered_map<GeoId, ConstantTuple> m_offsetScaleMap;

  ConstantTuple getOffsetAndScale(Acts::GeometryIdentifier geoId) const;

  template <typename T>
  std::vector<T>* readVec(TFile& tf, const char* key, size_t size_min,
                          size_t size_max) const;

  void readMap(const char* path);
};

}  // namespace ActsExamples
