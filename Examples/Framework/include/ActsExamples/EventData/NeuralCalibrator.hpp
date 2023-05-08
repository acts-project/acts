// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <ActsExamples/EventData/MeasurementCalibration.hpp>
#include <Acts/Plugins/Onnx/OnnxRuntimeBase.hpp>

#include <filesystem>

namespace ActsExamples {

class MDN_Model : public Acts::OnnxRuntimeBase {
public:
  MDN_Model(const std::filesystem::path& path);
  std::tuple<float, float, float, float> evaluate(
      const std::vector<float>& input) const;
};

class NeuralCalibrator : public MeasurementCalibrator {
 public:
  NeuralCalibrator(const std::filesystem::path& path);

  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const override;

  bool needsClusters() const override { return true; }

 private:
  MDN_Model m_mdn;
};

}  // namespace ActsExamples
