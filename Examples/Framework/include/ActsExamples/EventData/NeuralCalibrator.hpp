// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Plugins/Onnx/OnnxRuntimeBase.hpp>
#include <ActsExamples/EventData/MeasurementCalibration.hpp>

#include <filesystem>

namespace ActsExamples {

class NeuralCalibrator : public MeasurementCalibrator {
 public:
  size_t n_inputs = 57;  // TODO: set this dynamically

  NeuralCalibrator(const std::filesystem::path& modelPath);

  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const override;

  bool needsClusters() const override { return true; }

 private:
  Ort::Env m_env;
  Acts::OnnxRuntimeBase m_model;
};

}  // namespace ActsExamples
