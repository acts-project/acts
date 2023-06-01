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

// FIXME: Some concept of geometry (i.e. this is only for pixel!)
// Temp fix: pass a list of blessed volumes, defaulting to 7/8/9
// In subsequent MR,
// should move to GeometryHierarchyMap<MeasurementCalibrator> interface

class NeuralCalibrator : public MeasurementCalibrator {
 public:
  /// Measurement position calibration based on mixture density network
  /// (MDN) model. The model takes as input:
  ///
  /// - A 7x7 charge matrix centered on the center pixel of the cluster;
  /// - The volume and layer identifiers from
  ///   the GeometryIdentifier of the containing surface;
  /// - The bound phi and theta angles of the predicted track state;
  /// - The initial estimated position
  /// - The initial estimated variance
  ///
  /// Given these inputs, a mixture density network estimates
  /// the parameters of a gaussian mixture model:
  ///
  /// P(Y|X) = \sum_i P(Prior_i) N(Y|Mean_i(X), Variance_i(X))
  ///
  /// These are translated to single position + variance estimate by
  /// taking the most probable value based on the estimated priors.
  /// The measurements are assumed to be 2-dimensional.
  ///
  /// This class implements the MeasurementCalibrator interface, and
  /// therefore internally computes the network input and runs the
  /// inference engine itself.
  ///
  /// @param [in] modelPath The path to the .onnx model file
  /// @param [in] nComp The number of components in the gaussian mixture
  NeuralCalibrator(const std::filesystem::path& modelPath, size_t nComp = 1);

  /// The MeasurementCalibrator interface methods
  void calibrate(
      const MeasurementContainer& measurements,
      const ClusterContainer* clusters, const Acts::GeometryContext& /*gctx*/,
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy&
          trackState) const override;

  bool needsClusters() const override { return true; }

 private:
  Ort::Env m_env;
  Acts::OnnxRuntimeBase m_model;
  size_t m_nComp;
  size_t n_inputs =
      57;  // TODO make this configurable? e.g. for changing matrix size?
};

}  // namespace ActsExamples
