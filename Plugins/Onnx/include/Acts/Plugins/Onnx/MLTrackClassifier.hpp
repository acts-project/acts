// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <core/session/onnxruntime_cxx_api.h>

namespace Acts {

// Classify tracks as good/duplicate/fake using using a deep neural network.
class MLTrackClassifier {
 public:
  /// @brief The labels for track quality
  enum TrackLabels { good, duplicate, fake };

  /// @brief Default constructor
  MLTrackClassifier() = default;

  /// @brief Parametrized constructor
  ///
  /// @param env the ONNX runtime environment
  /// @param modelPath the path to the ML model in *.onnx format
  MLTrackClassifier(Ort::Env& env, const char* modelPath);

  /// @brief Default destructor
  ~MLTrackClassifier() = default;

  /// @brief Predict the track label
  ///
  /// @param inputFeatures The vector of input features for the trajectory to be
  /// classified
  /// @param decisionThreshProb The probability threshold used to predict the
  /// track label
  ///
  /// @return The predicted track label of the trajectory
  TrackLabels predictTrackLabel(std::vector<float>& inputFeatures,
                                const double& decisionThreshProb) const;

  /// @brief Run the ONNX inference function
  ///
  /// @param inputTensorValues The input feature values used for prediction
  ///
  /// @return Pointer to output values
  float* runONNXInference(std::vector<float>& inputTensorValues) const;

 private:
  /// ONNX runtime session / model properties
  std::unique_ptr<Ort::Session> m_session;
  std::vector<const char*> m_inputNodeNames;
  std::vector<int64_t> m_inputNodeDims;
  std::vector<const char*> m_outputNodeNames;
};

}  // namespace Acts