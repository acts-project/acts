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

// General class that sets up the ONNX runtime framework for loading a ML model
// and using it for inference.
class OnnxRuntimeBase {
 public:
  /// @brief Default constructor
  OnnxRuntimeBase() = default;

  /// @brief Parametrized constructor
  ///
  /// @param env the ONNX runtime environment
  /// @param modelPath the path to the ML model in *.onnx format
  OnnxRuntimeBase(Ort::Env& env, const char* modelPath);

  /// @brief Default destructor
  ~OnnxRuntimeBase() = default;

 protected:
  /// @brief Run the ONNX inference function
  ///
  /// @param inputTensorValues The input feature values used for prediction
  ///
  /// @return The output (predicted) values
  std::vector<float> runONNXInference(
      std::vector<float>& inputTensorValues) const;

 private:
  /// ONNX runtime session / model properties
  std::unique_ptr<Ort::Session> m_session;
  std::vector<const char*> m_inputNodeNames;
  std::vector<int64_t> m_inputNodeDims;
  std::vector<const char*> m_outputNodeNames;
  std::vector<int64_t> m_outputNodeDims;
};

}  // namespace Acts