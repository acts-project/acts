// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <core/session/onnxruntime_cxx_api.h>

namespace Acts {

/// Class implementing the Exa.TrkX track finding algorithm.
/// It holds the required ONNX objects.
class ExaTrkXTrackFinding {
 public:
  /// Configuration struct for the track finding.
  struct Config {
    // input model directory
    std::string inputMLModuleDir;

    // hyperparameters in the pipeline.
    int64_t spacepointFeatures = 3;
    int embeddingDim = 8;
    float rVal = 1.6;
    int knnVal = 500;
    float filterCut = 0.21;
  };

  /// Constructor of the track finding module
  ///
  /// @param config is the config struct to configure the module
  ExaTrkXTrackFinding(const Config& config);

  virtual ~ExaTrkXTrackFinding() = default;

  /// Do the track finding
  ///
  /// @param [in] input_values Packed spacepoints in the form
  /// [ r1, phi1, z1, r2, phi2, z2, ... ]
  /// @param [in] spacepointIDs corresponding spacepoint IDs to the input_values.
  /// @param [out] trackCandidates nested vector containing the spacepoint ids
  /// of the found tracks
  /// @note The input values are not const, because the underlying ONNX API
  /// takes only non-const pointers.
  void getTracks(std::vector<float>& input_values,
                 std::vector<uint32_t>& spacepointIDs,
                 std::vector<std::vector<uint32_t> >& trackCandidates) const;

  /// Return the configuration object of the track finding module
  const Config& config() const { return m_cfg; }

 private:
  void runSessionWithIoBinding(Ort::Session& sess,
                               std::vector<const char*>& inputNames,
                               std::vector<Ort::Value>& inputData,
                               std::vector<const char*>& outputNames,
                               std::vector<Ort::Value>& outputData) const;

  void buildEdges(std::vector<float>& embedFeatures,
                  std::vector<int64_t>& edgeList, int64_t numSpacepoints) const;

 private:
  Config m_cfg;
  std::unique_ptr<Ort::Env> m_env;
  std::unique_ptr<Ort::Session> e_sess;
  std::unique_ptr<Ort::Session> f_sess;
  std::unique_ptr<Ort::Session> g_sess;
};

}  // namespace Acts
