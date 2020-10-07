// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Performance/CKFPerformanceWriterOptions.hpp"

#include <string>
#include <boost/program_options.hpp>

void ActsExamples::Options::addCKFPerformanceWriterOptions(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::value;

  // construct path to demo model.onnx from current file path
  std::string currentFilePath = __FILE__;
  std::string delimiter = "Examples/";
  std::string commonPath =
      currentFilePath.substr(0, currentFilePath.find(delimiter));
  std::string demoModelPath =
      commonPath + "Examples/Run/Reconstruction/MLAmbiguityResolutionDemo.onnx";

  auto opt = desc.add_options();
  opt("use-ml-track-classification", value<bool>()->default_value(false),
      "Use trained neural network model to classify reconstructed tracks as "
      "good/duplicate/fake");
  opt("ml-model-path", value<std::string>()->default_value(demoModelPath),
      "Path to the trained neural network model (in onnx format)");
}

ActsExamples::CKFPerformanceWriter::Config
ActsExamples::Options::readCKFPerformanceWriterConfig(
    const ActsExamples::Options::Variables& variables) {
  using Config = typename ActsExamples::CKFPerformanceWriter::Config;
  Config ckfPerfWriterCfg;

  ckfPerfWriterCfg.useMLTrackClassifier =
      variables["use-ml-track-classification"].template as<bool>();
  ckfPerfWriterCfg.onnxModelFilename =
      variables["ml-model-path"].template as<std::string>();

  return ckfPerfWriterCfg;
}