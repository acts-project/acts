// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Onnx/OnnxRuntimeBase.hpp"

#include <cassert>
#include <stdexcept>

// parametrized constructor
Acts::OnnxRuntimeBase::OnnxRuntimeBase(Ort::Env& env, const char* modelPath) {
  // set the ONNX runtime session options
  Ort::SessionOptions sessionOptions;
  // set graph optimization level
  sessionOptions.SetGraphOptimizationLevel(
      GraphOptimizationLevel::ORT_ENABLE_BASIC);
  // create the Ort session
  m_session = std::make_unique<Ort::Session>(env, modelPath, sessionOptions);

  // default allocator
  Ort::AllocatorWithDefaultOptions allocator;

  // get the names of the input nodes of the model
  size_t numInputNodes = m_session->GetInputCount();
  m_inputNodeNames.resize(numInputNodes);

  // iterate over all input nodes and get the name
  for (size_t i = 0; i < numInputNodes; i++) {
    m_inputNodeNames[i] = m_session->GetInputName(i, allocator);

    // get the dimensions of the input nodes
    // here we assume that all input nodes have the dimensions
    Ort::TypeInfo inputTypeInfo = m_session->GetInputTypeInfo(i);
    auto tensorInfo = inputTypeInfo.GetTensorTypeAndShapeInfo();
    m_inputNodeDims = tensorInfo.GetShape();
    // fix for symbolic dim = -1 from python
    for (size_t j = 0; j < m_inputNodeDims.size(); j++) {
      if (m_inputNodeDims[j] < 0) {
        m_inputNodeDims[j] = 1;
      }
    }
  }

  // get the names of the output nodes
  size_t numOutputNodes = m_session->GetOutputCount();
  m_outputNodeNames.resize(numOutputNodes);

  // iterate over all output nodes and get the name
  for (size_t i = 0; i < numOutputNodes; i++) {
    m_outputNodeNames[i] = m_session->GetOutputName(i, allocator);

    // get the dimensions of the output nodes
    // here we assume that all output nodes have the dimensions
    Ort::TypeInfo outputTypeInfo = m_session->GetOutputTypeInfo(i);
    auto tensorInfo = outputTypeInfo.GetTensorTypeAndShapeInfo();
    m_outputNodeDims = tensorInfo.GetShape();
    // fix for symbolic dim = -1 from python
    for (size_t j = 0; j < m_outputNodeDims.size(); j++) {
      if (m_outputNodeDims[j] < 0) {
        m_outputNodeDims[j] = 1;
      }
    }
  }
}

// inference function using ONNX runtime
// the function assumes that the model has 1 input node and 1 output node
std::vector<float> Acts::OnnxRuntimeBase::runONNXInference(
    std::vector<float>& inputTensorValues) const {
  // create input tensor object from data values
  // note: this assumes the model has only 1 input node
  Ort::MemoryInfo memoryInfo =
      Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  Ort::Value inputTensor = Ort::Value::CreateTensor<float>(
      memoryInfo, inputTensorValues.data(), inputTensorValues.size(),
      m_inputNodeDims.data(), m_inputNodeDims.size());
  // double-check that inputTensor is a Tensor
  if (!inputTensor.IsTensor()) {
    throw std::runtime_error(
        "runONNXInference: conversion of input to Tensor failed. ");
  }

  // score model on input tensors, get back output tensors
  std::vector<Ort::Value> outputTensors =
      m_session->Run(Ort::RunOptions{nullptr}, m_inputNodeNames.data(),
                     &inputTensor, m_inputNodeNames.size(),
                     m_outputNodeNames.data(), m_outputNodeNames.size());
  // double-check that outputTensors contains Tensors and that the count matches
  // that of output nodes
  if (!outputTensors[0].IsTensor() ||
      (outputTensors.size() != m_outputNodeNames.size())) {
    throw std::runtime_error(
        "runONNXInference: calculation of output failed. ");
  }

  // get pointer to output tensor float values
  // note: this assumes the model has only 1 output node
  float* outputTensor = outputTensors.front().GetTensorMutableData<float>();

  // get the output values
  std::vector<float> outputTensorValues(m_outputNodeDims[1]);
  for (size_t i = 0; i < outputTensorValues.size(); i++) {
    outputTensorValues[i] = outputTensor[i];
  }
  return outputTensorValues;
}
