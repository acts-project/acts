// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Gnn/Tensor.hpp"

#include <memory>
#include <optional>
#include <string>

#include <torch/script.h>
#include <torch/torch.h>

namespace ActsPlugins::detail {

/// Translate a device description into the torch device the models run on.
///
/// The backend is checked here so that a device that torch cannot use is
/// reported when a stage is constructed, rather than silently running
/// somewhere else or failing at the first inference.
///
/// @param device the device to translate
/// @return the equivalent torch device
/// @throws std::runtime_error if the backend is not available
/// @throws std::invalid_argument if the device type is unknown
torch::Device toTorchDevice(const Device &device);

/// Load a TorchScript model for a pipeline stage.
///
/// The model is loaded on @p modelDevice if set; otherwise on @p device if it
/// is a CUDA device, else on the CPU (the historical behaviour).
///
/// @param modelPath path to the TorchScript file
/// @param device the device the pipeline tensors live on
/// @param modelDevice the optional device the model runs on
/// @return the loaded model in eval mode
/// @throws std::invalid_argument if @p device is MPS or the model fails to load
std::unique_ptr<torch::jit::Module> loadTorchModel(
    const std::string &modelPath, const Device &device,
    const std::optional<Device> &modelDevice);

}  // namespace ActsPlugins::detail
