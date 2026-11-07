// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/detail/TorchDevice.hpp"

#include <stdexcept>
#include <string>

namespace ActsPlugins::detail {

torch::Device toTorchDevice(const Device &device) {
  switch (device.type) {
    case Device::Type::eCPU:
      return torch::Device(torch::kCPU);
    case Device::Type::eCUDA: {
      if (!torch::cuda::is_available()) {
        throw std::runtime_error(
            "CUDA device requested but CUDA is not available");
      }
      if (device.index >=
          static_cast<std::size_t>(torch::cuda::device_count())) {
        throw std::runtime_error(
            "CUDA device index " + std::to_string(device.index) +
            " is out of range (" + std::to_string(torch::cuda::device_count()) +
            " devices available)");
      }
      return torch::Device(torch::kCUDA,
                           static_cast<c10::DeviceIndex>(device.index));
    }
    case Device::Type::eMPS: {
      if (!at::hasMPS()) {
        throw std::runtime_error(
            "MPS device requested but MPS is not available (needs a torch "
            "build with MPS support, running on Apple silicon)");
      }
      return torch::Device(torch::kMPS);
    }
  }
  throw std::invalid_argument("Unknown device type");
}

std::unique_ptr<torch::jit::Module> loadTorchModel(
    const std::string &modelPath, const Device &device,
    const std::optional<Device> &modelDevice) {
  if (device.isMps()) {
    throw std::invalid_argument(
        "MPS is not a device the pipeline tensors can live on: configure it as "
        "modelDevice instead of device");
  }
  const torch::Device loadDevice =
      modelDevice.has_value() ? toTorchDevice(*modelDevice)
                              : (device.isCuda() ? toTorchDevice(device)
                                                 : torch::Device(torch::kCPU));

  try {
    auto model = std::make_unique<torch::jit::Module>(
        torch::jit::load(modelPath, loadDevice));
    model->eval();
    return model;
  } catch (const c10::Error &e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

}  // namespace ActsPlugins::detail
