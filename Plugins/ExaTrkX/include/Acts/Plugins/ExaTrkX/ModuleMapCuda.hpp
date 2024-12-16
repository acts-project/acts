// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/Stages.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

template <typename T>
class CUDA_module_map_doublet;

template <typename T>
class CUDA_module_map_triplet;

namespace Acts {

namespace detail {
class GraphCreatorWrapperBase;

template <typename T>
struct CUDA_hit_data {
  std::size_t m_size;
  std::uint64_t *m_cuda_hit_id;
  T *m_cuda_x;
  T *m_cuda_y;
  T *m_cuda_R;
  T *m_cuda_phi;
  T *m_cuda_z;
  T *m_cuda_eta;

  std::size_t size() { return m_size; }
  T *cuda_x() { return m_cuda_x; }
  T *cuda_y() { return m_cuda_y; }
  T *cuda_z() { return m_cuda_z; }
  T *cuda_R() { return m_cuda_R; }
  T *cuda_phi() { return m_cuda_phi; }
  T *cuda_eta() { return m_cuda_eta; }
  std::uint64_t *cuda_hit_id() { return m_cuda_hit_id; }
};

template <typename T>
struct CUDA_edge_data {
  std::size_t size;

  int *cuda_graph_M1_hits;
  int *cuda_graph_M2_hits;

  T *cuda_graph_dR;
  T *cuda_graph_dz;
  T *cuda_graph_deta;
  T *cuda_graph_dphi;
  T *cuda_graph_phi_slope;
  T *cuda_graph_r_phi_slope;
};

}  // namespace detail

class ModuleMapCuda : public GraphConstructionBase {
 public:
  struct Config {
    std::string moduleMapPath;
    float rScale = 1.0;
    float phiScale = 1.0;
    float zScale = 1.0;
    float etaScale = 1.0;

    int maxEdgesAllocate = 3000;

    int gpuDevice = 0;
    int gpuBlocks = 512;
  };

 private:
  std::pair<Acts::detail::CUDA_edge_data<float>,
            Acts::detail::CUDA_hit_data<float>>
  makeEdges(Acts::detail::CUDA_hit_data<float> cuda_TThits,
            int *cuda_hit_indice) const;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::unique_ptr<CUDA_module_map_doublet<float>> m_cudaModuleMapDoublet;
  std::unique_ptr<CUDA_module_map_triplet<float>> m_cudaModuleMapTriplet;

  const auto &logger() const { return *m_logger; }

 public:
  ModuleMapCuda(const Config &cfg, std::unique_ptr<const Acts::Logger> logger);
  ~ModuleMapCuda();

  const auto &config() const { return m_cfg; }

  std::tuple<std::any, std::any, std::any> operator()(
      std::vector<float> &inputValues, std::size_t numNodes,
      const std::vector<std::uint64_t> &moduleIds,
      torch::Device device = torch::Device(torch::kCPU)) override;

  // TODO this returns nothing useful
  torch::Device device() const override { return torch::Device(torch::kCPU); }
};

}  // namespace Acts
