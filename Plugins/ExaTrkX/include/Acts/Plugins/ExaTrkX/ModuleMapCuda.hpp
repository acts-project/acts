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

struct CUstream_st;
typedef CUstream_st *cudaStream_t;

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
  std::size_t nEdges;
  int *cudaEdgePtr;
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

    bool moreParallel = true;
    int gpuDevice = 0;
    int gpuBlocks = 512;

    float epsilon = 1e-8f;
  };

 private:
  detail::CUDA_edge_data<float> makeEdges(
      detail::CUDA_hit_data<float> cuda_TThits, int *cuda_hit_indice,
      cudaStream_t &stream) const;

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  std::unique_ptr<CUDA_module_map_doublet<float>> m_cudaModuleMapDoublet;
  std::unique_ptr<CUDA_module_map_triplet<float>> m_cudaModuleMapTriplet;

  // TODO make this a managed storage soon
  std::uint64_t *m_cudaModuleMapKeys{};
  int *m_cudaModuleMapVals{};
  std::size_t m_cudaModuleMapSize{};

  const auto &logger() const { return *m_logger; }

 public:
  ModuleMapCuda(const Config &cfg, std::unique_ptr<const Acts::Logger> logger);
  ~ModuleMapCuda() override;

  const auto &config() const { return m_cfg; }

  PipelineTensors operator()(std::vector<float> &inputValues,
                             std::size_t numNodes,
                             const std::vector<std::uint64_t> &moduleIds,
                             const ExecutionContext &execContext = {}) override;
};

}  // namespace Acts
