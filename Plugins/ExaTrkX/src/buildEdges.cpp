// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "buildEdges.hpp"

#include <iostream>
#include <mutex>
#include <vector>

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <grid/counting_sort.h>
#include <grid/find_nbrs.h>
#include <grid/grid.h>
#include <grid/insert_points.h>
#include <grid/prefix_sum.h>
#include <torch/script.h>
#include <torch/torch.h>

torch::Tensor Acts::buildEdges(at::Tensor &embedFeatures,
                               int64_t numSpacepoints, int dim, float rVal,
                               int kVal, bool flipDirections) {
  using namespace torch::indexing;

  torch::Device device(torch::kCUDA);

  const int grid_params_size = 8;
  const int grid_delta_idx = 3;
  const int grid_total_idx = 7;
  const int grid_max_res = 128;
  const int grid_dim = 3;

  if (dim < 3) {
    throw std::runtime_error("DIM < 3 is not supported for now.\n");
  }

  const float radius_cell_ratio = 2.0;
  const int batch_size = 1;
  int G = -1;

  // Set up grid properties
  torch::Tensor grid_min;
  torch::Tensor grid_max;
  torch::Tensor grid_size;

  torch::Tensor embedTensor = embedFeatures.reshape({1, numSpacepoints, dim});
  torch::Tensor gridParamsCuda =
      torch::zeros({batch_size, grid_params_size}, device).to(torch::kFloat32);
  torch::Tensor r_tensor = torch::full({batch_size}, rVal, device);
  torch::Tensor lengths = torch::full({batch_size}, numSpacepoints, device);

  // build the grid
  for (int i = 0; i < batch_size; i++) {
    torch::Tensor allPoints =
        embedTensor.index({i, Slice(None, lengths.index({i}).item().to<long>()),
                           Slice(None, grid_dim)});
    grid_min = std::get<0>(allPoints.min(0));
    grid_max = std::get<0>(allPoints.max(0));
    gridParamsCuda.index_put_({i, Slice(None, grid_delta_idx)}, grid_min);

    grid_size = grid_max - grid_min;

    float cell_size =
        r_tensor.index({i}).item().to<float>() / radius_cell_ratio;

    if (cell_size < (grid_size.min().item().to<float>() / grid_max_res)) {
      cell_size = grid_size.min().item().to<float>() / grid_max_res;
    }

    gridParamsCuda.index_put_({i, grid_delta_idx}, 1 / cell_size);

    gridParamsCuda.index_put_({i, Slice(1 + grid_delta_idx, grid_total_idx)},
                              floor(grid_size / cell_size) + 1);

    gridParamsCuda.index_put_(
        {i, grid_total_idx},
        gridParamsCuda.index({i, Slice(1 + grid_delta_idx, grid_total_idx)})
            .prod());

    if (G < gridParamsCuda.index({i, grid_total_idx}).item().to<int>()) {
      G = gridParamsCuda.index({i, grid_total_idx}).item().to<int>();
    }
  }

  torch::Tensor pc_grid_cnt =
      torch::zeros({batch_size, G}, device).to(torch::kInt32);
  torch::Tensor pc_grid_cell =
      torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);
  torch::Tensor pc_grid_idx =
      torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);

  // put spacepoints into the grid
  InsertPointsCUDA(embedTensor, lengths.to(torch::kInt64), gridParamsCuda,
                   pc_grid_cnt, pc_grid_cell, pc_grid_idx, G);

  torch::Tensor pc_grid_off =
      torch::full({batch_size, G}, 0, device).to(torch::kInt32);
  torch::Tensor grid_params = gridParamsCuda.to(torch::kCPU);

  // for loop seems not to be necessary anymore
  pc_grid_off = PrefixSumCUDA(pc_grid_cnt, grid_params);

  torch::Tensor sorted_points =
      torch::zeros({batch_size, numSpacepoints, dim}, device)
          .to(torch::kFloat32);
  torch::Tensor sorted_points_idxs =
      torch::full({batch_size, numSpacepoints}, -1, device).to(torch::kInt32);

  CountingSortCUDA(embedTensor, lengths.to(torch::kInt64), pc_grid_cell,
                   pc_grid_idx, pc_grid_off, sorted_points, sorted_points_idxs);

  std::tuple<at::Tensor, at::Tensor> nbr_output = FindNbrsCUDA(
      sorted_points, sorted_points, lengths.to(torch::kInt64),
      lengths.to(torch::kInt64), pc_grid_off.to(torch::kInt32),
      sorted_points_idxs, sorted_points_idxs,
      gridParamsCuda.to(torch::kFloat32), kVal, r_tensor, r_tensor * r_tensor);
  torch::Tensor positiveIndices = std::get<0>(nbr_output) >= 0;

  torch::Tensor repeatRange = torch::arange(positiveIndices.size(1), device)
                                  .repeat({1, positiveIndices.size(2), 1})
                                  .transpose(1, 2);

  torch::Tensor stackedEdges =
      torch::stack({repeatRange.index({positiveIndices}),
                    std::get<0>(nbr_output).index({positiveIndices})});

  // Remove self-loops
  torch::Tensor selfLoopMask =
      stackedEdges.index({0}) != stackedEdges.index({1});
  stackedEdges = stackedEdges.index({Slice(), selfLoopMask});

  // Remove duplicates
  torch::Tensor duplicate_mask =
      stackedEdges.index({0}) > stackedEdges.index({1});
  stackedEdges = stackedEdges.index({Slice(), duplicate_mask});

  // Randomly flip direction
  if (flipDirections) {
    torch::Tensor random_cut_keep = torch::randint(2, {stackedEdges.size(1)});
    torch::Tensor random_cut_flip = 1 - random_cut_keep;
    torch::Tensor keep_edges =
        stackedEdges.index({Slice(), random_cut_keep.to(torch::kBool)});
    torch::Tensor flip_edges =
        stackedEdges.index({Slice(), random_cut_flip.to(torch::kBool)})
            .flip({0});
    stackedEdges = torch::cat({keep_edges, flip_edges}, 1);
  } else {
    stackedEdges = stackedEdges.toType(torch::kInt64).to(torch::kCPU);
  }

  return stackedEdges;
}

torch::Tensor buildEdgesBruteForce(at::Tensor &embedFeatures,
                                   int64_t numSpacepoints, int dim, float rVal,
                                   int) {
  auto distance = [](const at::Tensor &a, const at::Tensor &b) {
    return std::sqrt(((a - b) * (a - b)).sum().item().to<float>());
  };

  struct TwoRanges {
    int i_min, i_max, j_min, j_max;
  };

  torch::Tensor data =
      embedFeatures.reshape({numSpacepoints, dim}).to(torch::kCPU);
  std::cout << "data: " << data.size(0) << " " << data.size(1) << "\n";

  const int n_chunks = std::thread::hardware_concurrency();
  const int per_chunk = numSpacepoints / n_chunks;

  std::vector<TwoRanges> ranges;

  for (int i = 0; i < n_chunks; ++i) {
    for (int j = i; j < n_chunks; ++j) {
      TwoRanges r;

      r.i_min = i * per_chunk;

      if (i < n_chunks - 1) {
        r.i_max = (i + 1) * per_chunk;
      } else {
        r.i_max = numSpacepoints;
      }

      r.j_min = j * per_chunk;

      if (j < n_chunks - 1) {
        r.j_max = (j + 1) * per_chunk;
      } else {
        r.j_max = numSpacepoints;
      }

      ranges.push_back(r);
    }
  }
  std::cout << "#ranges: " << ranges.size() << "\n";

  std::vector<int> all_edges;
  all_edges.reserve(2 * numSpacepoints * 10);

  std::mutex res_mutex;
  std::mutex int_mutex;
  std::mutex progress_mutex;
  std::vector<float> progress(ranges.size());
  int index = -1;
  int print_id = 1;

  // This for loop can be easily parallelized, however this is not done here to
  // not mess up the dependencies. Possible options are:
  // clang-format off
  // std::for_each(std::execution::par, ranges.begin(), ranges.end(), [&](const TwoRanges &r) { 
  // tbb::parallel_for_each(ranges.begin(), ranges.end(), [&](const TwoRanges &r) {
  // clang-format on
  std::for_each(ranges.begin(), ranges.end(), [&](const TwoRanges &r) {
    const int my_id = [&]() {
      std::lock_guard<std::mutex> guard(int_mutex);
      return ++index;
    }();

    std::vector<int> edges;

    auto action = [&](int i, int j) {
      const auto d = distance(data[i], data[j]);

      if (d < rVal) {
        edges.push_back(i);
        edges.push_back(j);
      }
    };

    auto print_progress = [&](int i) {
      if (i % 50 == 0) {
        std::lock_guard<std::mutex> guard(progress_mutex);
        progress[my_id] = (100.0 * (i - r.i_min)) / (r.i_max - r.i_min);
        if (my_id == print_id) {
          const float p =
              std::accumulate(progress.begin(), progress.end(), 0.f) /
              progress.size();
          std::cout << "Average progress: " << p << "%           \n";
          print_id =
              std::distance(progress.begin(),
                            std::min_element(progress.begin(), progress.end()));
        }
      }
    };

    if (r.i_min == r.j_min && r.i_max == r.j_max) {
      for (int i = r.i_min; i < r.i_max; ++i) {
        print_progress(i);
        for (int j = i + 1; j < r.i_max; ++j) {
          action(i, j);
        }
      }
    } else {
      for (int i = r.i_min; i < r.i_max; ++i) {
        print_progress(i);
        for (int j = r.j_min; j < r.j_max; ++j) {
          action(i, j);
        }
      }
    }

    std::lock_guard<std::mutex> guard(res_mutex);
    for (auto &p : edges) {
      all_edges.emplace_back(std::move(p));
    }
  });

  auto edge_index = torch::tensor(all_edges)
                        .clone()
                        .reshape({static_cast<int>(all_edges.size() / 2), 2})
                        .transpose(0, 1);

  for (int i = 0; i < 5; ++i) {
    if ((edge_index[0][i].item<int>() != all_edges[2 * i]) ||
        (edge_index[1][i].item<int>() != all_edges[2 * i + 1])) {
      throw std::runtime_error("reshape error");
    }
  }

  return edge_index;
}
