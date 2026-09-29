// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Gnn/detail/buildEdges.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/KDTree.hpp"
#include "ActsPlugins/Gnn/detail/TensorVectorConversion.hpp"

#include <iostream>
#include <mutex>
#include <vector>

#include <torch/script.h>
#include <torch/torch.h>

using namespace torch::indexing;

using namespace Acts;

torch::Tensor ActsPlugins::detail::postprocessEdgeTensor(torch::Tensor edges,
                                                         bool removeSelfLoops,
                                                         bool removeDuplicates,
                                                         bool flipDirections) {
  // Remove self-loops
  if (removeSelfLoops) {
    torch::Tensor selfLoopMask = edges.index({0}) != edges.index({1});
    edges = edges.index({Slice(), selfLoopMask});
  }

  // Remove duplicates
  if (removeDuplicates) {
    torch::Tensor mask = edges.index({0}) > edges.index({1});
    edges.index_put_({Slice(), mask}, edges.index({Slice(), mask}).flip(0));
    edges = std::get<0>(torch::unique_dim(edges, -1, false));
  }

  // Randomly flip direction
  if (flipDirections) {
    torch::Tensor random_cut_keep = torch::randint(2, {edges.size(1)});
    torch::Tensor random_cut_flip = 1 - random_cut_keep;
    torch::Tensor keep_edges =
        edges.index({Slice(), random_cut_keep.to(torch::kBool)});
    torch::Tensor flip_edges =
        edges.index({Slice(), random_cut_flip.to(torch::kBool)}).flip({0});
    edges = torch::cat({keep_edges, flip_edges}, 1);
  }

  return edges.toType(torch::kInt64);
}

/// This is a very unsophisticated span implementation to avoid data copies in
/// the KDTree search.
/// Should be replaced with std::span when possible
template <typename T, std::size_t S>
struct Span {
  T *ptr;

  auto size() const { return S; }

  using const_iterator = T const *;
  const_iterator cbegin() const { return ptr; }
  const_iterator cend() const { return ptr + S; }

  auto &operator[](std::size_t i) const { return ptr[i]; }
};

template <std::size_t Dim>
float dist(const Span<float, Dim> &a, const Span<float, Dim> &b) {
  float s = 0.f;
  for (auto i = 0ul; i < Dim; ++i) {
    s += (a[i] - b[i]) * (a[i] - b[i]);
  }
  return std::sqrt(s);
};

template <std::size_t Dim>
struct BuildEdgesKDTree {
  static torch::Tensor invoke(torch::Tensor &embedFeatures, float rVal,
                              int kVal) {
    assert(embedFeatures.size(1) == Dim);
    embedFeatures = embedFeatures.to(torch::kCPU);

    ////////////////
    // Build tree //
    ////////////////
    using KDTree = KDTree<Dim, int, float, Span>;

    typename KDTree::vector_t features;
    features.reserve(embedFeatures.size(0));

    auto dataPtr = embedFeatures.data_ptr<float>();

    for (int i = 0; i < embedFeatures.size(0); ++i) {
      features.push_back({Span<float, Dim>{dataPtr + i * Dim}, i});
    }

    KDTree tree(std::move(features));

    /////////////////
    // Search tree //
    /////////////////
    std::vector<std::int32_t> edges;
    edges.reserve(2 * kVal * embedFeatures.size(0));

    for (int iself = 0; iself < embedFeatures.size(0); ++iself) {
      const Span<float, Dim> self{dataPtr + iself * Dim};

      RangeXD<Dim, float> range;
      for (auto j = 0ul; j < Dim; ++j) {
        range[j] = Range1D<float>(self[j] - rVal, self[j] + rVal);
      }

      tree.rangeSearchMapDiscard(
          range, [&](const Span<float, Dim> &other, const int &iother) {
            if (iself != iother && dist(self, other) <= rVal) {
              edges.push_back(iself);
              edges.push_back(iother);
            }
          });
    }

    // Transpose is necessary here, clone to get ownership
    return ActsPlugins::detail::vectorToTensor2D(edges, 2).t().clone();
  }
};

torch::Tensor ActsPlugins::detail::buildEdgesKDTree(
    torch::Tensor &embedFeatures, float rVal, int kVal, bool flipDirections) {
  auto tensor = template_switch<BuildEdgesKDTree, 1, 12>(
      embedFeatures.size(1), embedFeatures, rVal, kVal);

  return postprocessEdgeTensor(tensor, true, true, flipDirections);
}

torch::Tensor ActsPlugins::detail::buildEdges(torch::Tensor &embedFeatures,
                                              float rVal, int kVal,
                                              bool flipDirections) {
  if (embedFeatures.is_cuda()) {
    return detail::buildEdgesFRNN(embedFeatures, rVal, kVal, flipDirections);
  }
  return detail::buildEdgesKDTree(embedFeatures, rVal, kVal, flipDirections);
}
