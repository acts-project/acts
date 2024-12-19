// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <iostream>
#include <random>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <thrust/execution_policy.h>
#include <thrust/scan.h>

namespace Acts::detail {

template <typename T>
__device__ void swap(T &a, T &b) {
  T tmp = a;
  a = b;
  b = tmp;
}

// https://arxiv.org/abs/1910.05971
__global__ void connectedComponents(const int *sourceEdges,
                                    const int *targetEdges, int *labels,
                                    int *labelsNext, int numEdges,
                                    int numNodes) {
  for (int i = threadIdx.x; i < numNodes; i += blockDim.x) {
    labels[i] = i;
    labelsNext[i] = i;
  }

  __shared__ bool changed;

  while (true) {
    changed = false;
    __syncthreads();

    //printf("Iteration %i\n", n);

    // Tree hooking for each edge;
    for (int i = threadIdx.x; i < numEdges; i += blockDim.x) {
      auto u = sourceEdges[i];
      auto v = targetEdges[i];

      if (labels[u] == labels[labels[u]] && labels[v] < labels[u]) {
        labelsNext[labels[u]] = labels[v];
        changed = true;
        //printf("Edge (%i,%i): set labelsNext[%i] = labels[%i] = %i\n", u, v, labels[u], v, labels[v]);
      } else if (labels[v] == labels[labels[v]] && labels[u] < labels[v]) {
        labelsNext[labels[v]] = labels[u];
        changed = true;
        //printf("Edge (%i,%i): set labelsNext[%i] = labels[%i] = %i\n", u, v, labels[v], u, labels[u]);
      } else {
        //printf("Edge (%i,%i): no action\n", u, v);
      }
    }
    __syncthreads();

    for (int i = threadIdx.x; i < numNodes; i += blockDim.x) {
      labels[i] = labelsNext[i];
    }

    /*if(threadIdx.x == 0 ) {
      for(int i=0; i<numNodes; ++i) {
        printf("Vertex %i - label %i\n", i, labels[i]);
      }
    }*/

    // Shortcutting
    for (int i = threadIdx.x; i < numNodes; i += blockDim.x) {
      if (labels[i] != labels[labels[i]]) {
        labelsNext[i] = labels[labels[i]];
        //printf("Vertex %i: labelsNext[%i] = labels[%i] = %i\n", i, i, labels[i], labels[labels[i]]);
        changed = true;
      }
    }

    for (int i = threadIdx.x; i < numNodes; i += blockDim.x) {
      labels[i] = labelsNext[i];
    }

    /*if(threadIdx.x == 0 ) {
      for(int i=0; i<numNodes; ++i) {
        printf("Vertex after Shortcutting %i - label %i\n", i, labels[i]);
      }
    }*/

    __syncthreads();

    if (!changed) {
      //printf("break!\n");
      break;
    }
  }
}

__global__ void makeLabelMask(const int *labels, int *labelMask) {
  int i = threadIdx.x + blockDim.x * blockIdx.x;

  labelMask[labels[i]] = 1;
}

}  // namespace Acts::detail
