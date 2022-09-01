#pragma once

#include <torch/script.h>
#include <torch/torch.h>

namespace Acts {

torch::Tensor buildEdges(at::Tensor& embedFeatures, int64_t numSpacepoints,
                         int dim, float rVal, int kVal);

#if 0
torch::Tensor buildEdgesBruteForce(at::Tensor& embedFeatures,
                                   int64_t numSpacepoints, int dim, float rVal,
                                   int kVal);
#endif

}  // namespace Acts
