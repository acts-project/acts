#pragma once
#include <ATen/ATen.h>
#include <torch/extension.h>


void PrefixSumCUDA(
    const at::Tensor grid_cnt,
    int num_grids,
    at::Tensor grid_off);


void PrefixSumCPU(
    const at::Tensor grid_cnt,
    int num_grids,
    at::Tensor grid_off);

/***
PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.def("prefix_sum_cuda", &PrefixSumCUDA);
  m.def("prefix_sum_cpu", &PrefixSumCPU);
}
***/