#include <utilities/test_graphs.hpp>
#include <utilities/high_res_clock.h>
#include <utilities/base_fixture.hpp>
#include <utilities/test_utilities.hpp>

#include <cugraph/partition_manager.hpp>
#include <cugraph/utilities/error.hpp>

#include <cugraph/graph.hpp>
#include <cugraph/algorithms.hpp>
#include <cugraph/graph.hpp>
#include <cugraph/graph_functions.hpp>
#include <cugraph/graph_view.hpp>

#include <raft/cudart_utils.h>
#include <raft/handle.hpp>

// #include <thrust/sequence.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cerrno>
#include <cstring>
#include <tuple>
#include <boost/range/combine.hpp>

#ifndef CUDA_RT_CALL
#define CUDA_RT_CALL(call)                                                    \
  {                                                                           \
    cudaError_t cudaStatus = call;                                            \
    if (cudaSuccess != cudaStatus)                                            \
      fprintf(stderr,                                                         \
              "ERROR: CUDA RT call \"%s\" in line %d of file %s failed with " \
              "%s (%d).\n",                                                   \
              #call,                                                          \
              __LINE__,                                                       \
              __FILE__,                                                       \
              cudaGetErrorString(cudaStatus),                                 \
              cudaStatus);                                                    \
  }
#endif  // CUDA_RT_CALL


template <typename vertex_t, typename edge_t, typename weight_t>
__global__ void weakly_connected_components(
  std::vector<vertex_t>& rowIndices, std::vector<vertex_t>& colIndices,
  std::vector<weight_t>& edgeWeights,
  std::vector<vertex_t>& trackLabels
  )
{
    std::cout << "Weakly components Start" << std::endl;
    std::cout << "edge size: " << rowIndices.size() << " " << colIndices.size() << std::endl;
    raft::handle_t handle{};

    cugraph::graph_t<vertex_t, edge_t, weight_t, false, false> graph(handle);

    
    cudaStream_t stream;
    CUDA_RT_CALL(cudaStreamCreate(&stream));
    handle.set_stream(stream);

    //std::vector<std::tuple<size_t, size_t>>({{100, 0}})
    constexpr bool renumber = true;
    using store_transposed = bool;

    static int PERF = 0;
    HighResClock hr_clock{};

    if (PERF) {
      CUDA_TRY(cudaDeviceSynchronize());  // for consistent performance measurement
      hr_clock.start();
    }

    // learn from matrix_market_file_utilities.cu
    vertex_t maxVertexID_row = *std::max_element(rowIndices.begin(), rowIndices.end());
    vertex_t maxVertexID_col = *std::max_element(colIndices.begin(), colIndices.end());
    vertex_t maxVertex = std::max(maxVertexID_row, maxVertexID_col);

    
    vertex_t number_of_vertices = maxVertex;
    rmm::device_uvector<vertex_t> d_vertices(number_of_vertices, handle.get_stream());
    std::vector<vertex_t> vertex_idx(number_of_vertices);
    for(vertex_t idx=0; idx < number_of_vertices; idx++){
      vertex_idx[idx] = idx;
    }

    std::cout << "vertices size: " << vertex_idx.size() << std::endl;

    rmm::device_uvector<vertex_t> src_v(rowIndices.size(), handle.get_stream());
    rmm::device_uvector<vertex_t> dst_v(colIndices.size(), handle.get_stream());
    rmm::device_uvector<weight_t> weights_v(edgeWeights.size(), handle.get_stream());

    raft::update_device(src_v.data(), rowIndices.data(), rowIndices.size(), handle.get_stream());
    raft::update_device(dst_v.data(), colIndices.data(), colIndices.size(), handle.get_stream());
    raft::update_device(weights_v.data(), edgeWeights.data(), edgeWeights.size(), handle.get_stream());
    raft::update_device(d_vertices.data(), vertex_idx.data(),
                            vertex_idx.size(), handle.get_stream());
    //cugraph::graph_t<vertex_t, edge_t, weight_t, false, false> graph(handle);

    //std::tuple<cugraph::graph_t<vertex_t, edge_t, weight_t,store_transposed, multi_gpu>,rmm::device_uvector<vertex_t>> mygraph;
    std::cout << "create graph from edgelist" << std::endl;
    std::tie(graph, std::ignore) = cugraph::create_graph_from_edgelist<
      vertex_t, edge_t, weight_t, false, false>(
      handle,
      std::move(d_vertices),
      std::move(src_v),
      std::move(dst_v),
      std::move(weights_v),
      cugraph::graph_properties_t{true, false},
      false);

    // the last two booleans are:
    // store_transposed and multi-gpu
    //cugraph::graph_t<vertex_t, edge_t, weight_t, false, false> graph(handle);
    //rmm::device_uvector<vertex_t> d_renumber_map_labels(0, handle.get_stream());
    std::cout<<"After1: "<<std::endl;
    //std::tie(graph, d_renumber_map_labels) =
    //  input_usecase.template construct_graph<vertex_t, edge_t, weight_t, false, false>(
    //    handle, false, renumber);

    std::cout << "back from construct_graph" << std::endl;
    auto graph_view = graph.view();
    CUDA_TRY(cudaDeviceSynchronize());  // for consistent performance measurement
    
    
    int num_edges = graph_view.get_number_of_edges();
    int num_vert = graph_view.get_number_of_vertices();
    std::cout << "Number of Nodes:" << num_vert << std::endl;
    std::cout << "Number of Edges:" << num_edges << std::endl;
    //ASSERT_TRUE(graph_view.is_symmetric())
    //  << "Weakly connected components works only on undirected (symmetric) graphs.";
    std::cout << "Graph sym: "<< graph_view.is_symmetric() << std::endl;
    rmm::device_uvector<vertex_t> d_components(graph_view.get_number_of_vertices(),
                                                      handle.get_stream());
 
    std::cout << "2back from construct_graph" << std::endl;
    cugraph::weakly_connected_components(handle, graph_view, d_components.data());

    std::cout << "number of components: " << d_components.size() << std::endl;
    
    raft::update_host(trackLabels.data(),
      d_components.data(), d_components.size(), handle.get_stream());
}