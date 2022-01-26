#include <cugraph/algorithms.hpp>
#include <cugraph/graph.hpp>
#include <cugraph/graph_functions.hpp>
#include <cugraph/graph_view.hpp>

#include "utilities/test_graphs.hpp"
#include "utilities/high_res_clock.h"

#include <iostream>
#include <string>
#include <utility>

int main(int argc, char* argv[])
{
  std::string filename("polbooks.mtx");
  int opt;
  while ((opt = getopt(argc, argv, "f:")) != -1){
    switch (opt)
    {
    case 'f':
      filename = optarg;
      break;
    
    default:
      break;
    }
  }
  cugraph::test::File_Usecase const & input_usecase = cugraph::test::File_Usecase(filename);

  constexpr bool renumber = true;

  using vertex_t = int32_t;
  using edge_t = int32_t;
  using weight_t = float;

  static int PERF = 0;

  raft::handle_t handle{};
  HighResClock hr_clock{};

  if (PERF) {
    CUDA_TRY(cudaDeviceSynchronize());  // for consistent performance measurement
    hr_clock.start();
  }

  // the last two booleans are:
  // store_transposed and multi-gpu
  // cugraph::graph_t<vertex_t, edge_t, weight_t, false, false> graph(handle);
  // rmm::device_uvector<vertex_t> d_renumber_map_labels(0, handle.get_stream());

  auto [graph, d_renumber_map_labes] = 
    cugraph::test::construct_graph<vertex_t, edge_t, weight_t, false, false>(
      handle, input_usecase, false, renumber);

  auto graph_view = graph.view();
  CUDA_TRY(cudaDeviceSynchronize());  // for consistent performance measurement

  std::cout << "Number of Nodes:" << graph_view.get_number_of_vertices() << std::endl;
  std::cout << "Number of Edges:" << graph_view.get_number_of_edges() << std::endl;

  rmm::device_uvector<vertex_t> d_components(graph_view.get_number_of_vertices(),
                                              handle.get_stream());

  cugraph::weakly_connected_components(
        handle, graph_view, d_components.data());

  std::vector<vertex_t> h_cugraph_components(graph_view.get_number_of_vertices());
  raft::update_host(h_cugraph_components.data(),
                d_components.data(),
                d_components.size(),
                handle.get_stream());
  int idx = 0;
  int max_sp = 100;
  std::cout <<"print " << max_sp << " spacepoints" << std::endl;
  std::cout << "size of components: " << h_cugraph_components.size() << std::endl;
  for(auto trkid: h_cugraph_components){
    std::cout << idx << " " << trkid << std::endl;
    idx += 1;
    if (idx > max_sp) break;
  }
  // std::unordered_map<vertex_t, vertex_t> cuda_to_reference_map{};
  // for (size_t i = 0; i < h_reference_components.size(); ++i) {
  //   cuda_to_reference_map.insert({h_cugraph_components[i], h_reference_components[i]});
  // }


  return 0;
}
