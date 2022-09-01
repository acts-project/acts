#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace Acts {

template <typename vertex_t, typename edge_t, typename weight_t>
void weaklyConnectedComponents(vertex_t numNodes,
                               std::vector<vertex_t>& rowIndices,
                               std::vector<vertex_t>& colIndices,
                               std::vector<weight_t>& edgeWeights,
                               std::vector<vertex_t>& trackLabels,
                               float edge_cut) {
  typedef boost::adjacency_list<boost::vecS,         // edge list
                                boost::vecS,         // vertex list
                                boost::undirectedS,  // directedness
                                boost::no_property,  // property of vertices
                                float                // property of edges
                                >
      Graph;

  Graph g(numNodes);
  for (size_t idx = 0; idx < rowIndices.size(); ++idx) {
    if (edgeWeights[idx] > edge_cut) {
      boost::add_edge(rowIndices[idx], colIndices[idx], edgeWeights[idx], g);
    }
  }

  [[maybe_unused]] size_t num_components =
      boost::connected_components(g, &trackLabels[0]);
}

}  // namespace Acts
 
