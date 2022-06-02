#pragma once 

#include <torch/torch.h>
#include <torch/script.h>
using namespace torch::indexing;

#include <grid.h>
#include <insert_points.h>
#include <counting_sort.h>
#include <prefix_sum.h>
#include <find_nbrs.h>
#include "cuda.h"
#include "cuda_runtime_api.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

torch::Tensor buildEdges(
    at::Tensor& embedFeatures, int64_t numSpacepoints,
    int dim, float rVal, int kVal
);

void buildEdges(
    std::vector<float>& embedFeatures,
    std::vector<int64_t>& edgeList,
    int64_t numSpacepoints,
    int embeddingDim,    // dimension of embedding space
    float rVal, // radius of the ball
    int kVal    // number of nearest neighbors
);

template <typename vertex_t, typename edge_t, typename weight_t>
void weaklyConnectedComponents(
vertex_t numNodes,
std::vector<vertex_t>& rowIndices,
std::vector<vertex_t>& colIndices,
std::vector<weight_t>& edgeWeights,
std::vector<vertex_t>& trackLabels) 
{
typedef
    boost::adjacency_list<
    boost::vecS            // edge list
, boost::vecS            // vertex list
, boost::undirectedS     // directedness
, boost::no_property     // property associated with vertices
, float                  // property associated with edges
> Graph; 

Graph g(numNodes);
float edge_cut = 0.75;
for(size_t idx=0; idx < rowIndices.size(); ++idx) {
    if (edgeWeights[idx] > edge_cut) {
        boost::add_edge(rowIndices[idx], colIndices[idx], edgeWeights[idx], g);
    }
}
size_t num_components = boost::connected_components(g, &trackLabels[0]);

}