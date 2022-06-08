#include <cassert>
#include <iostream>

#include <Eigen/Core>

#include <torch/torch.h>
#include <ATen/ATen.h>

#include "Acts/Plugins/ExaTrkX/ExaTrkXUtils.hpp"


float distance(const at::Tensor &a, const at::Tensor &b) {
    assert(a.sizes() == b.sizes());
    assert(a.sizes().size() == 1);

    return std::sqrt(((a-b)*(a-b)).sum().item().to<float>());
}

int cantor_pair(int x, int y) {
    return y + ((x+y)*(x+y+1)) / 2;
}

std::pair<int, int> cantor_pair_inverse(int a) {
    auto f = [](int w) -> int { return (w*(w+1))/2; };
    auto q = [](int w) -> int { return std::floor( (std::sqrt(8*w + 1) - 1) / 2 ); };

    auto y = a - f(q(a));
    auto x = q(a) - y;

    return {x,y};
}

void test_random_graph(int emb_dim, int n_nodes, float r, int knn)
{
    // Create a random point cloud
    auto random_features = at::randn({n_nodes, emb_dim});

    // Generate the truth via brute-force
    Eigen::MatrixXf distance_matrix(n_nodes, n_nodes);

    std::vector<int> edges_ref_cantor;
    std::vector<int> edge_counts(n_nodes, 0);

    for(int i=0; i<n_nodes; ++i) {
        for(int j=i; j<n_nodes; ++j) {
            const auto d = distance(random_features[i], random_features[j]);
            distance_matrix(i,j) = d;
            distance_matrix(j,i) = d;

            if( d < r && i != j ) {
                edges_ref_cantor.push_back(cantor_pair(i,j));
                edge_counts[i]++;
            }
        }
    }

    // Check if knn can find all edges
    auto max_edges = *std::max_element(edge_counts.begin(), edge_counts.end());
    if( max_edges > knn) {
        std::cout << "Warning: edge_count per node can be higher than knn, test will fail\n";
    }

    std::cout << "Max edge count: " << max_edges << "\n";

    // Run the edge building
    auto random_features_cuda = random_features.to(torch::kCUDA);
    auto edges_test = buildEdges(random_features_cuda, n_nodes, emb_dim, r, knn);

    // Map the edges to cantor pairs
    std::vector<int> edges_test_cantor;

    for(int i=0; i<edges_test.size(1); ++i) {
        const auto a = edges_test[0][i].item<int>();
        const auto b = edges_test[1][i].item<int>();
        edges_test_cantor.push_back(
            a < b ? cantor_pair(a, b) : cantor_pair(b, a)
        );
    }

    std::sort(edges_ref_cantor.begin(), edges_ref_cantor.end());
    std::sort(edges_test_cantor.begin(), edges_test_cantor.end());

    // Check
    if( (edges_ref_cantor.size() != edges_test_cantor.size()) or
        not std::equal(edges_test_cantor.begin(), edges_test_cantor.end(), edges_ref_cantor.begin()) ) {

      auto print_in_a_but_not_in_b = [&](const auto a, const auto b) {
        std::vector<int> diff;

        std::set_difference(
            a.begin(), a.end(),
            b.begin(), b.end(),
            std::back_inserter(diff)
        );

        if( diff.empty() ) {
            std::cout << "  none\n";
            return;
        }

        for(auto c : diff) {
            const auto [e0, e1] = cantor_pair_inverse(c);
            std::cout << "  element (" << e0 << "," << e1 << ") with d = " << distance_matrix(e0, e1) << "\n";
        }
      };

      std::cout << "Elements in ref but not in test:\n";
      print_in_a_but_not_in_b(edges_ref_cantor, edges_test_cantor);
      std::cout << "Elements in test but not in ref:\n";
      print_in_a_but_not_in_b(edges_test_cantor, edges_ref_cantor);

      throw std::runtime_error("edges mismatch");
    }

    std::cout << "OK!" << std::endl;
}

int main(int argc, char ** argv) {
    int a = 345;
    int b = 23;
    auto [aa, bb] = cantor_pair_inverse(cantor_pair(a,b));
    assert(a == aa);
    assert(b == bb);

    if( a != aa || b != bb ) { throw; }


    std::vector<std::string> args(argv, argv+argc);

    int emb_dim = 3;
    int n_nodes = 20;
    float r = 1.5;
    int knn = 50;
    int seed = std::rand();

    try
    {
        if( argc > 1 ) {
            emb_dim = std::stoi(args.at(1));
            n_nodes = std::stoi(args.at(2));
            r = std::stof(args.at(3));
            knn = std::stoi(args.at(4));
        }

        if( argc == 6 ) {
            seed = std::stoi(args.at(5));
        }
    }
    catch(...) {
        std::cout << "Usage: " << argv[0] << " <emb_dim> <n_nodes> <r> <knn> (<seed>)\n";
        return 1;
    }

    torch::manual_seed(seed);

    std::cout << "Parameters:\n";
    std::cout << "emb_dim: " << emb_dim << "\n";
    std::cout << "n_nodes: " << n_nodes << "\n";
    std::cout << "r:       " << r << "\n";
    std::cout << "knn:     " << knn << "\n";
    std::cout << "seed:    " << seed << "\n";

    test_random_graph(emb_dim, n_nodes, r, knn);
}
