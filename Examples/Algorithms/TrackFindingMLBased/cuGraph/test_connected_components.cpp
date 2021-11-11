#include "mmio_read.h"

int main(int argc, char** argv) {
    
    std::cout << "Begin to read" << std::endl;
    auto file_path = "../../datanmodels/CutEdge.csv";
    std::ifstream f (file_path); 
    if (!f.is_open()) {
        perror (("error while opening file " + std::string(file_path)).c_str());
    }

    using vertex_t = int32_t;
    vertex_t a,b;
    std::vector<vertex_t> h_rows;
    std::vector<vertex_t> h_cols;
    std::vector<float> h_weights;
    std::vector<vertex_t> h_components;

    while (f >> a >> b) {
        h_rows.push_back(a);
        h_cols.push_back(b);
        h_weights.push_back(1.0);
    }    


    weakly_connected_components<int32_t,int32_t,float>(h_rows, h_cols, h_weights, h_components);
    return 0;
}