#include "ACTS/Digitization/Digitization.hpp"

const std::vector<std::vector<Acts::DigitizationCell>>
Acts::createClusters(const std::vector<Acts::DigitizationCell>& cells,
                     bool                                       commonCorner)
{
  // the output
  std::vector<std::vector<Acts::DigitizationCell>> mergedCells;
  // create the graph
  Graph graph;
  //  add the needed amount of vertices to the graph
  for (auto cell : cells) {
    // add vertex to graph for each cell
    add_vertex(graph);
  }

  // loop through cells and set the edges
  for (auto cellA = std::begin(cells); cellA != (std::end(cells) - 1);
       cellA++) {
    for (auto cellB = (cellA + 1); cellB != std::end(cells); cellB++) {
      // the channels
      auto c0A = cellA->channel0;
      auto c1A = cellA->channel1;
      auto c0B = cellB->channel0;
      auto c1B = cellB->channel1;
      // the conditions
      bool isNeighbour0   = (c0A == (c0B - 1) || c0A == (c0B + 1));
      bool isSameChannel0 = (c0A == c0B);
      bool isNeighbour1   = (c1A == (c1B - 1) || c1A == (c1B + 1));
      bool isSameChannel1 = (c1A == c1B);
      // distinguish between between merging cells, when they have a common
      // corner (8cell) or a common edge (4cell)
      if (commonCorner) {
        // 8cell
        if ((isNeighbour0 || isSameChannel0)
            && (isNeighbour1 || isSameChannel1)) {
          auto indexA = std::distance(cells.begin(), cellA);
          add_edge(indexA, indexA + 1, graph);
        }
      } else {
        // 4cell
        if ((isNeighbour0 && isSameChannel1)
            || (isNeighbour1 && isSameChannel0)) {
          auto indexA = std::distance(cells.begin(), cellA);
          add_edge(indexA, indexA + 1, graph);
        }
      }
    }  // cellsB
  }    // cellsA

  std::vector<size_t> component(num_vertices(graph));
  size_t              num = connected_components(graph, &component[0]);
  // copy the component map
  std::vector<size_t> keys = component;
  // sort
  std::sort(keys.begin(), keys.end());
  // get the keys
  auto it = std::unique(keys.begin(), keys.end());
  keys.erase(it, keys.end());

  for (auto& comp : keys) {
    std::vector<Acts::DigitizationCell> compCells;
    for (size_t i = 0; i != component.size(); ++i) {
      if (component[i] == comp) {
        compCells.push_back(cells.at(i));
      }
    }
    mergedCells.push_back(compCells);
  }
  return mergedCells;
}
