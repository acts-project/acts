#ifndef DIGITIZATION_DIGITIZATION_HPP
#define DIGITIZATION_DIGITIZATION_HPP

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "ACTS/Digitization/DigitizationCell.hpp"

namespace Acts {
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
    Graph;

const std::vector<std::vector<Acts::DigitizationCell>>
createClusters(const std::vector<Acts::DigitizationCell>& cells,
               bool                                       commonCorner = false);
}

#endif  // DIGITIZATION_DIGITIZATION_HPP
