
#include "Acts/TrackFinding/RoiDescriptor.hpp"

#include <cmath>
#include <sstream>

namespace Acts {


Acts::RoiDescriptor::RoiDescriptor(double eta, double etaMinus, double etaPlus, 
          double phi, double phiMinus, double phiPlus, 
          double zed, double zedMinus, double zedPlus) 
    :
    m_l1Id(0), m_roiId(0), m_roiWord(0)
{}

Acts::RoiDescriptor::~RoiDescriptor() { }

} 
