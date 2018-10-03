#pragma once

#include "Acts/Utilities/detail/Grid.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include <memory>


namespace Acts{

struct SeedingGridConfig{
    // magnetic field in kTesla
    float bFieldInZ;
    // minimum pT to be found by seedfinder in MeV
    float minPt;
    // maximum extension of sensitive detector layer relevant for seeding as distance from x=y=0 (i.e. in r) in mm
    float rMax;
    // maximum extension of sensitive detector layer relevant for seeding in positive direction in z in mm
    float zMax;
    // maximum extension of sensitive detector layer relevant for seeding in negative direction in z in mm
    float zMin;
    // maximum distance in r from middle space point to bottom or top spacepoint in mm
    float deltaRMax;
    // maximum forward direction expressed as cot(theta)
    float cotThetaMax;
};

using SPGrid = detail::Grid<std::vector<std::unique_ptr<const InternalSpacePoint> >,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed>,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>>;

class SPGridCreator{
public:

static std::unique_ptr<SPGrid> createGrid(const Acts::SeedingGridConfig& cfg);

};
}
