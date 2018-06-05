#pragma once

#include "Acts/Utilities/detail/Grid.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Seeding/SPForSeed.hpp"
#include <memory>


namespace Acts{

struct SeedingGridConfig{
    // magnetic field in kTesla
    float bFieldInZ = 0.00208;
    // minimum pT to be found by seedfinder in MeV
    float minPt = 400;
    // maximum extension of sensitive detector layer relevant for seeding as distance from x=y=0 (i.e. in r) in mm
    float rMax = 563;
    // maximum extension of sensitive detector layer relevant for seeding in positive direction in z in mm
    float zMax = 2750;
    // maximum extension of sensitive detector layer relevant for seeding in negative direction in z in mm
    float zMin = -2750;
    // maximum distance in r from middle space point to bottom or top spacepoint in mm
    float deltaRMax = 150;
    // maximum forward direction expressed as cot(theta)
    float cotThetaMax;
};

using SPGrid = detail::Grid<std::vector<std::shared_ptr<SPForSeed> >,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed>,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>>;

class SPGridCreator{
public:

static std::unique_ptr<SPGrid> createGrid(const Acts::SeedingGridConfig& cfg);

};
}
