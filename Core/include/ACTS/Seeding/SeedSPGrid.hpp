#ifndef SEEDSPGRIDCREATOR_H
#define SEEDSPGRIDCREATOR_H

#include "ACTS/Utilities/detail/Grid.hpp"
#include "ACTS/Seeding/SPForSeed.hpp"
#include "ACTS/Seeding/SeedmakerConfig.hpp"
#include <memory>


namespace Acts{
namespace Seeding{

using SPGrid = detail::Grid<std::vector<std::shared_ptr<SPForSeed> >,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Closed>,
detail::Axis<detail::AxisType::Equidistant, detail::AxisBoundaryType::Bound>>;

class SPGridCreator{
public:

static std::unique_ptr<SPGrid> createGrid(std::unique_ptr<Acts::Seeding::Config>& cfg);

};
}
}
#endif //SEEDSPGRIDCREATOR_H
