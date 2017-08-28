#include "ACTS/Utilities/BFieldMapUtils.hpp"
#include "ACTS/Utilities/detail/Axis.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"

Acts::InterpolatedBFieldMap::FieldMapper<2, 2> Acts::fieldMapperRZ(
    std::function<size_t(std::array<size_t, 2> binsRZ,
                         std::array<size_t, 2> nBinsRZ)> localToGlobalBin,
    std::vector<double>         rPos,
    std::vector<double>         zPos,
    std::vector<Acts::Vector2D> bField,
    double                      lengthUnit,
    double                      BFieldUnit,
    bool                        firstQuadrant)
{
  // [1] Create Grid
  // sort the values
  std::sort(rPos.begin(), rPos.end());
  std::sort(zPos.begin(), zPos.end());
  // Get unique values
  rPos.erase(std::unique(rPos.begin(), rPos.end()), rPos.end());
  zPos.erase(std::unique(zPos.begin(), zPos.end()), zPos.end());
  rPos.shrink_to_fit();
  zPos.shrink_to_fit();
  // get the number of bins
  size_t nBinsR = rPos.size();
  size_t nBinsZ = zPos.size();
  if (firstQuadrant) {
    nBinsZ = 2. * nBinsZ - 1;
  }
  // get the minimum and maximum
  auto   minMaxR = std::minmax_element(rPos.begin(), rPos.end());
  auto   minMaxZ = std::minmax_element(zPos.begin(), zPos.end());
  double rMin    = *minMaxR.first;
  double zMin    = *minMaxZ.first;
  if (firstQuadrant) {
    zMin = -(*minMaxZ.second);
  }
  // Create the axis for the grid
  Acts::detail::EquidistantAxis rAxis(
      rMin * lengthUnit, (*minMaxR.second) * lengthUnit, nBinsR);
  Acts::detail::EquidistantAxis zAxis(
      zMin * lengthUnit, (*minMaxZ.second) * lengthUnit, nBinsZ);
  // Create the grid
  typedef Acts::detail::Grid<Acts::Vector2D,
                             Acts::detail::EquidistantAxis,
                             Acts::detail::EquidistantAxis>
         Grid_t;
  Grid_t grid(std::make_tuple(std::move(rAxis), std::move(zAxis)));
  // [2] Set the bField values
  size_t nBinsPhi = 1000;

  for (size_t i = 0; i < nBinsR; ++i) {
    for (size_t j = 0; j < nBinsZ; ++j) {
      std::array<size_t, 2> nIndices = {{rPos.size(), zPos.size()}};
      Grid_t::index_t indices = {i, j};
      if (firstQuadrant) {
        size_t          n = std::abs(int(j) - (int(zPos.size()) - 1));
        Grid_t::index_t indicesFirstQuadrant = {i, n};
        grid.at(indices)
            = bField.at(localToGlobalBin(indicesFirstQuadrant, nIndices))
            * BFieldUnit;
      } else
        grid.at(indices)
            = bField.at(localToGlobalBin(indices, nIndices)) * BFieldUnit;
    }
  }

  // [3] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Acts::Vector3D& pos) {
    return Acts::Vector2D(pos.perp(), pos.z());
  };

  // [4] Create the transformation for the bfield
  // map (Br,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Acts::Vector2D& field,
                            const Acts::Vector3D& pos) {
    return Acts::Vector3D(
        field.x() * cos(pos.phi()), field.x() * sin(pos.phi()), field.y());
  };

  // [5] Create the mapper & BField Service
  // create field mapping
  return Acts::InterpolatedBFieldMap::FieldMapper<2, 2>(
      transformPos, transformBField, std::move(grid));
}

Acts::InterpolatedBFieldMap::FieldMapper<3, 3> Acts::fieldMapperXYZ(
    std::function<size_t(std::array<size_t, 3> binsXYZ,
                         std::array<size_t, 3> nBinsXYZ)> localToGlobalBin,
    std::vector<double>         xPos,
    std::vector<double>         yPos,
    std::vector<double>         zPos,
    std::vector<Acts::Vector3D> bField,
    double                      lengthUnit,
    double                      BFieldUnit,
    bool                        firstOctant)
{
  // [1] Create Grid
  // Sort the values
  std::sort(xPos.begin(), xPos.end());
  std::sort(yPos.begin(), yPos.end());
  std::sort(zPos.begin(), zPos.end());
  // Get unique values
  xPos.erase(std::unique(xPos.begin(), xPos.end()), xPos.end());
  yPos.erase(std::unique(yPos.begin(), yPos.end()), yPos.end());
  zPos.erase(std::unique(zPos.begin(), zPos.end()), zPos.end());
  xPos.shrink_to_fit();
  yPos.shrink_to_fit();
  zPos.shrink_to_fit();
  // get the number of bins
  size_t nBinsX = xPos.size();
  size_t nBinsY = yPos.size();
  size_t nBinsZ = zPos.size();
  if (firstOctant) {
    nBinsX = 2 * nBinsX - 1;
    nBinsY = 2 * nBinsY - 1;
    nBinsZ = 2 * nBinsZ - 1;
  }
  // get the minimum and maximum
  auto minMaxX = std::minmax_element(xPos.begin(), xPos.end());
  auto minMaxY = std::minmax_element(yPos.begin(), yPos.end());
  auto minMaxZ = std::minmax_element(zPos.begin(), zPos.end());
  // Create the axis for the grid
  double xMin = *minMaxX.first;
  double yMin = *minMaxY.first;
  double zMin = *minMaxZ.first;
  // If only the first octant is given
  if (firstOctant) {
    xMin = -*minMaxX.second;
    yMin = -*minMaxY.second;
    zMin = -*minMaxZ.second;
  }
  Acts::detail::EquidistantAxis xAxis(
      xMin * lengthUnit, *minMaxX.second * lengthUnit, nBinsX);
  Acts::detail::EquidistantAxis yAxis(
      yMin * lengthUnit, *minMaxY.second * lengthUnit, nBinsY);
  Acts::detail::EquidistantAxis zAxis(
      zMin * lengthUnit, *minMaxZ.second * lengthUnit, nBinsZ);
  // Create the grid
  typedef Acts::detail::Grid<Acts::Vector3D,
                             Acts::detail::EquidistantAxis,
                             Acts::detail::EquidistantAxis,
                             Acts::detail::EquidistantAxis>
         Grid_t;
  Grid_t grid(
      std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)));

  // [2] Set the bField values
  for (size_t i = 0; i < nBinsX; ++i) {
    for (size_t j = 0; j < nBinsY; ++j) {
      for (size_t k = 0; k < nBinsZ; ++k) {
        Grid_t::index_t indices = {i, j, k};
        std::array<size_t, 3> nIndices
            = {{xPos.size(), yPos.size(), zPos.size()}};
        if (firstOctant) {
          size_t          m = std::abs(int(i) - (int(xPos.size()) - 1));
          size_t          n = std::abs(int(j) - (int(yPos.size()) - 1));
          size_t          l = std::abs(int(k) - (int(zPos.size()) - 1));
          Grid_t::index_t indicesFirstOctant = {m, n, l};

          grid.at(indices)
              = bField.at(localToGlobalBin(indicesFirstOctant, nIndices))
              * BFieldUnit;

        } else
          grid.at(indices)
              = bField.at(localToGlobalBin(indices, nIndices)) * BFieldUnit;
      }
    }
  };
  // [3] Create the transformation for the position
  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Acts::Vector3D& pos) { return pos; };
  // [4] Create the transformation for the bfield
  // map (Bx,By,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Acts::Vector3D& field,
                            const Acts::Vector3D& pos) { return field; };

  // [5] Create the mapper & BField Service
  // create field mapping
  return Acts::InterpolatedBFieldMap::FieldMapper<3, 3>(
      transformPos, transformBField, std::move(grid));
}
