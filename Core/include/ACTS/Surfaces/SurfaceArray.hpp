// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SURFACES_SURFACEARRAY_H
#define ACTS_SURFACES_SURFACEARRAY_H

#include <iostream>
#include <type_traits>
#include <vector>
#include "ACTS/Surfaces/concept/AnySurfaceGridLookup.hpp"
#include "ACTS/Utilities/concept/AnyGrid.hpp"
#include "ACTS/Utilities/detail/Axis.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"
//#include <boost/any.hpp>
#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/any_cast.hpp>
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinningType.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

namespace bte = boost::type_erasure;

template <typename T, typename U>
constexpr bool
eqTypes()
{
  return std::is_same<T, U>::value;
}

using SurfaceVector = std::vector<const Surface*>;
template <class... Axes>
using SurfaceGrid = detail::Grid<SurfaceVector, Axes...>;

class SurfaceArray
{
  using StoredType = SurfaceVector;
  template <class Point, size_t DIM>
  using AnyGrid_t              = concept::AnyNDimGrid<StoredType, Point, DIM>;
  using AnySurfaceGridLookup_t = concept::AnySurfaceGridLookup<StoredType>;

public:
  template <size_t DIM>
  struct SurfaceLookup
  {
  public:
    SurfaceLookup(std::function<ActsVectorD<DIM>(const Vector3D&)> transformPos,
                  AnyGrid_t<ActsVectorD<DIM>, DIM> grid)
      : m_transformPos(std::move(transformPos)), m_grid(std::move(grid))
    {
    }

    StoredType&
    lookup(const Vector3D& pos)
    {
      // size_t bi = m_grid.getGlobalBinIndex(m_transformPos(pos));
      // std::cout << " -> bin " << bi << std::endl;
      return m_grid.at(m_transformPos(pos));
    }

    const StoredType&
    lookup(const Vector3D& pos) const
    {
      return m_grid.at(m_transformPos(pos));
    }

    std::vector<concept::AnyAxis<>>
    getAxes() const
    {
      auto arr = m_grid.getAxes();
      return std::vector<concept::AnyAxis<>>(arr.begin(), arr.end());
    }

    static constexpr size_t
    dimensions()
    {
      return DIM;
    }

  private:
    std::function<ActsVectorD<DIM>(const Vector3D&)> m_transformPos;
    AnyGrid_t<ActsVectorD<DIM>, DIM> m_grid;
  };

  SurfaceArray(AnySurfaceGridLookup_t gridLookup, SurfaceVector surfaces)
    : m_gridLookup(std::move(gridLookup)), m_surfaces(surfaces)
  {

    std::cout << "SurfaceArray go" << std::endl;

    // let's fill the grid
    for (const auto& srf : m_surfaces) {
      Vector3D pos = srf->binningPosition(binR);
      // std::cout << "fill" << std::endl;
      m_gridLookup.lookup(pos).push_back(srf);
    }
  }

  StoredType&
  at(const Vector3D& pos)
  {
    return m_gridLookup.lookup(pos);
  }

  const StoredType&
  at(const Vector3D& pos) const
  {
    return m_gridLookup.lookup(pos);
  }

  std::ostream&
  dump(std::ostream& sl) const
  {
    sl << "SurfaceArray:" << std::endl;
    sl << " - no surfaces: " << m_surfaces.size() << std::endl;
    sl << " - grid dim:    " << m_gridLookup.dimensions() << std::endl;

    auto axes = m_gridLookup.getAxes();

    for (size_t j = 0; j < axes.size(); ++j) {
      // auto boundaries = axes.boundaries();
      // std::cout << "axis is " << axis.isEquidistant() << std::endl;
      std::cout << " - axis " << (j + 1) << std::endl;
      std::cout << "   - bin edges: [ ";
      auto binEdges = axes.at(j).getBinEdges();
      for (size_t i = 0; i < binEdges.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << binEdges.at(i);
      }
      std::cout << " ]" << std::endl;
    }
    return sl;
  }

private:
  AnySurfaceGridLookup_t m_gridLookup;
  StoredType             m_surfaces;

  // void *p_grid;
  // GridType m_gridType;

  // template <size_t DIM>
  // void
  //_completeBinning()
  //{

  //// Copy and cast, this is relatively expensive, but we only do
  //// it once. Be sure to write to the original grid instead of the copyy!
  // AnyNDimGridT<DIM> nDimGrid = m_grid; // <-- this creates a copy!

  ////typename AnyNDimGridT<DIM>::index_t idx;
  ////typename AnyNDimGridT<DIM>::point_t binCtr;
  // std::array<size_t, DIM> idx;
  // std::array<double, DIM> binCtr_;
  // Vector3D binCtr;
  // double minDist;
  // double curDist;
  // const Surface* minSrf;
  // for (size_t i=0;i<m_grid.size();++i) {

  //// check if empty first, if so, skip
  // if (nDimGrid.at(i).size() > 0) continue;

  // idx = nDimGrid.getLocalBinIndices(i);
  // binCtr_ = nDimGrid.getBinCenter(idx);
  // binCtr << binCtr_[0], binCtr_[1], binCtr_[2];
  //// loop over all surfaces to find closest
  // minDist = std::numeric_limits<double>::max();
  // minSrf = nullptr;
  // for(const auto &srf : m_surfaces) {
  //// why does binningPosition take BinningValue?
  // curDist = (binCtr - srf->binningPosition(binR)).mag();
  // if (minDist < curDist) {
  // minSrf = srf;
  // minDist = curDist;
  //}
  //}

  //// add to bin
  // nDimGrid.at(i).push_back(minSrf);
  //}
  //}
};

}  // namespace Acts

#endif  // ACTS_SURFACES_SURFACEARRAY_H
