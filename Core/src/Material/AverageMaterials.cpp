// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/detail/AverageMaterials.hpp"

#include "Acts/Material/Material.hpp"

#include <cmath>

Acts::MaterialSlab Acts::detail::combineSlabs(const MaterialSlab& slab1,
                                              const MaterialSlab& slab2) {
  const auto& mat1 = slab1.material();
  const auto& mat2 = slab2.material();

  // NOTE 2020-08-26 msmk
  // the following computations provide purely geometrical averages of the
  // material properties. what we are really interested in are the material
  // properties that best describe the material interaction of the averaged
  // material. It is unclear if the latter can be computed in a general fashion,
  // so we stick to the purely geometric averages for now.

  // NOTE 2021-03-24 corentin
  // In the case of the atomic number, the averaging is based on the log
  // to properly account for the energy loss in multiple material.

  // use double for (intermediate) computations to avoid precision loss

  // the thickness properties are purely additive
  const double thickness = static_cast<double>(slab1.thickness()) +
                     static_cast<double>(slab2.thickness());

  // if the two materials are the same there is no need for additional
  // computation
  if (mat1 == mat2) {
    return {mat1, static_cast<float>(thickness)};
  }

  const double thicknessInX0 = static_cast<double>(slab1.thicknessInX0()) +
                         static_cast<double>(slab2.thicknessInX0());
  const double thicknessInL0 = static_cast<double>(slab1.thicknessInL0()) +
                         static_cast<double>(slab2.thicknessInL0());

  // molar amount-of-substance assuming a unit area, i.e. volume = thickness*1*1
  const double molarAmount1 = static_cast<double>(mat1.molarDensity()) *
                        static_cast<double>(slab1.thickness());
  const double molarAmount2 = static_cast<double>(mat2.molarDensity()) *
                        static_cast<double>(slab2.thickness());
  const double molarAmount = molarAmount1 + molarAmount2;

  // handle vacuum specially
  if (!(0.0 < molarAmount)) {
    return {Material(), static_cast<float>(thickness)};
  }

  // radiation/interaction length follows from consistency argument
  const float x0 = static_cast<float>(thickness / thicknessInX0);
  const float l0 = static_cast<float>(thickness / thicknessInL0);

  // compute average molar density by dividing the total amount-of-substance by
  // the total volume for the same unit area, i.e. volume = totalThickness*1*1
  const float molarDensity = static_cast<float>(molarAmount / thickness);

  // assume two slabs of materials with N1,N2 atoms/molecules each with atomic
  // masses A1,A2 and nuclear charges. We have a total of N = N1 + N2
  // atoms/molecules and should have average atomic masses and nuclear charges
  // of
  //
  //     A = (N1*A1 + N2*A2) / (N1+N2) = (N1/N)*A1 + (N2/N)*A2 = W1*A1 + W2*A2
  //
  // the number of atoms/molecules in a given volume V with molar density rho is
  //
  //     N = V * rho * Na
  //
  // where Na is the Avogadro constant. the weighting factors follow as
  //
  //     Wi = (Vi*rhoi*Na) / (V1*rho1*Na + V2*rho2*Na)
  //        = (Vi*rhoi) / (V1*rho1 + V2*rho2)
  //
  // which can be computed from the molar amount-of-substance above.

  const double molarWeight1 = molarAmount1 / molarAmount;
  const double molarWeight2 = molarAmount2 / molarAmount;
  const float ar = static_cast<float>(molarWeight1 * static_cast<double>(mat1.Ar()) + molarWeight2 * static_cast<double>(mat2.Ar()));

  // In the case of the atomic number, its main use is the computation
  // of the mean excitation energy approximated in ATL-SOFT-PUB-2008-003 as :
  //     I = 16 eV * Z^(0.9)
  // This mean excitation energy will then be used to compute energy loss
  // which will be proportional to :
  //     Eloss ~ ln(1/I)*thickness
  // In the case of two successive material :
  //     Eloss = El1 + El2
  //           ~ ln(Z1)*t1 + ln(Z2)*t2
  //           ~ ln(Z)*t
  // To respect this the average atomic number thus need to be defined as :
  //     ln(Z)*t = ln(Z1)*t1 + ln(Z2)*t2
  //           Z = Exp( ln(Z1)*t1/t + ln(Z2)*t2/t )

  const double thicknessWeight1 = static_cast<double>(slab1.thickness()) / thickness;
  const double thicknessWeight2 = static_cast<double>(slab2.thickness()) / thickness;
  float z = 0.f;
  if (mat1.Z() != 0.f && mat2.Z() != 0.f) {
    z = std::exp(thicknessWeight1 * std::log(static_cast<double>(mat1.Z())) +
            thicknessWeight2 * std::log(static_cast<double>(mat2.Z())));
  } else {
    z = static_cast<float>(thicknessWeight1 * static_cast<double>(mat1.Z()) + thicknessWeight2 * static_cast<double>(mat2.Z()));
  }

  return {Material::fromMolarDensity(x0, l0, ar, z, molarDensity),
          static_cast<float>(thickness)};
}
