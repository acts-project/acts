// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <filesystem>
#include <vector>

namespace ActsPlugins {

/// solution method.
/// This covers all the solution
/// strategies supported by MP-2.
enum class MillePedeSolutionStrategy {
  Inversion,         ///> built-in inversion, calculates errors
  Diagonalization,   ///> diagonalisation - slower, provides eigenvector list
  Decomposition,     ///> built-in root-free cholesky - does not calc. errors
  FullMINRES,        ///> fast-approximate MINRES, full storage
  SparseMINRES,      ///> fast-approximate MINRES, sparse storage
  FullMINRES_QLP,    ///> improved MINRES, full storage
  SparseMINRES_QLP,  ///> improved MINRES, sparse storage
  FullLAPACK,  ///> LAPACK cholesky factorisation, packed triangular storage -
               // requires MKL/OpenBLAS
  UnpackedLAPACK,  ///> LAPACK cholesky, with unpacked storage for LAPACK-based
                   // constraint elimination - requires MKL/OpenBLAS
  SparsePARDISO    ///> Intel PARDISO direct sparse solver (targets very sparse
                   // global matries) - requires oneMKL PARDISO
};

/// encodes an equality constraint in Millepede, of the shape
/// f_1 x_1 + f_2 x_2 + ... + f_n x_n = c.
/// Will be applied through either lagrange multipliers
/// or householder transformations (depending on strategy).
struct MillePedeEqualityConstraint {
  std::vector<std::pair<int, double>>
      labelsAndWeights;  /// left hand side: labels and associated weights
  double constraint;     /// right hand side.
};

/// Configuration options.
/// The full set of options is described at
/// https://millepede.pages.desy.de/millepede-ii/option_page.html and
/// https://millepede.pages.desy.de/millepede-ii/changes_page.html. Here, only
/// a short summary is given to the function of each, please refer to the
/// official doc. Not all options are exposed - please use the "extraLines"
/// string flag to add any other desired options directly as config lines.
struct MillePedeSteeringConfig {
  using enum MillePedeSolutionStrategy;
  MillePedeSolutionStrategy strategy = Inversion;  /// solution method
  int minIterations = 3;          /// minimum iterations for solution
  double convergenceLimit = 0.8;  /// convergence limit for iteration
  std::tuple<int, int, int> entriesCut = {
      100, 10, 2};  /// entries cut - before (0) and after (1) chi2 cut, and
                    /// scale factor for cases with partial rejected DoF. See
                    /// detailed doc for more.
  int outlierDownweighting =
      3;  // number of iterations over which outlier downweighting takes
          // places (incremental reduction of )
  double downweightFractionCut =
      0.1;  // fraction of outliers allowed before cases are rejected
  int nOMPthreads = 1;  // number of openMP threads used for calculations
  int nIOthreads = 1;   // number of I/O threads used for reading binaries
  int matIter = 1;      // iteration up to which the full matrix is recalculated
  int printCounts = 2;  // enable printout of counts in result file
  std::pair<double, double> chi2Cut = {
      30.,
      6.};  // 3-sigma chi2 cutoff - first and later iterations of local fits
  bool monitorResiduals = true;  // enable built-in residual monitoring
  bool monitorPulls = true;      // enable internal pull monitoring
  bool skipEmptyCons =
      true;  // skip empty constraints (with no variable ali pars)
  bool countRecords = true;  // set parameter counting to record-level
  std::vector<std::string> extraLines =
      {};  // any other config lines to add to the steering file
  std::vector<MillePedeEqualityConstraint> constraints =
      {};                                    // equality constraints to apply
  std::vector<std::string> inputFiles = {};  // input files to use
};

/// @brief attempt to generate a steering file at the specified location.
/// Will process the arguments passed in the config object to determine the
/// content.
/// @return If successful, a valid path to the file. Else an empty path.
Acts::Result<std::filesystem::path> generateMillePedeSteeringFile(
    const std::filesystem::path& destination,
    const MillePedeSteeringConfig& config, const Acts::Logger& logger);

}  // namespace ActsPlugins
