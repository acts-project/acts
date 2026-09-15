// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <vector>

namespace ActsPlugins::ActsToMille {

class MillePedeSteering {
 public:
  /// solution method
  enum class Strategy {
    inversion,         // built-in inversion, calculates errors
    diagonalization,   // diagonalisation - slower, provides eigenvector list
    decomposition,     // built-in root-free cholesky - does not calc. errors
    fullMINRES,        // fast-approximate MINRES, full storage
    sparseMINRES,      // fast-approximate MINRES, sparse storage
    fullMINRES_QLP,    // improved MINRES, full storage
    sparseMINRES_QLP,  // improved MINRES, sparse storage
    fullLAPACK,  // LAPACK cholesky factorisation, packed triangular storage -
                 // requires MKL/OpenBLAS
    unpackedLAPACK,  // LAPACK cholesky, with unpacked storage for LAPACK-based
                     // constraint elimination - requires MKL/OpenBLAS
    sparsePARDISO    // Intel PARDISO direct sparse solver (targets very sparse
                     // global matries) - requires oneMKL PARDISO
  };

  /// encodes an equality constraint in Millepede, of the shape
  /// f_1 x_1 + f_2 x_2 + ... + f_n x_n = c.
  /// Will be respected through either lagrange multipliers
  /// or householder transformations.
  struct equalityConstraint {
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
  struct Config {
    Strategy strategy = Strategy::inversion;  /// solution method
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
    int matIter = 1;  // iteration up to which the full matrix is recalculated
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
    std::vector<equalityConstraint> constraints =
        {};                                    // equality constraints to apply
    std::vector<std::string> inputFiles = {};  // input files to use
  };

  /// constructor
  explicit MillePedeSteering(Acts::Logging::Level level = Acts::Logging::INFO)
      : m_logger(Acts::getDefaultLogger("MillePedeSteering", level)) {}

  /// @brief attempt to generate a steering file at the specified location.
  /// Will process the arguments passed in the config object to determine the
  /// content.
  /// @return If successful, a valid path to the file. Else an empty path.
  std::filesystem::path generateSteeringFile(
      const std::filesystem::path& destination, const Config& config) const;

 private:
  std::unique_ptr<const Acts::Logger> m_logger;  /// logger

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  void addInputBlock(std::ofstream& out, const Config& conf) const;
  void addSolutionMethod(std::ofstream& out, const Config& conf) const;
  void addConstraintsBlock(std::ofstream& out, const Config& conf) const;
  void addConstraint(std::ofstream& out,
                     const equalityConstraint& constraint) const;
  void addOthers(std::ofstream& out, const Config& conf) const;
};

}  // namespace ActsPlugins::ActsToMille
