#!/usr/bin/env python3

# /// script
# dependencies = [
#   "particle==0.24.0",
# ]
# ///

#
# use scikit-hep/particle to generate c++ code for the particle data table.
#

import io
import sys
import subprocess
import argparse
from pathlib import Path

from particle import Particle


CODE_HEADER = """\
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The entries within this file have been automatically created using the
// particle data files from Review of Particle Physics
// by the Berkeley Particle Data Group.

#pragma once

#include <cstdint>
#include <map>
#include <limits>

"""

# Only needed if the user wants to add additional particles at runtime from a CSV file, otherwise it can be ignored
CSV_FILE_READER = """
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

inline float readCsvFloatAt(const std::string& csvPath, // path to the CSV file
                            int lineIndex,              // line number to read, one-based
                            int columnIndex) {          // column number to read, one-based
  if (lineIndex == 0 || columnIndex == 0) {
    throw std::invalid_argument("Line and column numbers must be one-based and greater than 0");
  }

  std::ifstream input(csvPath);
  if (!input.is_open()) {
    throw std::runtime_error("Failed to open CSV file");
  }

  std::string targetLine;
  for (int line = 1; line <= lineIndex; ++line) {
    if (!std::getline(input, targetLine)) {
      throw std::runtime_error("Failed to read line from CSV file");
    }
  }

  std::vector<std::string> columns;
  {
    std::stringstream lineStream(targetLine);
    std::string value;
    while (std::getline(lineStream, value, ',')) {
      columns.push_back(value);
    }
  }
  if (static_cast<std::size_t>(columnIndex) > columns.size()) {
    throw std::invalid_argument("Column number is out of bounds");
  }

  try {
    return std::stof(columns[columnIndex - 1]);
  } catch (...) {
    std::string errorMsg = "Failed to convert column value to float: " + columns[columnIndex - 1] + " at line " + std::to_string(lineIndex) + " column " + std::to_string(columnIndex) + " in file " + csvPath;
    throw std::runtime_error(errorMsg.c_str());
  }
  return 0.f;
}

"""


def convert_mass_to_GeV(mass_MeV):
    """
    Convert the mass from MeV to GeV.
    """
    # return mass_MeV
    return mass_MeV / 1000.0


def convert_charge_to_e(charge):
    """
    Convert the charge from three-charge to electron charge.
    """
    return charge / 3.0


def identity(v):
    """
    Return the value unchanged.
    """
    return v


def generate_code(additional_particles_csv: str | None = None):
    """
    Generate C++ source text containing particle lookup maps.

    By default, the generated code contains all particles provided by
    ``particle.Particle.all()`` and creates one ``std::map`` per property:
    charge, mass, and name, keyed by signed PDG ID.

    If ``additional_particles_csv`` is provided, the function also appends
    extra entries from that CSV to each generated map and emits a helper C++
    CSV reader. Those extra values are read from the CSV file at runtime in the
    generated C++ code.

    Parameters
    ----------
    additional_particles_csv : str | None, optional
        Path to a CSV file containing additional particles.

        Expected row format (no header required):
            ``ID,Name,Mass,Charge``

        Where:
        - ``ID`` is a signed integer PDG ID.
        - ``Name`` is a particle name string.
        - ``Mass`` is in MeV.
        - ``Charge`` is in three-charge units.

        Empty lines and lines starting with ``#`` are ignored.
        Invalid rows are skipped with a warning to stderr.

    Returns
    -------
    str
        Generated C++ code as a single string, ending with a newline.
    """

    lines = CODE_HEADER.splitlines()

    if additional_particles_csv is not None:
        lines.append(CSV_FILE_READER)

    lines.append(
        f"static constexpr std::uint32_t kParticlesCount = {len(Particle.all())}u;"
    )

    # ensure the rows are sorted by the signed pdg number (first column)
    # name, c++ type, and output format for each column
    columns = [
        ("three_charge", "Charge", "float", "{}", convert_charge_to_e),
        ("mass", "Mass", "float", "{}f", convert_mass_to_GeV),
        ("name", "Name", " const  char* const", '"{}"', identity),
    ]
    # build a separate array for each column
    for variable_name, target_name, type_name, value_format, transform in columns:
        cpp_name = f"particlesMap{target_name}"
        map_type = f"std::map<std::int32_t, {type_name}>"

        lines.append(f"inline const {map_type}& {cpp_name}() {{")
        lines.append(f"static const {map_type} map = {{")

        for p in Particle.all():
            value = getattr(p, variable_name)

            lines.append(f"// {p.name}")

            if value is None and type_name == "float":
                lines.append(
                    f"{{ {int(p.pdgid)} , std::numeric_limits<float>::quiet_NaN() }},"
                )
            else:
                lines.append(
                    f"{{ {int(p.pdgid)}, {value_format.format(transform(value))} }},"
                )

        # If additional particles are provided, add them to the table that is dynamically generated at runtime
        if additional_particles_csv is not None:
            with open(additional_particles_csv, "r", encoding="utf-8") as f:
                for line_number, line in enumerate(f, start=1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue  # skip empty lines and comments
                    parts = line.split(",")
                    if len(parts) != 4:
                        print("Warning: Skipping invalid line in additional particles CSV:", line, file=sys.stderr)
                        continue
                    pdgid_str, name, mass_str, charge_str = parts
                    try:
                        pdgid = int(pdgid_str)
                        float(mass_str)
                        float(charge_str)
                    except ValueError:
                        print(
                            "Warning: Skipping line with invalid number format in additional particles CSV: ", line,
                            file=sys.stderr)
                        continue
                    lines.append(f"// {name}")

                    if variable_name == "name":
                        lines.append(
                            f"{{ {pdgid}, {value_format.format(name)} }},"
                        )
                    elif variable_name == "mass":
                        lines.append(
                            f"{{ {pdgid}, readCsvFloatAt(\"{additional_particles_csv}\", {line_number}, 3) }},"
                        )
                    elif variable_name == "three_charge":
                        lines.append(
                            f"{{ {pdgid}, readCsvFloatAt(\"{additional_particles_csv}\", {line_number}, 4) }},"
                        )
                    else:
                        raise ValueError(f"Unsupported variable name: {variable_name}")

        lines.append("};")
        lines.append("return map;")
        lines.append("}")

    # ensure we end with a newline
    lines.append("")
    return "\n".join(lines)


def clang_format(content):
    """
    Format the given content using clang-format and return it.
    """
    args = [
        "clang-format",
        "--assume-filename=ParticleData.hpp",
        "-",
    ]
    process = subprocess.run(
        args,
        input=content,
        capture_output=True,
        check=True,
        encoding="utf-8",
        text=True,
    )
    return process.stdout


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate the particle data table.")
    p.add_argument("output", type=Path, nargs="?", default=None, help="Output file.")
    p.add_argument(
        "--format", action="store_true", help="Run clang-format on the output."
    )
    p.add_argument("--additional-particles", type=Path, help="CSV file with additional particles.", default=None)

    args = p.parse_args()
    if args.output is None:
        output_file = sys.stdout
    else:
        # will overwrite existing file
        output_file = io.open(args.output, mode="wt", encoding="utf-8")

    code = generate_code(additional_particles_csv=args.additional_particles)
    if args.format:
        code = clang_format(code)
    output_file.write(code)
