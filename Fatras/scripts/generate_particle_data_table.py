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


def generate_code():
    """
    Generate
    """
    # ensure the rows are sorted by the signed pdg number (first column)
    # name, c++ type, and output format for each column
    columns = [
        ("three_charge", "Charge", "float", "{}", convert_charge_to_e),
        ("mass", "Mass", "float", "{}f", convert_mass_to_GeV),
        ("name", "Name", " const  char* const", '"{}"', identity),
    ]
    lines = []
    lines = [
        CODE_HEADER,
        f"static constexpr uint32_t kParticlesCount = {len(Particle.all())}u;",
    ]
    # build a separate array for each column
    for variable_name, target_name, type_name, value_format, transform in columns:

        cpp_name = f"kParticlesMap{target_name}"

        lines.append(
            f"static const std::map<std::int32_t, {type_name}> {cpp_name} = {{"
        )

        for p in Particle.all():
            value = getattr(p, variable_name)

            if variable_name != "Name":
                lines.append(f"// {p.name}")

            if value is None and type_name == "float":
                lines.append(
                    f"{{ {int(p.pdgid)} , std::numeric_limits<float>::quiet_NaN() }},"
                )
            else:
                lines.append(
                    f"{{ {int(p.pdgid)}, {value_format.format(transform(value))} }},"
                )

        lines.append("};")

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

    args = p.parse_args()
    if args.output is None:
        output_file = sys.stdout
    else:
        # will overwrite existing file
        output_file = io.open(args.output, mode="wt", encoding="utf-8")

    code = generate_code()
    if args.format:
        code = clang_format(code)
    output_file.write(code)
