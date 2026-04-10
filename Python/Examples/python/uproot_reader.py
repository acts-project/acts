# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Uproot-based reader for simulated particles and hits.

Reads the ROOT files written by RootParticleWriter and RootSimHitWriter
using uproot, without requiring the ROOT plugin. Suitable for use in the
PyPI distribution.
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import uproot

import acts
import acts.examples


class UprootReader(acts.examples.IReader):
    """Reads simulated particles and hits from ROOT files using uproot.

    Reads the file formats produced by RootParticleWriter (one entry per event,
    vector-valued branches) and RootSimHitWriter (one row per hit, scalar
    branches grouped by event_id).

    Either or both file paths may be supplied; omitting one skips that
    collection entirely.  All data is loaded upfront for thread-safe concurrent
    access from the sequencer.
    """

    def __init__(
        self,
        particleFilePath: Optional[Union[Path, str]] = None,
        simHitFilePath: Optional[Union[Path, str]] = None,
        outputParticles: str = "particles",
        outputSimHits: str = "simhits",
        level: acts.logging.Level = acts.logging.INFO,
    ):
        acts.examples.IReader.__init__(self, "UprootReader", level)

        self._outputParticles = outputParticles
        self._outputSimHits = outputSimHits

        self._particleHandle = acts.examples.WriteDataHandle(
            self, acts.examples.SimParticleContainer, "OutputParticles"
        )
        self._hitHandle = acts.examples.WriteDataHandle(
            self, acts.examples.SimHitContainer, "OutputSimHits"
        )

        event_id_sets = []

        self._particle_data = None
        if particleFilePath is not None:
            with uproot.open(str(particleFilePath)) as f:
                tree = f["particles"]
                self._particle_data = tree.arrays(library="np")
                event_ids = self._particle_data["event_id"]
                self._particle_entry_map = {
                    int(eid): i for i, eid in enumerate(event_ids)
                }
                event_id_sets.append(set(event_ids.tolist()))

        self._hit_data = None
        if simHitFilePath is not None:
            with uproot.open(str(simHitFilePath)) as f:
                tree = f["hits"]
                self._hit_data = tree.arrays(library="np")
                hit_event_ids = self._hit_data["event_id"]
                self._hit_event_range_map = self._build_hit_event_map(hit_event_ids)
                self._has_barcode_branch = "barcode" in tree.keys()
                event_id_sets.append(set(self._hit_event_range_map.keys()))

        all_ids = set().union(*event_id_sets)
        self._min_event = int(min(all_ids))
        self._max_event = int(max(all_ids))

    @staticmethod
    def _build_hit_event_map(event_ids: np.ndarray) -> dict:
        """Build a {event_id: (start, end)} map from a sorted event_id array."""
        event_map = {}
        if len(event_ids) == 0:
            return event_map
        current_id = int(event_ids[0])
        start = 0
        for i in range(1, len(event_ids)):
            eid = int(event_ids[i])
            if eid != current_id:
                event_map[current_id] = (start, i)
                current_id = eid
                start = i
        event_map[current_id] = (start, len(event_ids))
        return event_map

    def name(self) -> str:
        return "UprootReader"

    def availableEvents(self):
        return (self._min_event, self._max_event + 1)

    def initialize(self):
        if self._particle_data is not None:
            self._particleHandle.initialize(self._outputParticles)
        if self._hit_data is not None:
            self._hitHandle.initialize(self._outputSimHits)
        return acts.examples.ProcessCode.SUCCESS

    def read(self, context):
        event_number = context.eventNumber

        if self._particle_data is not None:
            particles = acts.examples.SimParticleContainer()
            if event_number in self._particle_entry_map:
                idx = self._particle_entry_map[event_number]
                d = self._particle_data
                n = len(d["particle_type"][idx])
                u = acts.UnitConstants
                for i in range(n):
                    barcode = acts.examples.SimBarcode()
                    barcode.vertexPrimary = int(d["vertex_primary"][idx][i])
                    barcode.vertexSecondary = int(d["vertex_secondary"][idx][i])
                    barcode.particle = int(d["particle"][idx][i])
                    barcode.generation = int(d["generation"][idx][i])
                    barcode.subParticle = int(d["sub_particle"][idx][i])

                    p = acts.examples.SimParticle(
                        barcode,
                        acts.PdgParticle(int(d["particle_type"][idx][i])),
                        float(d["q"][idx][i]) * u.e,
                        float(d["m"][idx][i]) * u.GeV,
                    )
                    p.process = acts.examples.ProcessType(int(d["process"][idx][i]))

                    px = float(d["px"][idx][i])
                    py = float(d["py"][idx][i])
                    pz = float(d["pz"][idx][i])
                    p.setInitialPosition4(
                        acts.Vector4(
                            float(d["vx"][idx][i]) * u.mm,
                            float(d["vy"][idx][i]) * u.mm,
                            float(d["vz"][idx][i]) * u.mm,
                            float(d["vt"][idx][i]) * u.mm,
                        )
                    )
                    p.setInitialDirection(acts.Vector3(px, py, pz))
                    p.setInitialAbsoluteMomentum(float(d["p"][idx][i]) * u.GeV)

                    p.setFinalMaterialPassed(
                        float(d["total_x0"][idx][i]) * u.mm,
                        float(d["total_l0"][idx][i]) * u.mm,
                    )
                    p.setFinalNumberOfHits(int(d["number_of_hits"][idx][i]))
                    p.setFinalOutcome(
                        acts.examples.ParticleOutcome(int(d["outcome"][idx][i]))
                    )

                    particles.insert(p)

            self._particleHandle(context, particles)

        if self._hit_data is not None:
            hits = acts.examples.SimHitContainer()
            if event_number in self._hit_event_range_map:
                start, end = self._hit_event_range_map[event_number]
                d = self._hit_data
                u = acts.UnitConstants
                for i in range(start, end):
                    geoid = acts.GeometryIdentifier(int(d["geometry_id"][i]))

                    if self._has_barcode_branch:
                        bc_data = d["barcode"][i]
                        barcode = acts.examples.SimBarcode()
                        barcode.vertexPrimary = int(bc_data[0])
                        barcode.vertexSecondary = int(bc_data[1])
                        barcode.particle = int(bc_data[2])
                        barcode.generation = int(bc_data[3])
                        barcode.subParticle = int(bc_data[4])
                    else:
                        barcode = acts.examples.SimBarcode()
                        barcode.vertexPrimary = int(d["barcode_vertex_primary"][i])
                        barcode.vertexSecondary = int(
                            d["barcode_vertex_secondary"][i]
                        )
                        barcode.particle = int(d["barcode_particle"][i])
                        barcode.generation = int(d["barcode_generation"][i])
                        barcode.subParticle = int(d["barcode_sub_particle"][i])

                    pos4 = acts.Vector4(
                        float(d["tx"][i]) * u.mm,
                        float(d["ty"][i]) * u.mm,
                        float(d["tz"][i]) * u.mm,
                        float(d["tt"][i]) * u.mm,
                    )
                    before4 = acts.Vector4(
                        float(d["tpx"][i]) * u.GeV,
                        float(d["tpy"][i]) * u.GeV,
                        float(d["tpz"][i]) * u.GeV,
                        float(d["te"][i]) * u.GeV,
                    )
                    after4 = acts.Vector4(
                        (float(d["tpx"][i]) + float(d["deltapx"][i])) * u.GeV,
                        (float(d["tpy"][i]) + float(d["deltapy"][i])) * u.GeV,
                        (float(d["tpz"][i]) + float(d["deltapz"][i])) * u.GeV,
                        (float(d["te"][i]) + float(d["deltae"][i])) * u.GeV,
                    )
                    index = int(d["index"][i])

                    hit = acts.examples.SimHit(
                        geoid, barcode, pos4, before4, after4, index
                    )
                    hits.insert(hit)

            self._hitHandle(context, hits)

        return acts.examples.ProcessCode.SUCCESS
