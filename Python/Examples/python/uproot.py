# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Uproot-based readers for simulated particles and hits.

Reads the ROOT files written by RootParticleWriter and RootSimHitWriter
using uproot, without requiring the ROOT plugin. Suitable for use in the
PyPI distribution.
"""

from pathlib import Path
from typing import Union

import numpy as np
import uproot

import acts
import acts.examples


class UprootParticleReader(acts.examples.IReader):
    """Reads simulated particles from a ROOT file using uproot.

    Reads the file format produced by RootParticleWriter (one entry per event,
    vector-valued branches). All data is loaded upfront for thread-safe
    concurrent access from the sequencer.
    """

    def __init__(
        self,
        filePath: Union[Path, str],
        outputParticles: str = "particles",
        level: acts.logging.Level = acts.logging.INFO,
    ):
        acts.examples.IReader.__init__(self, "UprootParticleReader", level)

        self._outputParticles = outputParticles

        self._particleHandle = acts.examples.WriteDataHandle(
            self, acts.examples.SimParticleContainer, "OutputParticles"
        )

        with uproot.open(str(filePath)) as f:
            tree = f["particles"]
            self._data = tree.arrays(library="np")
            event_ids = self._data["event_id"]
            self._entry_map = {int(eid): i for i, eid in enumerate(event_ids)}

        self._min_event = min(self._entry_map.keys())
        self._max_event = max(self._entry_map.keys())

    def availableEvents(self):
        return (self._min_event, self._max_event + 1)

    def initialize(self):
        self._particleHandle.initialize(self._outputParticles)
        return acts.examples.ProcessCode.SUCCESS

    def read(self, context):
        event_number = context.eventNumber
        particles = acts.examples.SimParticleContainer()

        if event_number in self._entry_map:
            idx = self._entry_map[event_number]
            d = self._data
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
                p.process = acts.examples.GenerationProcess(int(d["process"][idx][i]))

                p.fourPosition = acts.Vector4(
                    *[float(d[k][idx][i]) * u.mm for k in ("vx", "vy", "vz", "vt")]
                )
                p.direction = acts.Vector3(
                    *[float(d[k][idx][i]) for k in ("px", "py", "pz")]
                )
                p.absoluteMomentum = float(d["p"][idx][i]) * u.GeV

                p.setFinalMaterialPassed(
                    float(d["total_x0"][idx][i]) * u.mm,
                    float(d["total_l0"][idx][i]) * u.mm,
                )
                p.numberOfHits = int(d["number_of_hits"][idx][i])
                p.outcome = acts.examples.SimulationOutcome(int(d["outcome"][idx][i]))

                particles.insert(p)

        self._particleHandle(context, particles)
        return acts.examples.ProcessCode.SUCCESS


class UprootSimHitReader(acts.examples.IReader):
    """Reads simulated hits from a ROOT file using uproot.

    Reads the file format produced by RootSimHitWriter (one row per hit,
    scalar branches grouped by event_id). All data is loaded upfront for
    thread-safe concurrent access from the sequencer.
    """

    def __init__(
        self,
        filePath: Union[Path, str],
        outputSimHits: str = "simhits",
        level: acts.logging.Level = acts.logging.INFO,
    ):
        acts.examples.IReader.__init__(self, "UprootSimHitReader", level)

        self._outputSimHits = outputSimHits

        self._hitHandle = acts.examples.WriteDataHandle(
            self, acts.examples.SimHitContainer, "OutputSimHits"
        )

        with uproot.open(str(filePath)) as f:
            tree = f["hits"]
            self._data = tree.arrays(library="np")
            self._event_range_map = self._build_event_range_map(self._data["event_id"])

        all_ids = set(self._event_range_map.keys())
        self._min_event = int(min(all_ids))
        self._max_event = int(max(all_ids))

    @staticmethod
    def _build_event_range_map(event_ids: np.ndarray) -> dict:
        """Build a {event_id: (start, end)} map from a sorted event_id array."""
        if len(event_ids) == 0:
            return {}
        unique_ids, starts = np.unique(event_ids, return_index=True)
        ends = np.append(starts[1:], len(event_ids))
        return {
            int(uid): (int(s), int(e)) for uid, s, e in zip(unique_ids, starts, ends)
        }

    def availableEvents(self):
        return (self._min_event, self._max_event + 1)

    def initialize(self):
        self._hitHandle.initialize(self._outputSimHits)
        return acts.examples.ProcessCode.SUCCESS

    def read(self, context):
        event_number = context.eventNumber
        hits = acts.examples.SimHitContainer()

        if event_number in self._event_range_map:
            start, end = self._event_range_map[event_number]
            d = self._data
            u = acts.UnitConstants
            for i in range(start, end):
                geoid = acts.GeometryIdentifier(int(d["geometry_id"][i]))

                barcode = acts.examples.SimBarcode()
                barcode.vertexPrimary = int(d["barcode_vertex_primary"][i])
                barcode.vertexSecondary = int(d["barcode_vertex_secondary"][i])
                barcode.particle = int(d["barcode_particle"][i])
                barcode.generation = int(d["barcode_generation"][i])
                barcode.subParticle = int(d["barcode_sub_particle"][i])

                pos4 = acts.Vector4(
                    *[float(d[k][i]) * u.mm for k in ("tx", "ty", "tz", "tt")]
                )
                before4 = acts.Vector4(
                    *[float(d[k][i]) * u.GeV for k in ("tpx", "tpy", "tpz", "te")]
                )
                after4 = acts.Vector4(
                    (float(d["tpx"][i]) + float(d["deltapx"][i])) * u.GeV,
                    (float(d["tpy"][i]) + float(d["deltapy"][i])) * u.GeV,
                    (float(d["tpz"][i]) + float(d["deltapz"][i])) * u.GeV,
                    (float(d["te"][i]) + float(d["deltae"][i])) * u.GeV,
                )

                hit = acts.examples.SimHit(
                    geoid, barcode, pos4, before4, after4, int(d["index"][i])
                )
                hits.insert(hit)

        self._hitHandle(context, hits)
        return acts.examples.ProcessCode.SUCCESS
