"""
guntam_transformer_seeder.py

Integrates the GUNTAM ONNX transformer into ACTS as a PythonCallable seeder.

Requires ``acts.SpacePointContainer2``, ``acts.SeedContainer2``, and the
associated Python bindings for mutable space points and seeds.

Usage:
    from guntam_transformer_seeder import guntam_transformer_seeder

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.PythonCallable,
        customSeeder=guntam_transformer_seeder,
        customSeederConfig={"model_path": "/path/to/model.onnx"},
        ...
    )
"""

import numpy as np
import onnxruntime as ort

import acts
import acts.examples

# Geometry acceptance window matching the GUNTAM training domain.
_R_MAX = 500.0   # mm
_Z_MAX = 1000.0  # mm

# Avoid z-degenerate triplets, which can make track-parameter estimation ill-defined.
_Z_EPS = 0.01    # mm


def _apply_model_acceptance(sp, xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (hits_filtered, original_indices) after applying the acceptance cut.

    Indices in original_indices map positions in hits_filtered back to rows in xyz.
    sp.r is read directly from the stored column rather than recomputed from x, y.
    """
    r = np.asarray(sp.r)
    mask = (r < _R_MAX) & (np.abs(xyz[:, 2]) < _Z_MAX)
    return xyz[mask].astype(np.float32), np.where(mask)[0].astype(np.int64)


def _filter_valid_seeds(
    seeds: np.ndarray, scores: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Keep only seeds with a finite score and all three hit indices filled."""
    valid = np.isfinite(scores) & (np.sum(seeds >= 0, axis=1) == 3)
    return seeds[valid], scores[valid]


def _filter_degenerate_triplets(
    seeds: np.ndarray, hits: np.ndarray
) -> np.ndarray:
    """Return a boolean mask for triplets that pass geometric validity checks.

    Rejects seeds with out-of-bounds indices, duplicate hit indices, or any
    two z-coordinates within _Z_EPS mm of each other.
    """
    n = len(hits)

    in_bounds = np.all((seeds >= 0) & (seeds < n), axis=1)
    seeds_safe = seeds[in_bounds]

    if len(seeds_safe) == 0:
        return in_bounds

    i0, i1, i2 = seeds_safe[:, 0], seeds_safe[:, 1], seeds_safe[:, 2]
    no_dups = (i0 != i1) & (i1 != i2) & (i0 != i2)

    z = hits[:, 2].astype(np.float64)
    z_ok = (
        (np.abs(z[i1] - z[i0]) > _Z_EPS)
        & (np.abs(z[i2] - z[i1]) > _Z_EPS)
        & (np.abs(z[i2] - z[i0]) > _Z_EPS)
    )

    good = np.zeros(len(seeds), dtype=bool)
    good[in_bounds] = no_dups & z_ok
    return good


class _GuntamAlgorithm(acts.examples.IAlgorithm):
    """Per-event algorithm: reads spacepoints, runs GUNTAM ONNX, writes seeds.

    The ONNX session is created once at construction. Configure num_threads
    conservatively when running the sequencer with multiple event threads.
    """

    def __init__(
        self,
        model_path: str,
        sp_key: str,
        seeds_key: str,
        log_level,
        num_threads: int = 1,
    ):
        acts.examples.IAlgorithm.__init__(self, "GuntamTransformerSeeder", log_level)

        self._sp_handle = acts.examples.ReadDataHandle(
            self, acts.SpacePointContainer2, "InputSpacePoints"
        )
        self._sp_handle.initialize(sp_key)

        self._seeds_handle = acts.examples.WriteDataHandle(
            self, acts.SeedContainer2, "OutputSeeds"
        )
        self._seeds_handle.initialize(seeds_key)

        sess_options = ort.SessionOptions()
        sess_options.intra_op_num_threads = num_threads
        sess_options.inter_op_num_threads = num_threads

        self._session = ort.InferenceSession(
            model_path,
            sess_options=sess_options,
            providers=["CPUExecutionProvider"],
        )

    def execute(self, ctx) -> acts.examples.ProcessCode:
        sp = self._sp_handle(ctx.eventStore)
        xyz = np.stack(
            [np.asarray(sp.x), np.asarray(sp.y), np.asarray(sp.z)], axis=1
        )

        hits, orig_idx = _apply_model_acceptance(sp, xyz)

        seeds_raw, scores_raw = self._session.run(
            ["seeds", "seed_scores"],
            {"hits": hits},
        )

        seeds_v, scores_v = _filter_valid_seeds(seeds_raw, scores_raw)
        geo_mask = _filter_degenerate_triplets(seeds_v, hits)
        seeds_v = seeds_v[geo_mask]
        scores_v = scores_v[geo_mask]

        container = acts.SeedContainer2()
        container.assignSpacePointContainer(sp)
        seed = None
        for triplet, score in zip(seeds_v, scores_v):
            sp_indices = [int(orig_idx[j]) for j in triplet]
            seed = container.createSeed()
            seed.quality = float(score)
            seed.vertexZ = 0.0
            seed.assignSpacePointIndices(sp_indices)

        # MutableSeedProxy2 holds raw pointers into the container; drop it before
        # the whiteboard write which transfers container ownership to C++.
        del seed

        self._seeds_handle(ctx, container)
        return acts.examples.ProcessCode.SUCCESS


def guntam_transformer_seeder(
    s,
    spacePoints: str,
    outputSeeds: str,
    config: dict,
    **kwargs,
) -> str:
    """PythonCallable entry point for the GUNTAM ONNX seeder.

    config keys:
        ``model_path``  (str, required) — path to the ONNX model file.
        ``num_threads`` (int, default 1) — ONNX intra- and inter-op thread count.
            When running events in parallel (Sequencer.numThreads > 1), ensure
            numThreads * num_threads does not exceed available CPU cores.

    **kwargs absorbs trackingGeometry, logLevel, and any future addSeeding additions.
    Returns the whiteboard key for the output seeds.
    """
    try:
        model_path = config["model_path"]
    except KeyError as exc:
        raise ValueError("customSeederConfig must contain 'model_path'") from exc
    num_threads = int(config.get("num_threads", 1))
    log_level = kwargs.get("logLevel", acts.logging.INFO)
    s.addAlgorithm(
        _GuntamAlgorithm(model_path, spacePoints, outputSeeds, log_level, num_threads)
    )
    return outputSeeds
