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
        customSeederConfig={
            "model_path": "/path/to/model.onnx",
            "r_max": 500.0,
            "z_max": 1000.0,
            "score_threshold": 0.35,
            "providers": ["CUDAExecutionProvider", "CPUExecutionProvider"],
        },
        ...
    )
"""

import numpy as np
import onnxruntime as ort

import acts
import acts.examples

# Numerical stability guard: seeds where two spacepoints are identical in the transverse
# plane crash the CKF. This is not a physics threshold, it exists solely to prevent
# crashes from pathological seeds.
_SP_DEDUP_GUARD_MM = 1e-3


def _apply_model_acceptance(
    sp, xyz: np.ndarray, r_max: float, z_max: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return (filtered_sp_xyz, original_indices) after applying the acceptance cut.

    Indices in original_indices map positions in filtered_sp_xyz back to rows in xyz.
    sp.r is read directly from the stored column rather than recomputed from x, y.
    """
    r = np.asarray(sp.r)
    mask = (r < r_max) & (np.abs(xyz[:, 2]) < z_max)
    return xyz[mask].astype(np.float32), np.where(mask)[0].astype(np.uint32)


def _filter_valid_seeds(
    seeds: np.ndarray, scores: np.ndarray, sp_xyz: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Keep only seeds with a finite score, all three spacepoint indices filled, and no
    duplicate indices. When sp_xyz is provided, also drops seeds with degenerate
    pairwise transverse separation (crash guard, see _SP_DEDUP_GUARD_MM)."""
    valid = np.isfinite(scores) & (np.sum(seeds >= 0, axis=1) == 3)
    seeds_v = seeds[valid]
    scores_v = scores[valid]
    if len(seeds_v) > 0:
        i0, i1, i2 = seeds_v[:, 0], seeds_v[:, 1], seeds_v[:, 2]
        no_dups = (i0 != i1) & (i1 != i2) & (i0 != i2)
        seeds_v = seeds_v[no_dups]
        scores_v = scores_v[no_dups]
    if sp_xyz is not None and len(seeds_v) > 0:
        i0, i1, i2 = seeds_v[:, 0], seeds_v[:, 1], seeds_v[:, 2]
        xy = sp_xyz[:, :2]
        d01 = np.linalg.norm(xy[i0] - xy[i1], axis=1)
        d12 = np.linalg.norm(xy[i1] - xy[i2], axis=1)
        d02 = np.linalg.norm(xy[i0] - xy[i2], axis=1)
        non_degen = np.minimum(np.minimum(d01, d12), d02) >= _SP_DEDUP_GUARD_MM
        seeds_v = seeds_v[non_degen]
        scores_v = scores_v[non_degen]
    return seeds_v, scores_v


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
        r_max: float = 500.0,
        z_max: float = 1000.0,
        score_threshold: float = 0.35,
        providers: list[str] | None = None,
    ):
        acts.examples.IAlgorithm.__init__(self, "GuntamTransformerSeeder", log_level)

        self._r_max = r_max
        self._z_max = z_max
        self._score_threshold = score_threshold

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
            providers=providers or ["CPUExecutionProvider"],
        )

    def execute(self, ctx) -> acts.examples.ProcessCode:
        sp = self._sp_handle(ctx.eventStore)
        xyz = np.stack([np.asarray(sp.x), np.asarray(sp.y), np.asarray(sp.z)], axis=1)

        filtered_sp_xyz, orig_idx = _apply_model_acceptance(
            sp, xyz, self._r_max, self._z_max
        )

        seeds_raw, scores_raw = self._session.run(
            ["seeds", "seed_scores"],
            {"hits": filtered_sp_xyz},
        )

        seeds_v, scores_v = _filter_valid_seeds(seeds_raw, scores_raw, filtered_sp_xyz)

        score_mask = scores_v > self._score_threshold
        seeds_v = seeds_v[score_mask]
        scores_v = scores_v[score_mask]

        container = acts.SeedContainer2()
        container.assignSpacePointContainer(sp)

        sp_indices_v = orig_idx[seeds_v]

        seed_proxies = [container.createSeed() for _ in range(len(seeds_v))]
        for seed, sp_indices, score in zip(seed_proxies, sp_indices_v, scores_v):
            seed.quality = float(score)
            seed.vertexZ = 0.0
            seed.assignSpacePointIndices(sp_indices.tolist())

        # MutableSeedProxy2 holds raw pointers into the container; drop them before
        # the whiteboard write which transfers container ownership to C++.
        del seed_proxies

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
        ``model_path``      (str, required)            — path to the ONNX model file.
        ``num_threads``     (int, default 1)            — ONNX intra/inter-op thread count.
            When running events in parallel (Sequencer.numThreads > 1), ensure
            numThreads * num_threads does not exceed available CPU cores.
        ``r_max``           (float, default 500.0)     — radial acceptance cut in mm.
        ``z_max``           (float, default 1000.0)    — longitudinal acceptance cut in mm.
        ``score_threshold`` (float, default 0.35)      — minimum score to keep a seed.
        ``providers``       (list[str], default None)  — ONNX execution providers in
            priority order. None falls back to ["CPUExecutionProvider"]. Pass
            ["CUDAExecutionProvider", "CPUExecutionProvider"] for GPU with CPU fallback.

    **kwargs absorbs trackingGeometry, logLevel, and any future addSeeding additions.
    Returns the whiteboard key for the output seeds.
    """
    try:
        model_path = config["model_path"]
    except KeyError as exc:
        raise ValueError("customSeederConfig must contain 'model_path'") from exc
    num_threads = int(config.get("num_threads", 1))
    r_max = float(config.get("r_max", 500.0))
    z_max = float(config.get("z_max", 1000.0))
    score_threshold = float(config.get("score_threshold", 0.35))
    providers = config.get("providers", None)
    log_level = kwargs.get("logLevel", acts.logging.INFO)
    s.addAlgorithm(
        _GuntamAlgorithm(
            model_path,
            spacePoints,
            outputSeeds,
            log_level,
            num_threads=num_threads,
            r_max=r_max,
            z_max=z_max,
            score_threshold=score_threshold,
            providers=providers,
        )
    )
    return outputSeeds
