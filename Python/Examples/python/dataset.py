import os
from pathlib import Path
from typing import Optional

import acts
import acts.examples

#: HuggingFace dataset id of the ColliderML Release 1 dataset, matching
#: `colliderml.core.hf_download.DEFAULT_DATASET_ID` upstream.
COLLIDERML_DATASET_ID = "CERN/ColliderML-Release-1"


def getColliderMLDirectory() -> Path:
    """Return the local ColliderML cache directory.

    Uses `$COLLIDERML_DATA_DIR` if set, otherwise falls back to
    `~/.cache/colliderml`. This matches the canonical `colliderml` Python
    library (`colliderml.core.hf_download.default_data_dir`), so a cache
    populated with `colliderml download ...` is found without any further
    configuration.

    Raises `RuntimeError` if the resulting directory does not exist.
    """
    env = os.environ.get("COLLIDERML_DATA_DIR")
    if env:
        path = Path(env).expanduser().resolve()
    else:
        path = (Path.home() / ".cache" / "colliderml").resolve()

    if not path.is_dir():
        raise RuntimeError(
            f"ColliderML data directory not found at {path}. "
            f"Set $COLLIDERML_DATA_DIR to point at an existing ColliderML "
            f"cache, or populate the default location with "
            f"'colliderml download ...'."
        )
    return path


def _sanitizedDatasetId(datasetId: str) -> str:
    return datasetId.replace("/", "__")


def getColliderMLObjectDirectory(
    object: str,
    channel: str = "ttbar",
    pileup: str = "pu0",
    dataDir: Optional[Path] = None,
    datasetId: str = COLLIDERML_DATASET_ID,
) -> Path:
    """Resolve the local directory holding one ColliderML "object" (e.g.
    `particles`, `tracker_hits`, `tracks`) for a given channel/pileup.

    Parameters
    ----------
    object : str
        Object name, e.g. "particles", "tracker_hits", "tracks".
    channel : str
        Physics channel/process, e.g. "ttbar".
    pileup : str
        Pileup token, e.g. "pu0" or "pu200".
    dataDir : Path, None
        Base cache directory. If not given, uses `getColliderMLDirectory()`.
    datasetId : str
        HuggingFace dataset id the config was downloaded from.

    A directory downloaded via the `colliderml` library's `download_config`
    lands at `<dataDir>/<datasetId sanitized>/<config>/data/<config>`. Some
    hand-prepared samples (e.g. the ACTS CI sample) omit the dataset-id
    level, so both layouts are tried.
    """
    if dataDir is None:
        dataDir = getColliderMLDirectory()
    else:
        dataDir = Path(dataDir)

    config = f"{channel}_{pileup}_{object}"
    candidates = [
        dataDir / _sanitizedDatasetId(datasetId) / config / "data" / config,
        dataDir / config / "data" / config,
    ]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate

    searched = "\n".join(f"  {c}" for c in candidates)
    raise RuntimeError(
        f"ColliderML config '{config}' not found. Searched:\n{searched}\n"
        f"Download it with:\n"
        f"  colliderml download --channels {channel} --pileup {pileup} --objects {object}"
    )


def addColliderML(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    channel: str = "ttbar",
    pileup: str = "pu0",
    dataDir: Optional[Path] = None,
    particlesDir: Optional[Path] = None,
    hitsDir: Optional[Path] = None,
    tracksDir: Optional[Path] = None,
    readTracks: bool = True,
    geoIdMapPath: Optional[Path] = None,
    geoIdMapSourcePrefix: str = "gen1",
    geoIdMapTargetPrefix: str = "gen3",
    outputParticles: str = "particles",
    outputSimHits: str = "simhits",
    outputMeasurements: str = "measurements",
    outputMeasurementSubset: str = "measurement_subset",
    outputMeasSimHitsMap: str = "measurement_simhits_map",
    outputMeasParticlesMap: str = "measurement_particles_map",
    outputParticleMeasurementsMap: str = "particle_measurements_map",
    outputTracks: str = "colliderml_tracks",
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.Sequencer:
    """Read a ColliderML sample and convert it to ACTS EDM collections.

    Parameters
    ----------
    s : Sequencer
        the sequencer to add the ColliderML reading steps to
    trackingGeometry : TrackingGeometry
        tracking geometry to build measurements against
    channel, pileup : str
        selects the sample, e.g. channel="ttbar", pileup="pu0". Only used to
        resolve `particlesDir`/`hitsDir`/`tracksDir` via
        `getColliderMLObjectDirectory` when those are not given explicitly.
    dataDir : Path, None
        base ColliderML cache directory, forwarded to
        `getColliderMLObjectDirectory`. None triggers `getColliderMLDirectory()`.
    particlesDir, hitsDir, tracksDir : Path, None
        raw overrides for the resolved directories, bypassing channel/pileup/
        dataDir resolution entirely.
    readTracks : bool
        if True, also read and convert the dataset's own published tracks
        (output under `outputTracks`, default "colliderml_tracks"). Not all
        samples ship these.
    geoIdMapPath : Path, None
        CSV mapping detector geometry ids to `trackingGeometry` geometry ids,
        e.g. produced by `generate_geoid_map.py`. None uses `trackingGeometry`
        directly.
    """
    from acts.examples.arrow import (
        ColliderMLRelease1InputConverter,
        ParquetReader,
    )

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if particlesDir is None:
        particlesDir = getColliderMLObjectDirectory(
            "particles", channel=channel, pileup=pileup, dataDir=dataDir
        )
    if hitsDir is None:
        hitsDir = getColliderMLObjectDirectory(
            "tracker_hits", channel=channel, pileup=pileup, dataDir=dataDir
        )
    if readTracks and tracksDir is None:
        tracksDir = getColliderMLObjectDirectory(
            "tracks", channel=channel, pileup=pileup, dataDir=dataDir
        )

    readerCollections = {
        "cml_particles": str(particlesDir),
        "cml_hits": str(hitsDir),
    }
    readerSchemas = {
        "cml_particles": ColliderMLRelease1InputConverter.particleSchema(),
        "cml_hits": ColliderMLRelease1InputConverter.hitSchema(),
    }
    if readTracks:
        readerCollections["cml_tracks"] = str(tracksDir)
        readerSchemas["cml_tracks"] = ColliderMLRelease1InputConverter.tracksSchema()

    s.addReader(
        ParquetReader(
            level=customLogLevel(),
            collections=readerCollections,
            expectedSchemas=readerSchemas,
        )
    )

    converterKWArgs = dict(
        inputParticlesTable="cml_particles",
        inputHitsTable="cml_hits",
        outputParticles=outputParticles,
        outputSimHits=outputSimHits,
        outputMeasurements=outputMeasurements,
        outputMeasurementSubset=outputMeasurementSubset,
        outputMeasSimHitsMap=outputMeasSimHitsMap,
        outputMeasParticlesMap=outputMeasParticlesMap,
        outputParticleMeasurementsMap=outputParticleMeasurementsMap,
        trackingGeometry=trackingGeometry,
    )
    if geoIdMapPath is not None:
        converterKWArgs["geoIdMapPath"] = geoIdMapPath
        converterKWArgs["geoIdMapSourcePrefix"] = geoIdMapSourcePrefix
        converterKWArgs["geoIdMapTargetPrefix"] = geoIdMapTargetPrefix
    if readTracks:
        converterKWArgs["inputTracksTable"] = "cml_tracks"
        converterKWArgs["outputTracks"] = outputTracks

    s.addAlgorithm(
        ColliderMLRelease1InputConverter(level=customLogLevel(), **converterKWArgs)
    )

    return s
