from pathlib import Path
from math import sqrt
import sys, os
import acts
import acts.examples


def getOpenDataDetector(
    odd_dir: Path,
    mdecorator=None,
    logLevel=acts.logging.INFO,
):

    import acts.examples.dd4hep

    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    odd_xml = odd_dir / "xml" / "OpenDataDetector.xml"
    if not odd_xml.exists():
        raise RuntimeError(f"OpenDataDetector.xml not found at {odd_xml}")

    env_vars = []
    map_name = "libOpenDataDetector.components"
    lib_name = None
    if sys.platform == "linux":
        env_vars = ["LD_LIBRARY_PATH"]
        lib_name = "libOpenDataDetector.so"
    elif sys.platform == "darwin":
        env_vars = ["DYLD_LIBRARY_PATH", "DD4HEP_LIBRARY_PATH"]
        lib_name = "libOpenDataDetector.dylib"

    if lib_name is not None and len(env_vars) > 0:
        found = False
        for env_var in env_vars:
            for lib_dir in os.environ.get(env_var, "").split(":"):
                lib_dir = Path(lib_dir)
                if (lib_dir / map_name).exists() and (lib_dir / lib_name).exists():
                    found = True
                    break
        if not found:
            msg = (
                "Unable to find OpenDataDetector factory library. "
                f"You might need to point {'/'.join(env_vars)} to build/thirdparty/OpenDataDetector/factory or other ODD install location"
            )
            raise RuntimeError(msg)

    volumeRadiusCutsMap = {
        28: [850.0],  # LStrip negative z
        30: [850.0],  # LStrip positive z
        23: [400.0, 550.0],  # SStrip negative z
        25: [400.0, 550.0],  # SStrip positive z
        16: [100.0],  # Pixels negative z
        18: [100.0],  # Pixels positive z
    }

    def geoid_hook(geoid, surface):
        if geoid.volume() in volumeRadiusCutsMap:
            r = sqrt(surface.center()[0] ** 2 + surface.center()[1] ** 2)

            geoid.setExtra(1)
            for cut in volumeRadiusCutsMap[geoid.volume()]:
                if r > cut:
                    geoid.setExtra(geoid.extra() + 1)

        return geoid

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_xml)],
        logLevel=customLogLevel(),
        dd4hepLogLevel=customLogLevel(),
        geometryIdentifierHook=acts.GeometryIdentifierHook(geoid_hook),
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    if mdecorator is None:
        mdecorator = acts.JsonMaterialDecorator(
            rConfig=config,
            jFileName=str(odd_dir / "config/odd-material-mapping-config.json"),
            level=customLogLevel(minLevel=acts.logging.WARNING),
        )

    trackingGeometry, deco = detector.finalize(dd4hepConfig, mdecorator)

    return detector, trackingGeometry, deco
