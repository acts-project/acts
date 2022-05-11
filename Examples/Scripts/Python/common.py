import acts
import acts.examples
import os
import sys


from pathlib import Path


def getOpenDataDetectorDirectory():
    return (
        Path(__file__).parent.parent.parent.parent / "thirdparty" / "OpenDataDetector"
    )


def getOpenDataDetector(mdecorator=None):
    odd_dir = getOpenDataDetectorDirectory()

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
        lib_name = "libOpenDataDetector.dyld"

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
                f"You might need to point {'/'.join(env_vars)} at it"
            )
            raise RuntimeError(msg)

    import acts.examples.dd4hep

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_xml)]
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    if mdecorator is None:
        mdecorator = acts.JsonMaterialDecorator(
            rConfig=config,
            jFileName=str(odd_dir / "config/odd-material-mapping-config.json"),
            level=acts.logging.WARNING,
        )

    trackingGeometry, deco = detector.finalize(dd4hepConfig, mdecorator)

    return detector, trackingGeometry, deco


def addPythia8(
    sequencer: acts.examples.Sequencer,
    rnd: acts.examples.RandomNumbers,
    nhard=1,
    npileup=200,
    beam0=acts.PdgParticle.eProton,
    beam1=acts.PdgParticle.eProton,
    cmsEnergy=14 * acts.UnitConstants.TeV,
    vertexStddev: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
    vertexMean: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
):
    """This function steers the particle generation using Pythia8

    NB. this version is included here only for compatibility. Please use pythia8.addPythia8 instead.
    """
    import pythia8

    evGen = pythia8.addPythia8(
        sequencer,
        rnd=rnd,
        nhard=nhard,
        npileup=npileup,
        beam=(beam0, beam1),
        cmsEnergy=cmsEnergy,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=vertexStddev, mean=vertexMean
        ),
        returnEvGen=True,
    )

    return evGen
