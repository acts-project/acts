import acts
import acts.examples

from pathlib import Path


def getOpenDataDetectorDirectory():
    return (
        Path(__file__).parent.parent.parent.parent / "thirdparty" / "OpenDataDetector"
    )


def getOpenDataDetector(mdecorator=None):
    import acts.examples.dd4hep

    odd_dir = getOpenDataDetectorDirectory()

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[str(odd_dir / "xml/OpenDataDetector.xml")]
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
