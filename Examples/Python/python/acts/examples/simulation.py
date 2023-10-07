from typing import Optional, Union, Any, List
from pathlib import Path
from collections import namedtuple
from collections.abc import Iterable

import acts
from acts.examples import (
    RandomNumbers,
    EventGenerator,
    FixedMultiplicityGenerator,
    CsvParticleWriter,
    ParticlesPrinter,
    RootParticleWriter,
)

# Defaults (given as `None` here) use class defaults defined in
# Examples/Algorithms/Generators/ActsExamples/Generators/ParametricParticleGenerator.hpp
MomentumConfig = namedtuple(
    "MomentumConfig",
    ["min", "max", "transverse"],
    defaults=[None, None, None],
)
EtaConfig = namedtuple(
    "EtaConfig", ["min", "max", "uniform"], defaults=[None, None, None]
)
PhiConfig = namedtuple("PhiConfig", ["min", "max"], defaults=[None, None])
ParticleConfig = namedtuple(
    "ParticleConfig",
    ["num", "pdg", "randomizeCharge", "charge", "mass"],
    defaults=[None, None, None, None, None],
)
ParticleSelectorConfig = namedtuple(
    "ParticleSelectorConfig",
    [
        "rho",  # (min,max)
        "absZ",  # (min,max)
        "time",  # (min,max)
        "phi",  # (min,max)
        "eta",  # (min,max)
        "absEta",  # (min,max)
        "pt",  # (min,max)
        "removeCharged",  # bool
        "removeNeutral",  # bool
    ],
    defaults=[(None, None)] * 7 + [None] * 2,
)


@acts.examples.NamedTypeArgs(
    momentumConfig=MomentumConfig,
    etaConfig=EtaConfig,
    phiConfig=PhiConfig,
    particleConfig=ParticleConfig,
)
def addParticleGun(
    s: acts.examples.Sequencer,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    momentumConfig: MomentumConfig = MomentumConfig(),
    etaConfig: EtaConfig = EtaConfig(),
    phiConfig: PhiConfig = PhiConfig(),
    particleConfig: ParticleConfig = ParticleConfig(),
    multiplicity: int = 1,
    vtxGen: Optional[EventGenerator.VertexGenerator] = None,
    printParticles: bool = False,
    rnd: Optional[RandomNumbers] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the particle generation using the particle gun

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    momentumConfig : MomentumConfig(min, max, transverse)
        momentum configuration: minimum momentum, maximum momentum, transverse
    etaConfig : EtaConfig(min, max, uniform)
        pseudorapidity configuration: eta min, eta max, uniform
    phiConfig : PhiConfig(min, max)
        azimuthal angle configuration: phi min, phi max
    particleConfig : ParticleConfig(num, pdg, randomizeCharge, charge, mass)
        particle configuration: number of particles, particle type, charge flip
    multiplicity : int, 1
        number of generated vertices
    vtxGen : VertexGenerator, None
        vertex generator module
    printParticles : bool, False
        print generated particles
    rnd : RandomNumbers, None
        random number generator
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Preliminaries
    rnd = rnd or RandomNumbers(seed=228)

    # Input
    evGen = EventGenerator(
        level=customLogLevel(),
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=multiplicity),
                vertex=vtxGen
                or acts.examples.GaussianVertexGenerator(
                    mean=acts.Vector4(0, 0, 0, 0),
                    stddev=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    **acts.examples.defaultKWArgs(
                        p=(momentumConfig.min, momentumConfig.max),
                        pTransverse=momentumConfig.transverse,
                        eta=(etaConfig.min, etaConfig.max),
                        phi=(phiConfig.min, phiConfig.max),
                        etaUniform=etaConfig.uniform,
                        numParticles=particleConfig.num,
                        pdg=particleConfig.pdg,
                        randomizeCharge=particleConfig.randomizeCharge,
                        charge=particleConfig.charge,
                        mass=particleConfig.mass,
                    )
                ),
            )
        ],
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    s.addWhiteboardAlias("particles", evGen.config.outputParticles)

    if printParticles:
        s.addAlgorithm(
            ParticlesPrinter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            CsvParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            RootParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "particles.root"),
            )
        )

    return s


def addPythia8(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    nhard: int = 1,
    npileup: int = 200,
    beam: Optional[
        Union[acts.PdgParticle, Iterable]
    ] = None,  # default: acts.PdgParticle.eProton
    cmsEnergy: Optional[float] = None,  # default: 14 * acts.UnitConstants.TeV
    hardProcess: Optional[Iterable] = None,  # default: ["HardQCD:all = on"]
    pileupProcess: Iterable = ["SoftQCD:all = on"],
    vtxGen: Optional[acts.examples.EventGenerator.VertexGenerator] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    printParticles: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the particle generation using Pythia8

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    rnd : RandomNumbers, None
        random number generator
    nhard, npileup : int, 1, 200
        Number of hard-scatter and pileup vertices
    beam : PdgParticle|[PdgParticle,PdgParticle], eProton
        beam particle(s)
    cmsEnergy : float, 14 TeV
        CMS energy
    hardProcess, pileupProcess : [str], ["HardQCD:all = on"], ["SoftQCD:all = on"]
        hard and pileup processes
    vtxGen : VertexGenerator, None
        vertex generator module
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    printParticles : bool, False
        print generated particles
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Preliminaries
    rnd = rnd or acts.examples.RandomNumbers()
    vtxGen = vtxGen or acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
    )
    if not isinstance(beam, Iterable):
        beam = (beam, beam)

    generators = []
    if nhard is not None and nhard > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=customLogLevel(),
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=hardProcess,
                    ),
                ),
            )
        )
    if npileup > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=customLogLevel(),
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=pileupProcess,
                    ),
                ),
            )
        )

    # Input
    evGen = acts.examples.EventGenerator(
        level=customLogLevel(),
        generators=generators,
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    s.addWhiteboardAlias("particles", evGen.config.outputParticles)

    if printParticles:
        s.addAlgorithm(
            acts.examples.ParticlesPrinter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            acts.examples.RootParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "pythia8_particles.root"),
            )
        )

    return s


def addParticleSelection(
    s: acts.examples.Sequencer,
    config: ParticleSelectorConfig,
    inputParticles: str,
    outputParticles: str,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    preselectedParticles: ParticleSelectorConfig
        the particle selection configuration
    inputParticles: str
        the identifier for the input particles to be selected
    outputParticles: str
        the identifier for the final selected particle collection
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    s.addAlgorithm(
        acts.examples.ParticleSelector(
            **acts.examples.defaultKWArgs(
                rhoMin=config.rho[0],
                rhoMax=config.rho[1],
                absZMin=config.absZ[0],
                absZMax=config.absZ[1],
                timeMin=config.time[0],
                timeMax=config.time[1],
                phiMin=config.phi[0],
                phiMax=config.phi[1],
                etaMin=config.eta[0],
                etaMax=config.eta[1],
                absEtaMin=config.absEta[0],
                absEtaMax=config.absEta[1],
                ptMin=config.pt[0],
                ptMax=config.pt[1],
                removeCharged=config.removeCharged,
                removeNeutral=config.removeNeutral,
            ),
            level=customLogLevel(),
            inputParticles=inputParticles,
            outputParticles=outputParticles,
        )
    )


def addFatras(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    preSelectParticles: Optional[ParticleSelectorConfig] = ParticleSelectorConfig(),
    postSelectParticles: Optional[ParticleSelectorConfig] = None,
    enableInteractions=False,
    inputParticles: str = "particles_input",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the detector simulation using Fatras

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Fatras steps (returned from addFatras)
    trackingGeometry : tracking geometry
    field : magnetic field
    rnd : RandomNumbers
        random number generator
    preSelectParticles : ParticleSelectorConfig(rho, absZ, time, phi, eta, absEta, pt, removeCharged, removeNeutral), None
        ParticleSelector configuration to select particles as input to Fatras. Each range is specified as a tuple of (min,max).
        Default of no selections specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSelector.hpp
        Specify preSelectParticles=None to inhibit ParticleSelector altogether.
    postSelectParticles : ParticleSelectorConfig(rho, absZ, time, phi, eta, absEta, pt, removeCharged, removeNeutral), None
        Similar to preSelectParticles but applied after simulation to "particles_initial", therefore also filters secondaries.
    enableInteractions : Enable the particle interactions in the simulation
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Selector
    if preSelectParticles is not None:
        particles_selected = "particles_selected"
        addParticleSelection(
            s,
            preSelectParticles,
            inputParticles=inputParticles,
            outputParticles=particles_selected,
        )
    else:
        particles_selected = inputParticles

    # Simulation
    alg = acts.examples.FatrasSimulation(
        level=customLogLevel(),
        inputParticles=particles_selected,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
        emScattering=enableInteractions,
        emEnergyLossIonisation=enableInteractions,
        emEnergyLossRadiation=enableInteractions,
        emPhotonConversion=enableInteractions,
    )

    # Sequencer
    s.addAlgorithm(alg)

    # Selector
    if postSelectParticles is not None:
        particlesInitial = "particles_initial_selected"
        addParticleSelection(
            s,
            postSelectParticles,
            inputParticles=alg.config.outputParticlesInitial,
            outputParticles=particlesInitial,
        )

        particlesFinal = "particles_final_selected"
        addParticleSelection(
            s,
            postSelectParticles,
            inputParticles=alg.config.outputParticlesFinal,
            outputParticles=particlesFinal,
        )
    else:
        particlesInitial = alg.config.outputParticlesInitial
        particlesFinal = alg.config.outputParticlesFinal

    # Only add alias for 'particles_initial' as this is the one we use most
    s.addWhiteboardAlias("particles", particlesInitial)

    # Output
    addSimWriters(
        s,
        alg.config.outputSimHits,
        outputDirCsv,
        outputDirRoot,
        logLevel=logLevel,
        particlesInitial=particlesInitial,
        particlesFinal=particlesFinal,
    )

    return s


def addSimWriters(
    s: acts.examples.Sequencer,
    inputSimHits: Optional[str] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    particlesInitial="particles_initial",
    particlesFinal="particles_final",
) -> None:

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=customLogLevel(),
                outputDir=str(outputDirCsv),
                inputParticles=particlesInitial,
                outputStem="particles_initial",
            )
        )
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=customLogLevel(),
                outputDir=str(outputDirCsv),
                inputParticles=particlesFinal,
                outputStem="particles_final",
            )
        )
        s.addWriter(
            acts.examples.CsvSimHitWriter(
                level=customLogLevel(),
                inputSimHits=inputSimHits,
                outputDir=str(outputDirCsv),
                outputStem="hits",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=customLogLevel(),
                inputParticles="particles_initial",
                filePath=str(outputDirRoot / "particles_initial.root"),
            )
        )
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=customLogLevel(),
                inputParticles="particles_final",
                filePath=str(outputDirRoot / "particles_final.root"),
            )
        )
        s.addWriter(
            acts.examples.RootSimHitWriter(
                level=customLogLevel(),
                inputSimHits=inputSimHits,
                filePath=str(outputDirRoot / "hits.root"),
            )
        )


def getG4DetectorContruction(
    detector: Any,
) -> Any:
    try:
        from acts.examples import TelescopeDetector
        from acts.examples.geant4 import TelescopeG4DetectorConstruction

        if type(detector) is TelescopeDetector:
            return TelescopeG4DetectorConstruction(detector)
    except Exception as e:
        print(e)

    try:
        from acts.examples.dd4hep import DD4hepDetector
        from acts.examples.geant4.dd4hep import DDG4DetectorConstruction

        if type(detector) is DD4hepDetector:
            return DDG4DetectorConstruction(detector)
    except Exception as e:
        print(e)

    raise AttributeError(f"cannot find a suitable detector construction for {detector}")


def addGeant4(
    s: acts.examples.Sequencer,
    detector: Optional[Any],
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    g4detectorConstruction: Optional[Any] = None,
    volumeMappings: List[str] = [],
    materialMappings: List[str] = [],
    inputParticles: str = "particles_input",
    preSelectParticles: Optional[ParticleSelectorConfig] = ParticleSelectorConfig(),
    postSelectParticles: Optional[ParticleSelectorConfig] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the detector simulation using Geant4

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Geant4 steps (returned from addGeant4)
    trackingGeometry : tracking geometry
    field : magnetic field
    rnd : RandomNumbers, None
        random number generator
    preSelectParticles : ParticleSelectorConfig(rho, absZ, time, phi, eta, absEta, pt, removeCharged, removeNeutral), None
        ParticleSelector configuration to select particles as input to Geant4. Each range is specified as a tuple of (min,max).
        Default of no selections specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/ParticleSelector.hpp
        Specify preSelectParticles=None to inhibit ParticleSelector altogether.
    postSelectParticles : ParticleSelectorConfig(rho, absZ, time, phi, eta, absEta, pt, removeCharged, removeNeutral), None
        Similar to preSelectParticles but applied after simulation to "particles_initial", therefore also filters secondaries.
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    """

    from acts.examples.geant4 import Geant4Simulation, makeGeant4SimulationConfig

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Selector
    if preSelectParticles is not None:
        particles_selected = "particles_selected"
        addParticleSelection(
            s,
            preSelectParticles,
            inputParticles=inputParticles,
            outputParticles=particles_selected,
        )
    else:
        particles_selected = inputParticles

    if g4detectorConstruction is None:
        if detector is None:
            raise AttributeError("detector not given")
        g4detectorConstruction = getG4DetectorContruction(detector)

    g4conf = makeGeant4SimulationConfig(
        level=customLogLevel(),
        detector=g4detectorConstruction,
        randomNumbers=rnd,
        inputParticles=particles_selected,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        volumeMappings=volumeMappings,
        materialMappings=materialMappings,
    )
    g4conf.outputSimHits = "simhits"
    g4conf.outputParticlesInitial = "particles_initial"
    g4conf.outputParticlesFinal = "particles_final"

    # Simulation
    alg = Geant4Simulation(
        level=customLogLevel(),
        config=g4conf,
    )

    # Sequencer
    s.addAlgorithm(alg)

    # Selector
    if postSelectParticles is not None:
        particlesInitial = "particles_initial_selected"
        addParticleSelection(
            s,
            postSelectParticles,
            inputParticles=g4conf.outputParticlesInitial,
            outputParticles=particlesInitial,
        )

        particlesFinal = "particles_final_selected"
        addParticleSelection(
            s,
            postSelectParticles,
            inputParticles=g4conf.outputParticlesFinal,
            outputParticles=particlesFinal,
        )
    else:
        particlesInitial = alg.config.outputParticlesInitial
        particlesFinal = alg.config.outputParticlesFinal

    # Only add alias for 'particles_initial' as this is the one we use most
    s.addWhiteboardAlias("particles", particlesInitial)

    # Output
    addSimWriters(
        s,
        alg.config.outputSimHits,
        outputDirCsv,
        outputDirRoot,
        logLevel=logLevel,
        particlesInitial=particlesInitial,
        particlesFinal=particlesFinal,
    )

    return s


def addDigitization(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Union[Path, str],
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    doMerge: Optional[bool] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.Sequencer:
    """This function steers the digitization step

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Digitization steps (returned from addDigitization)
    trackingGeometry : tracking geometry
    field : magnetic field
    digiConfigFile : Path|str, path
        Configuration (.json) file for digitization or smearing description
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    rnd : RandomNumbers, None
        random number generator
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    # Preliminaries
    rnd = rnd or acts.examples.RandomNumbers()

    # Digitization
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(
            str(digiConfigFile),
        ),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits="simhits",
        outputSourceLinks="sourcelinks",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="measurement_particles_map",
        outputMeasurementSimHitsMap="measurement_simhits_map",
        doMerge=doMerge,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, customLogLevel())

    s.addAlgorithm(digiAlg)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        rmwConfig = acts.examples.RootMeasurementWriter.Config(
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputSimHits=digiAlg.config.inputSimHits,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            filePath=str(outputDirRoot / f"{digiAlg.config.outputMeasurements}.root"),
            trackingGeometry=trackingGeometry,
        )
        rmwConfig.addBoundIndicesFromDigiConfig(digiAlg.config)
        s.addWriter(acts.examples.RootMeasurementWriter(rmwConfig, customLogLevel()))

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            acts.examples.CsvMeasurementWriter(
                level=customLogLevel(),
                inputMeasurements=digiAlg.config.outputMeasurements,
                inputClusters=digiAlg.config.outputClusters,
                inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
                outputDir=str(outputDirCsv),
            )
        )

    return s
