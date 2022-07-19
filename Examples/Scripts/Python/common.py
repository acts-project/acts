# import acts
# import acts.examples
# import os
# import sys


# from pathlib import Path


# def addPythia8(
#     sequencer: acts.examples.Sequencer,
#     rnd: acts.examples.RandomNumbers,
#     nhard=1,
#     npileup=200,
#     beam0=acts.PdgParticle.eProton,
#     beam1=acts.PdgParticle.eProton,
#     cmsEnergy=14 * acts.UnitConstants.TeV,
#     vertexStddev: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
#     vertexMean: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
# ):
#     """This function steers the particle generation using Pythia8

#     NB. this version is included here only for compatibility. Please use pythia8.addPythia8 instead.
#     """
#     import acts.examples.sumulation

#     evGen = acts.examples.simulation.addPythia8(
#         sequencer,
#         rnd=rnd,
#         nhard=nhard,
#         npileup=npileup,
#         beam=(beam0, beam1),
#         cmsEnergy=cmsEnergy,
#         vtxGen=acts.examples.GaussianVertexGenerator(
#             stddev=vertexStddev, mean=vertexMean
#         ),
#         returnEvGen=True,
#     )

#     return evGen
