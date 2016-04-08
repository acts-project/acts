#include "GaudiKernel/DeclareFactoryEntries.h"

#include "FatrasTools/PhotonConversionSampler.h"
#include "FatrasTools/HadronicInteractionParametricSampler.h"
#include "FatrasTools/EnergyLossSampler.h"
#include "FatrasTools/ElectronEnergyLossSampler.h"
#include "FatrasTools/MultipleScatteringSamplerHighland.h"
#include "FatrasTools/MultipleScatteringSamplerGaussianMixture.h"
#include "FatrasTools/MultipleScatteringSamplerGeneralMixture.h"


using namespace Acts;

DECLARE_TOOL_FACTORY( PhotonConversionSampler )
DECLARE_TOOL_FACTORY( HadronicInteractionParametricSampler )
DECLARE_TOOL_FACTORY( EnergyLossSampler )
DECLARE_TOOL_FACTORY( ElectronEnergyLossSampler )
DECLARE_TOOL_FACTORY( MultipleScatteringSamplerHighland )
DECLARE_TOOL_FACTORY( MultipleScatteringSamplerGaussianMixture )
DECLARE_TOOL_FACTORY( MultipleScatteringSamplerGeneralMixture )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( FatrasTools )
{
  DECLARE_TOOL( PhotonConversionSampler )
  DECLARE_TOOL( HadronicInteractionParametricSampler )
  DECLARE_TOOL( EnergyLossSampler )
  DECLARE_TOOL( ElectronEnergyLossSampler )
  DECLARE_TOOL( MultipleScatteringSamplerHighland )
  DECLARE_TOOL( MultipleScatteringSamplerGaussianMixture )
  DECLARE_TOOL( MultipleScatteringSamplerGeneralMixture )
}

