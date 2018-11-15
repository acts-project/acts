
// THIS IS JUST A TEST TO SEE IF RESULT DIFFERENCES BETWEEN ACTS AND ATLAS ARE DUE TO HARDCODED BINFINDING
template <typename SpacePoint>
std::set<size_t>
Acts::ATLASTopBinFinder<SpacePoint>::findBins(size_t phiBin, size_t zBin,const Acts::SpacePointGrid<SpacePoint>* binnedSP){

  std::set<size_t> neighbourBins = binnedSP->neighborHoodIndices({phiBin,zBin},1);  
  if(zBin == 5){
    return neighbourBins;
  }
// <10 not necessary because grid doesn't return bins that don't exist (only works for 11 zbins)
  if(zBin >5){
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin,zBin-1  }));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin-1,zBin-1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin+1,zBin-1}));
  }else{
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin,zBin+1  }));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin-1,zBin+1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin+1,zBin+1}));
  }
  return neighbourBins;
}

