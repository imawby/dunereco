////////////////////////////////////////////////////////////////////////
/// \file    TrueEnergyCalc.h
/// \brief   A class to calculate true energies from sim edeps
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef TRUEENERGYCALC_H
#define TRUEENERGYCALC_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace FDSelection
{

class TrueEnergyCalc
{
  public:
    TrueEnergyCalc(const fhicl::ParameterSet& pset);

    void InitialiseCalc(const art::Event &evt);
    double GetTrueNuEnergy();
    double GetTrueParticleEnergy(const int trackID);
    void Reset();

  private:
    lar_pandora::MCParticleMap m_mcParticleMap;
    std::map<int, std::vector<art::Ptr<sim::SimEnergyDeposit>>> m_simEDepMap;

    std::string m_simEDepLabel;
    std::string m_generatorLabel;
    std::string m_geantLabel;
};

}

#endif  // TRUEENERGYCALC
