////////////////////////////////////////////////////////////////////////
/// \file    TrueEnergyCalc.cxx
/// \brief   A class to calculate true energies from sim edeps
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/FDSelections/EnergyReco/TrueEnergyCalc.h"

#include "nusimdata/SimulationBase/MCTruth.h"

namespace FDSelection
{

/////////////////////////////////////////////////////////////

TrueEnergyCalc::TrueEnergyCalc(const fhicl::ParameterSet& pset) :
    m_simEDepLabel(pset.get<std::string>("SimEDepModuleLabel")),
    m_generatorLabel(pset.get<std::string>("GeneratorModuleLabel")),
    m_geantLabel(pset.get<std::string>("GeantModuleLabel"))
{
}

/////////////////////////////////////////////////////////////

void TrueEnergyCalc::InitialiseCalc(const art::Event &evt)
{
    Reset();

    // Find neutrino
    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;

    if (evt.getByLabel(m_generatorLabel, mcTruthHandle))
        art::fill_ptr_vector(mcTruthVector, mcTruthHandle);

    if (mcTruthVector.empty() || (mcTruthVector.at(0)->Origin() != simb::kBeamNeutrino))
        return;

    // Find neutrino children
    std::vector<art::Ptr<simb::MCParticle>> mcNuChildren;

    art::FindManyP<simb::MCParticle> mcParticleAssoc = art::FindManyP<simb::MCParticle>(mcTruthHandle, evt, m_geantLabel);
    const std::vector<art::Ptr<simb::MCParticle>> mcParticleVector = mcParticleAssoc.at(mcTruthVector.at(0).key());
    lar_pandora::LArPandoraHelper::BuildMCParticleMap(mcParticleVector, m_mcParticleMap);

    // Finally fill out MCParticle -> SimEDeps
    art::Handle<std::vector<sim::SimEnergyDeposit>> simEDepHandle;
    std::vector<art::Ptr<sim::SimEnergyDeposit>> simEDepVector;

    if (evt.getByLabel(m_simEDepLabel, simEDepHandle))
        art::fill_ptr_vector(simEDepVector, simEDepHandle);

    for (const art::Ptr<sim::SimEnergyDeposit> &simEDep : simEDepVector)
    {
        const int ownerTrackID = std::abs(simEDep->TrackID());

        if (m_mcParticleMap.find(ownerTrackID) == m_mcParticleMap.end())
            continue;

        m_simEDepMap[ownerTrackID].push_back(simEDep);
    }
}

/////////////////////////////////////////////////////////////

double TrueEnergyCalc::GetTrueNuEnergy()
{
    double totalEnergy = 0.0;

    for (const auto &entry : m_simEDepMap)
    {
        for (const art::Ptr<sim::SimEnergyDeposit> &simEDep : entry.second)
        {
            totalEnergy += simEDep->Energy();
        }
    }

    return totalEnergy;
}

/////////////////////////////////////////////////////////////

double TrueEnergyCalc::GetTrueParticleEnergy(const int trackID)
{
    double totalEnergy = 0.0;

    if (m_simEDepMap.find(trackID) == m_simEDepMap.end())
        return -999.9;

    const std::vector<art::Ptr<sim::SimEnergyDeposit>> &simEDeps = m_simEDepMap.at(trackID);

    for (const art::Ptr<sim::SimEnergyDeposit> &simEDep : simEDeps)
        totalEnergy += simEDep->Energy();

    return totalEnergy;
}

/////////////////////////////////////////////////////////////

void TrueEnergyCalc::Reset()
{
    // Clear our algorithm's maps
    m_mcParticleMap.clear();
    m_simEDepMap.clear();
}

/////////////////////////////////////////////////////////////

}
