////////////////////////////////////////////////////////////////////////
/// \file    ShowerVarManager.cxx
/// \brief   A class to manage the Ivysaurus 2D shower variable input 
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

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/Ivysaurus/Managers/ShowerVarManager.h"


namespace ivysaurus
{

ShowerVarManager::ShowerVars::ShowerVars() : 
    m_displacement(-1.f),
    m_initialGapSize(-1.f),
    m_largestGapSize(-1.f),
    m_pathwayLength(-1.f),
    m_pathwayScatteringAngle2D(-1.f),
    m_nShowerHits(-1.f),
    m_foundHitRatio(-1.f),
    m_scatterAngle(-1.f),
    m_openingAngle(-1.f),
    m_nuVertexEnergyAsymmetry(-1.f),
    m_nuVertexEnergyWeightedMeanRadialDistance(-1.f),
    m_showerStartEnergyAsymmetry(-1.f),
    m_showerStartMoliereRadius(-1.f),
    m_nAmbiguousViews(-1.f),
    m_unaccountedEnergy(-999.f)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ShowerVarManager::ShowerVarManager(const fhicl::ParameterSet& pset) :
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel"))
{
}

/////////////////////////////////////////////////////////////

ShowerVarManager::~ShowerVarManager()
{
}

/////////////////////////////////////////////////////////////

bool ShowerVarManager::EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    FillDisplacement(evt, shower, showerVars);
    FillConnectionPathwayVars(evt, pfparticle, showerVars);

    return true;
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillDisplacement(const art::Event &evt, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    const art::Ptr<recob::PFParticle> &nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());
    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float displacement = (nuVertexPosition - showerStart).Mag();

    showerVars.SetDisplacement(displacement);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillConnectionPathwayVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const larpandoraobj::PFParticleMetadata::PropertiesMap &propertyMap = metadata->GetPropertiesMap();

    if (propertyMap.find("MaxInitialGapSize") != propertyMap.end()) showerVars.SetInitialGapSize(propertyMap.at("MaxInitialGapSize"));
    if (propertyMap.find("MinLargestProjectedGapSize") != propertyMap.end()) showerVars.SetLargestGapSize(propertyMap.at("MinLargestProjectedGapSize"));
    if (propertyMap.find("PathwayLengthMin") != propertyMap.end()) showerVars.SetPathwayLength(propertyMap.at("PathwayLengthMin"));
    if (propertyMap.find("MaxShowerStartPathwayScatteringAngle2D") != propertyMap.end()) showerVars.SetPathwayScatteringAngle2D(propertyMap.at("MaxShowerStartPathwayScatteringAngle2D"));
    if (propertyMap.find("MaxNPostShowerStartHits") != propertyMap.end()) showerVars.SetNShowerHits(propertyMap.at("MaxNPostShowerStartHits"));
    if (propertyMap.find("MaxFoundHitRatio") != propertyMap.end()) showerVars.SetFoundHitRatio(propertyMap.at("MaxFoundHitRatio"));
    if (propertyMap.find("MaxPostShowerStartScatterAngle") != propertyMap.end()) showerVars.SetScatterAngle(propertyMap.at("MaxPostShowerStartScatterAngle"));
    if (propertyMap.find("MaxPostShowerStartOpeningAngle") != propertyMap.end()) showerVars.SetOpeningAngle(propertyMap.at("MaxPostShowerStartOpeningAngle"));
    if (propertyMap.find("MaxPostShowerStartNuVertexEnergyAsymmetry") != propertyMap.end()) showerVars.SetNuVertexEnergyAsymmetry(propertyMap.at("MaxPostShowerStartNuVertexEnergyAsymmetry"));
    if (propertyMap.find("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance") != propertyMap.end()) 
        showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(propertyMap.at("MaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance"));
    if (propertyMap.find("MaxPostShowerStartShowerStartEnergyAsymmetry") != propertyMap.end()) showerVars.SetShowerStartEnergyAsymmetry(propertyMap.at("MaxPostShowerStartShowerStartEnergyAsymmetry"));
    if (propertyMap.find("MinPostShowerStartShowerStartMoliereRadius") != propertyMap.end()) showerVars.SetShowerStartMoliereRadius(propertyMap.at("MinPostShowerStartShowerStartMoliereRadius"));
    if (propertyMap.find("NViewsWithAmbiguousHits") != propertyMap.end()) showerVars.SetNAmbiguousViews(propertyMap.at("NViewsWithAmbiguousHits"));
    if (propertyMap.find("AmbiguousHitMaxUnaccountedEnergy") != propertyMap.end()) showerVars.SetUnaccountedEnergy(propertyMap.at("AmbiguousHitMaxUnaccountedEnergy"));
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::Reset(ShowerVarManager::ShowerVars &showerVars) const
{
    showerVars.SetDisplacement(-1.f);
    showerVars.SetInitialGapSize(-1.f);
    showerVars.SetLargestGapSize(-1.f);
    showerVars.SetPathwayLength(-1.f);
    showerVars.SetPathwayScatteringAngle2D(-1.f);
    showerVars.SetNShowerHits(-1.f);
    showerVars.SetFoundHitRatio(-1.f);
    showerVars.SetScatterAngle(-1.f);
    showerVars.SetOpeningAngle(-1.f);
    showerVars.SetNuVertexEnergyAsymmetry(-1.f);
    showerVars.SetNuVertexEnergyWeightedMeanRadialDistance(-1.f);
    showerVars.SetShowerStartEnergyAsymmetry(-1.f);
    showerVars.SetShowerStartMoliereRadius(-1.f);
    showerVars.SetNAmbiguousViews(-1.f);
    showerVars.SetUnaccountedEnergy(-999.f);
}

/////////////////////////////////////////////////////////////

}
