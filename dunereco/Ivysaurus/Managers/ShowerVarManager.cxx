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
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/Ivysaurus/Managers/ShowerVarManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

namespace ivysaurus
{

ShowerVarManager::ShowerVars::ShowerVars() : 
    m_displacement(-1.f),
    m_DCA(-1.f),
    m_trackStubLength(-1.f),
    m_nuVertexAvSeparation(-1.f),
    m_nuVertexChargeAsymmetry(-1.f),
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
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel"))
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
    FillTrackStub(evt, shower, showerVars);
    FillNuVertexAvSeparation(evt, pfparticle, showerVars);
    FillNuVertexChargeAsymmetry(evt, shower, showerVars);

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

    // DCA
    const double alpha = std::fabs((shower->ShowerStart() - nuVertexPosition).Dot(shower->Direction()));
    const TVector3 r = shower->ShowerStart() - (alpha * shower->Direction());
    const float dca = (r - nuVertexPosition).Mag();

    showerVars.SetDisplacement(displacement);
    showerVars.SetDCA(dca);
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

void ShowerVarManager::FillTrackStub(const art::Event &evt, const art::Ptr<recob::Shower> shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    art::Handle< std::vector<recob::Shower> > showerListHandle;
    evt.getByLabel(m_showerModuleLabel, showerListHandle);

    art::FindManyP<recob::Track> initialTrackAssn(showerListHandle, evt, m_showerModuleLabel);
    std::vector<art::Ptr<recob::Track>> initialTrackStub = initialTrackAssn.at(shower.key());

    if (initialTrackStub.empty())
        return;

    const art::Ptr<recob::Track> &trackStub = initialTrackStub.front();

    const float trackStubLength = std::sqrt((trackStub->Start() - trackStub->End()).Mag2());
    showerVars.SetTrackStubLength(trackStubLength);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillNuVertexAvSeparation(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return;

    const art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    if (!dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);

    const art::Ptr<recob::PFParticle> &nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());
    const TVector3 showerStart = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const TVector3 displacement = (showerStart - nuVertexPosition).Unit();

    // Calc charge weighted axis-hit separation
    float totalCharge = 0.f;
    float numeratorSum = 0.f;

    const std::vector<art::Ptr<recob::Hit>> collectionHits = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfparticle, evt, m_recoModuleLabel, 2);

    for (const art::Ptr<recob::Hit> &hit : collectionHits)
    {
        const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacepoints.empty())
            continue;

        const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.front();

        const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z()) - nuVertexPosition;
        const float transverse = std::sqrt(displacement.Cross(spacepointPos).Mag2());
        const double charge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});

        totalCharge += charge;
        numeratorSum += (transverse * charge);
    }

    const float nuVertexAvSeparation = (totalCharge < std::numeric_limits<float>::epsilon()) ? 0.f : (numeratorSum / totalCharge);
    showerVars.SetNuVertexAvSeparation(nuVertexAvSeparation);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::FillNuVertexChargeAsymmetry(const art::Event &evt, const art::Ptr<recob::Shower> &shower, 
    ShowerVarManager::ShowerVars &showerVars) const
{
    if (!dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel))
        return;

    const art::Ptr<recob::PFParticle> &nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, m_recoModuleLabel);
    const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, m_recoModuleLabel);
    const TVector3 nuVertexPosition = TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z());

    // Get initial direction of the shower from the track stub
    const TVector3 showerDirection = TVector3(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());

    const std::vector<art::Ptr<recob::Hit>> &hits = dune_ana::DUNEAnaShowerUtils::GetHits(shower, evt, m_showerModuleLabel);

    std::vector<art::Ptr<recob::Hit>> hitsU, hitsV, hitsW;

    for (const art::Ptr<recob::Hit> &hit : hits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;

        const IvysaurusUtils::PandoraView pandoraView = IvysaurusUtils::GetPandora2DView(hit);

        std::vector<art::Ptr<recob::Hit>> &hitVector = pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_U ? hitsU :
            pandoraView == IvysaurusUtils::PandoraView::TPC_VIEW_V ? hitsV : hitsW;

        hitVector.push_back(hit);
    }

    const float nuVertexChargeAsymmetryU = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsU, IvysaurusUtils::PandoraView::TPC_VIEW_U);
    const float nuVertexChargeAsymmetryV = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsV, IvysaurusUtils::PandoraView::TPC_VIEW_V);
    const float nuVertexChargeAsymmetryW = GetViewNuVertexChargeAsymmetry(evt, nuVertexPosition, showerDirection, hitsW, IvysaurusUtils::PandoraView::TPC_VIEW_W);
    const float maxNuVertexChargeAsymmetry = std::max(std::max(nuVertexChargeAsymmetryU, nuVertexChargeAsymmetryV), nuVertexChargeAsymmetryW);

    showerVars.SetNuVertexChargeAsymmetry(maxNuVertexChargeAsymmetry);
}

/////////////////////////////////////////////////////////////

float ShowerVarManager::GetViewNuVertexChargeAsymmetry(const art::Event &evt, const TVector3 &nuVertexPosition, const TVector3 &showerDirection, 
    const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const TVector3 viewNuVertex = IvysaurusUtils::ProjectIntoPandoraView(nuVertexPosition, pandoraView);
    const TVector3 point2 = nuVertexPosition + (showerDirection * 10.f);
    const TVector3 viewPoint2 = IvysaurusUtils::ProjectIntoPandoraView(point2, pandoraView);
    const TVector3 centralAxis = (viewPoint2 - viewNuVertex).Unit();
    const TVector3 yAxis = TVector3(0.f, 1.f, 0.f);
    const TVector3 orthAxis = centralAxis.Cross(yAxis).Unit();

    float chargeAsymmetry = 0.f;
    float totalCharge = 0.f;

    for (const art::Ptr<recob::Hit> &viewHit : viewHits)
    {
        const TVector3 hitPosition = IvysaurusUtils::ObtainPandoraHitPosition(evt, viewHit, pandoraView);
        const float l = orthAxis.Dot(hitPosition - viewNuVertex);
        float charge = std::fabs(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {viewHit}));

        totalCharge += charge;

        charge *= (l < 0.f) ? -1.0 : 1.0;
        chargeAsymmetry += charge;
    }

    chargeAsymmetry = totalCharge < std::numeric_limits<float>::epsilon() ? 0.f : (chargeAsymmetry / totalCharge);

    return std::fabs(chargeAsymmetry);
}

/////////////////////////////////////////////////////////////

void ShowerVarManager::Reset(ShowerVarManager::ShowerVars &showerVars) const
{
    showerVars.SetDisplacement(-1.f);
    showerVars.SetDCA(-1.f);
    showerVars.SetTrackStubLength(-1.f);
    showerVars.SetNuVertexAvSeparation(-1.f);
    showerVars.SetNuVertexChargeAsymmetry(-1.f);
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
