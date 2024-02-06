////////////////////////////////////////////////////////////////////////
/// \file    GridManager.cxx
/// \brief   A class to manage the Ivysaurus 2D grid input 
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

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "dunereco/Ivysaurus/Managers/GridManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaSpacePointUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"


namespace ivysaurus
{

GridManager::Grid::Grid(const TVector3 origin, const float driftSpan, const float wireSpan, 
    const unsigned int dimensions, const float maxGridEntry, const unsigned int nSigmaConsidered, const float integralStep,
    const IvysaurusUtils::PandoraView pandoraView, const bool isInitialised) : 
        m_axisDimensions(dimensions),
        m_maxGridEntry(maxGridEntry),
        m_nSigmaConsidered(nSigmaConsidered),
        m_integralStep(integralStep),
        m_pandoraView(pandoraView),
        m_isInitialised(isInitialised)
{
    m_gridValues = std::vector<std::vector<float>>(m_axisDimensions, std::vector<float>(m_axisDimensions, 0.0));

    const float driftInterval = driftSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_driftBoundaries.push_back(origin.X() + (i * driftInterval));

    const float wireInterval = wireSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_wireBoundaries.push_back(origin.Z() + (i * wireInterval));

    m_isNormalised = false;
}

/////////////////////////////////////////////////////////////

bool GridManager::Grid::IsInsideGrid(const TVector3 &position, const float width) const
{
    //////////////////////////////////////
    // Wire axis
    const float gridMinWireCoord(std::min(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float gridMaxWireCoord(std::max(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float hitWireCoord(position.Z());

    if (std::fabs(gridMinWireCoord - hitWireCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if (std::fabs(gridMaxWireCoord - hitWireCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if ((hitWireCoord < gridMinWireCoord) || (hitWireCoord > gridMaxWireCoord))
        return false;
    //////////////////////////////////////

    //////////////////////////////////////
    // Drift axis
    const float gridMinDriftCoord(std::min(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float gridMaxDriftCoord(std::max(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float hitMinDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() - (m_nSigmaConsidered * (width / 2.0));
    const float hitMaxDriftCoord = (width < std::numeric_limits<float>::epsilon()) ? position.X() : position.X() + (m_nSigmaConsidered * (width / 2.0));

    if (std::fabs(hitMaxDriftCoord - gridMinDriftCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if (std::fabs(hitMinDriftCoord - gridMaxDriftCoord) < std::numeric_limits<float>::epsilon())
        return false;

    if ((hitMaxDriftCoord < gridMinDriftCoord) || (hitMinDriftCoord > gridMaxDriftCoord))
        return false;
    //////////////////////////////////////

    return true;
}

/////////////////////////////////////////////////////////////

void GridManager::Grid::AddToGrid(const TVector3 &position, const float width, const float energy)
{
    //////////////////////////////////////
    // Get wire bin
    const float wireInterval = std::fabs(m_wireBoundaries.at(0) - m_wireBoundaries.at(1));
    int wireBin = std::floor((position.Z() - m_wireBoundaries.front()) / wireInterval); 

    if (m_wireBoundaries.back() < m_wireBoundaries.front())
        wireBin = std::floor((m_wireBoundaries.front() - position.Z()) / wireInterval);

    if (wireBin < 0)
        return;

    if (wireBin >= static_cast<int>(m_axisDimensions))
        return;
    //////////////////////////////////////

    //////////////////////////////////////
    // Now fill assuming hits are Gaussian...
    const float driftInterval = std::fabs(m_driftBoundaries.at(0) - m_driftBoundaries.at(1));
    const float hitLowEdge = position.X() - (m_nSigmaConsidered * (width / 2.f));
    const float hitHighEdge = position.X() + (m_nSigmaConsidered * (width / 2.f));

    // Link up 'low X' and 'high X' with grid's 'start' and 'end' definitions
    float hitStartEdge = hitLowEdge;
    float hitEndEdge = hitHighEdge;
    int startDriftBin = std::floor((hitStartEdge - m_driftBoundaries.front()) / driftInterval); 
    int endDriftBin = std::floor((hitEndEdge - m_driftBoundaries.front()) / driftInterval);

    if (m_driftBoundaries.back() < m_driftBoundaries.front())
    {
        hitStartEdge = hitHighEdge;
        hitEndEdge = hitLowEdge;
        startDriftBin = std::floor((m_driftBoundaries.front() - hitStartEdge) / driftInterval);
        endDriftBin = std::floor((m_driftBoundaries.front() - hitEndEdge) / driftInterval);
    }

    // Loop over the drift bings, and fill grid
    for (int iDriftBin = startDriftBin; iDriftBin <= endDriftBin; ++iDriftBin)
    {
        if (iDriftBin < 0)
            continue;

        if (iDriftBin >= static_cast<int>(m_axisDimensions))
            continue;

        const float integralStartX = (iDriftBin == startDriftBin) ? hitStartEdge : m_driftBoundaries.at(iDriftBin);
        const float integralEndX = (iDriftBin == endDriftBin) ? hitEndEdge : m_driftBoundaries.at(iDriftBin + 1);

        // Integrate the Gaussian area...
        const float chargeFraction = IvysaurusUtils::IntegrateGaussian(integralStartX, integralEndX, position.X(), (width / 2.f), m_integralStep); 
        const float entryEnergy = energy * chargeFraction;

        // Now fill grid
        m_gridValues[iDriftBin][wireBin] += entryEnergy;
    }
    //////////////////////////////////////
}

/////////////////////////////////////////////////////////////

void GridManager::Grid::NormaliseGrid()
{
    if (m_isNormalised)
        throw cet::exception("ivysaur::GridManager") << "the entries are already normalised!";

    for (unsigned int driftIndex = 0; driftIndex < m_axisDimensions; ++driftIndex)
    {
        for (unsigned int wireIndex = 0; wireIndex < m_axisDimensions; ++wireIndex)
        {
            float gridEntry = m_gridValues[driftIndex][wireIndex];

            if (gridEntry > m_maxGridEntry)
                gridEntry = m_maxGridEntry;

            gridEntry /= m_maxGridEntry;

            m_gridValues[driftIndex][wireIndex] = gridEntry;
        }
    }

    m_isNormalised = true;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GridManager::GridManager(const fhicl::ParameterSet& pset) :
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_gridSize3D(pset.get<float>("GridSize3D")),
    m_dimensions(pset.get<float>("GridDimensions")),
    m_addChildrenToGrid(pset.get<bool>("AddChildrenToGrid")),
    m_maxGridEntry(pset.get<float>("MaxGridEntry")),
    m_nSigmaConsidered(pset.get<unsigned int>("NSigmaConsidered")),
    m_integralStep(pset.get<float>("IntegralStep")),
    m_recombFactor(pset.get<float>("RecombFactor")),
    m_calorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
}

/////////////////////////////////////////////////////////////

GridManager::~GridManager()
{
}

/////////////////////////////////////////////////////////////

GridManager::Grid GridManager::ObtainViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    const IvysaurusUtils::PandoraView pandoraView, const bool isStart) const
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return Grid(TVector3(0.f, 0.f, 0.f), 0.f, 0.f, 0, m_maxGridEntry, m_nSigmaConsidered, m_integralStep, pandoraView, false);

    const float trackScore = metaMap.at("TrackScore");

    // Find the extremal diagonal.. 
    TVector3 position1 = TVector3(0.f, 0.f, 0.f);
    TVector3 position2 = TVector3(0.f, 0.f, 0.f); // Along the particle direction from position1

    if (isStart)
    {
        if (!GetStartExtremalPoints(evt, pfparticle, position1, position2))
            return Grid(TVector3(0.f, 0.f, 0.f), 0.f, 0.f, 0, m_maxGridEntry, m_nSigmaConsidered, m_integralStep, pandoraView, false);
    }
    else
    {
        if (trackScore > 0.5f)
        {
            if (!GetEndExtremalPointsTrack(evt, pfparticle, position1, position2))
            {
                if (!GetEndExtremalPointsShower(evt, pfparticle, position1, position2))
                {
                    return Grid(TVector3(0.f, 0.f, 0.f), 0.f, 0.f, 0, m_maxGridEntry, m_nSigmaConsidered, m_integralStep, pandoraView, false);
                }
            }
        }
        else
        {
            if (!GetEndExtremalPointsShower(evt, pfparticle, position1, position2))
            {
                if (!GetEndExtremalPointsTrack(evt, pfparticle, position1, position2))
                {
                    return Grid(TVector3(0.f, 0.f, 0.f), 0.f, 0.f, 0, m_maxGridEntry, m_nSigmaConsidered, m_integralStep, pandoraView, false);
                }
            }
        }
    }

    // Now need to project these things into the 'Pandora view'
    const TVector3 projectedPosition1 = ProjectIntoPandoraView(position1, pandoraView);
    const TVector3 projectedPosition2 = ProjectIntoPandoraView(position2, pandoraView);
    const float driftSpan = projectedPosition2.X() - projectedPosition1.X();
    const float wireSpan = projectedPosition2.Z() - projectedPosition1.Z();

    Grid grid = Grid(projectedPosition1, driftSpan, wireSpan, m_dimensions, m_maxGridEntry, m_nSigmaConsidered, m_integralStep, pandoraView, true);

    FindHitsInGrid(evt, pfparticle, grid);

    return grid;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetStartExtremalPoints(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel); 
    const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);

    int nBins = 180;
    float angleMin = 0.f, angleMax = 2.f * M_PI;
    float binWidth = (angleMax - angleMin) / static_cast<float>(nBins);

    std::vector<std::vector<int>> spatialDist(nBins, std::vector<int>(nBins, 0));
    std::vector<std::vector<float>> energyDist(nBins, std::vector<float>(nBins, 0.f));

    // theta0YZ then theta0XZ
    // measure from Y to Z, and Z to X? fool.
    int highestSP = 0;
    float highestEnergy = 0.f;
    int bestTheta0YZBin = -1;
    int bestTheta0XZBin = -1; 

    for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
    {
        const TVector3 vertexPos = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
        const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z());
        const TVector3 displacement = spacepointPos - vertexPos;
        const float mag = sqrt((displacement.X() * displacement.X()) + (displacement.Y() * displacement.Y()) + (displacement.Z() * displacement.Z()));

        if (mag > m_gridSize3D)
            continue;

        const float magXZ = sqrt((displacement.X() * displacement.X()) + (displacement.Z() * displacement.Z()));

        float theta0YZ = (mag < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.Y() / mag) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f : 
            std::acos(displacement.Y() / mag);

        float theta0XZ = (magXZ < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.Z() / magXZ) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f :
            std::acos(displacement.Z() / magXZ);

        // try do signed-ness
        if (displacement.Z() < 0.f)
            theta0YZ += M_PI;

        if (displacement.X() < 0.f)
            theta0XZ += M_PI;

        const int bin0YZ = std::floor(theta0YZ / binWidth);
        const int bin0XZ = std::floor(theta0XZ / binWidth);

        const std::vector<art::Ptr<recob::Hit>> assocHits = dune_ana::DUNEAnaSpacePointUtils::GetHits(spacepoint, evt, m_recoModuleLabel);

        if (assocHits.empty())
            continue;

        spatialDist[bin0YZ][bin0XZ] += 1;
        energyDist[bin0YZ][bin0XZ] += ObtainHitEnergy(evt, assocHits.front());

        if (((spatialDist[bin0YZ][bin0XZ] == highestSP) && (energyDist[bin0YZ][bin0XZ] > highestEnergy)) ||
            (spatialDist[bin0YZ][bin0XZ] > highestSP))
        {
            highestSP = spatialDist[bin0YZ][bin0XZ];
            highestEnergy = energyDist[bin0YZ][bin0XZ];
            bestTheta0YZBin = bin0YZ;
            bestTheta0XZBin = bin0XZ;
        }
    }

    if ((bestTheta0YZBin < 0) || (bestTheta0XZBin < 0))
        return false;

    const float bestTheta0YZ = angleMin + ((static_cast<float>(bestTheta0YZBin) + 0.5f) * binWidth);
    const float bestTheta0XZ = angleMin + ((static_cast<float>(bestTheta0XZBin) + 0.5f) * binWidth);

    /*
    const float y = std::cos(bestTheta0YZ);
    const float xzMag = std::sin(bestTheta0YZ); 
    const float x = xzMag * std::sin(bestTheta0XZ);
    const float z = xzMag * std::cos(bestTheta0XZ);
    */

    TVector3 direction = TVector3(std::fabs(std::sin(bestTheta0YZ) * std::sin(bestTheta0XZ)), std::fabs(std::cos(bestTheta0YZ)), 
        std::fabs(std::sin(bestTheta0YZ) * std::cos(bestTheta0XZ)));

    if (bestTheta0XZ > M_PI)
        direction.SetX(direction.X() * -1.f);

    if (bestTheta0YZ > M_PI)
        direction.SetZ(direction.Z() * -1.f);

    if ((bestTheta0YZ > (M_PI / 2.f)) && (bestTheta0YZ < (M_PI * 3.f / 2.f)))
        direction.SetY(direction.Y() * -1.f);

    position1 = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position2 = position1 + (direction * diagonalLength);

    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetStartExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    const TVector3 direction = TVector3(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
    position1 = TVector3(track->Start().X(), track->Start().Y(), track->Start().Z());
    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position2 = position1 + (direction * diagonalLength);

    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetStartExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    // Get the track stub
    art::Handle<std::vector<recob::Shower>> showerHandle;
    evt.getByLabel(m_showerModuleLabel, showerHandle);

    art::FindManyP<recob::Track> initialTrackAssoc(showerHandle, evt, m_showerModuleLabel);
    std::vector<art::Ptr<recob::Track>> initialTrackStubVector = initialTrackAssoc.at(shower.key());

    if (initialTrackStubVector.size() != 1)
        return false;

    art::Ptr<recob::Track> initialTrackStub = initialTrackStubVector.at(0);

    const TVector3 direction = TVector3(initialTrackStub->StartDirection().X(), initialTrackStub->StartDirection().Y(), initialTrackStub->StartDirection().Z());
    position1 = TVector3(initialTrackStub->Start().X(), initialTrackStub->Start().Y(), initialTrackStub->Start().Z());
    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position2 = position1 + (direction * diagonalLength);

    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetEndExtremalPointsTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    // Is track...
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    const TVector3 direction = TVector3(track->EndDirection().X(), track->EndDirection().Y(), track->EndDirection().Z());
    const TVector3 end = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position1 = end - (direction * (diagonalLength / 2.0)); 
    position2 = end + (direction * (diagonalLength / 2.0));

    return true;
}

/////////////////////////////////////////////////////////////

bool GridManager::GetEndExtremalPointsShower(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    TVector3 &position1, TVector3 &position2) const
{
    // Is shower...
    if (!dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
        return false;

    const art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

    const TVector3 direction = TVector3(shower->Direction().X(), shower->Direction().Y(), shower->Direction().Z());
    const TVector3 start = TVector3(shower->ShowerStart().X(), shower->ShowerStart().Y(), shower->ShowerStart().Z());
    const float length = shower->Length();
    const TVector3 end = start + (length * direction);

    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    position1 = end - (direction * (diagonalLength / 2.0)); 
    position2 = end + (direction * (diagonalLength / 2.0));

    return true;
}

/////////////////////////////////////////////////////////////

void GridManager::FindHitsInGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    Grid &grid) const
{
    IvysaurusUtils::PandoraView pandoraView = grid.GetPandoraView();

    // Get all 2D hits
    std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);

    // Add in children hits...
    if (m_addChildrenToGrid)
    {
        const std::vector<art::Ptr<recob::PFParticle>> pfpChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfparticle, evt, m_recoModuleLabel);

        for (const art::Ptr<recob::PFParticle> &childPFP : pfpChildren)
        {
            const art::Ptr<recob::Vertex> &childVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(childPFP, evt, m_recoModuleLabel);
            const TVector3 vertex3D = TVector3(childVertex->position().X(), childVertex->position().Y(), childVertex->position().Z());
            const TVector3 projectedVertex = IvysaurusUtils::ProjectIntoPandoraView(vertex3D, pandoraView);

            if (!grid.IsInsideGrid(projectedVertex, 0.f))
                continue;

            const std::vector<art::Ptr<recob::Hit>> &childHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(childPFP, evt, m_recoModuleLabel);
            pfpHits.insert(pfpHits.begin(), childHits.begin(), childHits.end());
        }
    }

    // Add hits to hit list
    for (const art::Ptr<recob::Hit> hit : pfpHits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;
                     
        const IvysaurusUtils::PandoraView thisPandoraView = IvysaurusUtils::GetPandora2DView(hit);

        if (thisPandoraView != pandoraView)
            continue;

        // Get 2D hit position and width
        float hitWidth = 0.f;
        TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, pandoraView, pandoraHitPosition, hitWidth);

        // Check hit is inside grid
        if (!grid.IsInsideGrid(pandoraHitPosition, hitWidth))
            continue;

        grid.AddToGridHitList(hit);
    }
}

/////////////////////////////////////////////////////////////

void GridManager::FillViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    GridManager::Grid &grid) const
{
    IvysaurusUtils::PandoraView pandoraView = grid.GetPandoraView();

    // Get hits in grid
    std::vector<art::Ptr<recob::Hit>> gridHitList = grid.GetGridHitList();

    // Fill grid
    for (const art::Ptr<recob::Hit> hit : gridHitList)
    {
        // Get 2D hit position and width
        float hitWidth = 0.f;
        TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, pandoraView, pandoraHitPosition, hitWidth);

        // Check hit is inside grid
        if (!grid.IsInsideGrid(pandoraHitPosition, hitWidth))
            continue;

        // Add its energy to the grid
        const float energy = ObtainHitEnergy(evt, hit);
        grid.AddToGrid(pandoraHitPosition, hitWidth, energy);
    }

    grid.NormaliseGrid();
}

/////////////////////////////////////////////////////////////

float GridManager::ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const double charge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});
    const double nElectrons = m_calorimetryAlg.ElectronsFromADCArea(charge, hit->WireID().Plane);
    const double hitEnergy = nElectrons / m_recombFactor / util::kGeVToElectrons;

    return hitEnergy;
}

/////////////////////////////////////////////////////////////

GridManager::Grid GridManager::ObtainViewDisplacementGrid(const art::Event &evt, const TVector3 &nuVertex3D, const GridManager::Grid &caloGrid) const
{
    // Now fill the displacement grid
    Grid dispGrid = caloGrid;
    dispGrid.ResetGridEntries();

    IvysaurusUtils::PandoraView pandoraView = dispGrid.GetPandoraView();
    const TVector3 nuVertex2D = IvysaurusUtils::ProjectIntoPandoraView(nuVertex3D, pandoraView);
    const std::vector<float> &driftBoundaries = dispGrid.GetDriftBoundaries();
    const std::vector<float> &wireBoundaries = dispGrid.GetWireBoundaries();
    const float driftIntervalSigned = driftBoundaries.at(1) - driftBoundaries.at(0);
    const float wireIntervalSigned = wireBoundaries.at(1) - wireBoundaries.at(0);

    // Get hits in grid
    std::vector<art::Ptr<recob::Hit>> gridHitList = dispGrid.GetGridHitList();

    // Fill grid
    for (const art::Ptr<recob::Hit> hit : gridHitList)
    {
        // Get 2D hit position and width
        float hitWidth = 0.f;
        TVector3 pandoraHitPosition = TVector3(0.f, 0.f, 0.f);
        IvysaurusUtils::ObtainPandoraHitPositionAndWidth(evt, hit, pandoraView, pandoraHitPosition, hitWidth);

        // Check hit is inside grid
        if (!dispGrid.IsInsideGrid(pandoraHitPosition, 0.f))
            continue;

        // Get wire bin
        int wireBin = std::floor((pandoraHitPosition.Z() - wireBoundaries.front()) / std::fabs(wireIntervalSigned));

        if (wireBoundaries.back() < wireBoundaries.front())
            wireBin = std::floor((wireBoundaries.front() - pandoraHitPosition.Z()) / std::fabs(wireIntervalSigned));

        // floating-point precision
        if ((wireBin < 0) || (wireBin >= static_cast<int>(m_dimensions)))
            continue;

        // Get drift bin
        int driftBin = std::floor((pandoraHitPosition.X() - driftBoundaries.front()) / std::fabs(driftIntervalSigned));

        if (driftBoundaries.back() < driftBoundaries.front())
            driftBin = std::floor((driftBoundaries.front() - pandoraHitPosition.X()) / std::fabs(driftIntervalSigned));

        // floating-point precision
        if ((driftBin < 0) || (driftBin >= static_cast<int>(m_dimensions)))
            continue;

        // Has that entry been set?
        if (dispGrid.GetGridEntry(driftBin, wireBin) > std::numeric_limits<float>::epsilon())
            continue;

        // Drift Coordinate
        const float driftCoord = driftBoundaries.at(0) + ((driftBin + 0.5) * driftIntervalSigned); 

        // WireCoordinate
        const float wireCoord = wireBoundaries.at(0) + ((wireBin + 0.5) * wireIntervalSigned); 

        // NuVertex Separation
        const float deltaX = driftCoord - nuVertex2D.X();
        const float deltaZ = wireCoord - nuVertex2D.Z();
        const float nuSeparation = sqrt((deltaX * deltaX) + (deltaZ * deltaZ));

        dispGrid.SetGridEntry(driftBin, wireBin, nuSeparation);
    }

    return dispGrid;
}

/////////////////////////////////////////////////////////////

}
