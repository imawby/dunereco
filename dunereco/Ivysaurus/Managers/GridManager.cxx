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

#include "dunereco/Ivysaurus/Managers/GridManager.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"


namespace ivysaurus
{

GridManager::Grid::Grid(const TVector3 origin, const float driftSpan, const float wireSpan, 
    const unsigned int dimensions, const GridManager::PandoraView pandoraView, const bool isInitialised) : 
        m_axisDimensions(dimensions),
        m_pandoraView(pandoraView),
        m_isInitialised(isInitialised)
{
    m_gridValues = std::vector<std::vector<float>>(m_axisDimensions, std::vector<float>(m_axisDimensions, 0.0));
    m_countValues = std::vector<std::vector<float>>(m_axisDimensions, std::vector<float>(m_axisDimensions, 0.0));

    const float driftInterval = driftSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_driftBoundaries.push_back(origin.X() + (i * driftInterval));

    const float wireInterval = wireSpan / m_axisDimensions;

    for (unsigned int i = 0; i <= m_axisDimensions; ++i)
        m_wireBoundaries.push_back(origin.Z() + (i * wireInterval));

    m_isNormalised = false;
}

/////////////////////////////////////////////////////////////

bool GridManager::Grid::IsInsideGrid(const TVector3 &position) const
{
    const float lowerDriftCoord(std::min(m_driftBoundaries.front(), m_driftBoundaries.back()));
    const float higherDriftCoord(std::max(m_driftBoundaries.front(), m_driftBoundaries.back()));

    if ((position.X() < lowerDriftCoord) || (position.X() > higherDriftCoord))
        return false;

    const float lowerWireCoord(std::min(m_wireBoundaries.front(), m_wireBoundaries.back()));
    const float higherWireCoord(std::max(m_wireBoundaries.front(), m_wireBoundaries.back()));

    if ((position.Z() < lowerWireCoord) || (position.Z() > higherWireCoord))
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

void GridManager::Grid::AddToGrid(const TVector3 &position, const float energy, const float weight)
{
    // Get drift bin
    const float driftInterval = std::fabs(m_driftBoundaries.at(0) - m_driftBoundaries.at(1));
    const unsigned int driftBin = std::floor(std::fabs(position.X() - m_driftBoundaries.front()) / driftInterval); 

    if (driftBin >= m_axisDimensions)
        throw cet::exception("ivysaur::GridManager") << "the bin does not exist!";

    // Get wire bin
    const float wireInterval = std::fabs(m_wireBoundaries.at(0) - m_wireBoundaries.at(1));
    const unsigned int wireBin = std::floor(std::fabs(position.Z() - m_wireBoundaries.front()) / wireInterval); 

    if (wireBin >= m_axisDimensions)
        throw cet::exception("ivysaur::GridManager") << "the bin does not exist!";

    // Now fill grid
    m_gridValues[driftBin][wireBin] += energy;
    m_countValues[driftBin][wireBin] += weight;
}


/////////////////////////////////////////////////////////////

void GridManager::Grid::NormaliseGrid()
{
    if (m_isNormalised)
        throw cet::exception("ivysaur::GridManager") << "this grid is already normalised!";

    for (unsigned int driftIndex = 0; driftIndex < m_axisDimensions; ++driftIndex)
    {
        for (unsigned int wireIndex = 0; wireIndex < m_axisDimensions; ++wireIndex)
        {
            m_gridValues[driftIndex][wireIndex] /= m_countValues[driftIndex][wireIndex];
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
    m_gridSize3D(pset.get<float>("GridSize3D")),
    m_dimensions(pset.get<float>("GridDimensions")),
    m_uWireAngle(pset.get<float>("UWireAngle")), //radians
    m_vWireAngle(pset.get<float>("VWireAngle")), // radians
    m_wWireAngle(pset.get<float>("WWireAngle")) // radians
{
}

/////////////////////////////////////////////////////////////

GridManager::~GridManager()
{
}

/////////////////////////////////////////////////////////////

const GridManager::Grid GridManager::ObtainViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    const GridManager::PandoraView pandoraView) const
{
    // We need a pfparticle direction, let's take the track direction
    // I'm assuming that pfps have been fitted as tracks & showers
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_recoModuleLabel))
        return Grid(TVector3(0.f, 0.f, 0.f), 0.f, 0.f, 0, pandoraView, false);

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_recoModuleLabel);
    const TVector3 startDirection = TVector3(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
    const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);
    const TVector3 startPosition = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    const float diagonalLength = sqrt(2.0 * (m_gridSize3D * m_gridSize3D));
    const TVector3 extrapolatedPosition = startPosition + (startDirection * diagonalLength);

    // Now need to project these things into the 'Pandora view'
    const TVector3 projectedStart = ProjectIntoPandoraView(startPosition, pandoraView);
    const TVector3 projectedExtrapolated = ProjectIntoPandoraView(extrapolatedPosition, pandoraView);
    const float driftSpan = projectedExtrapolated.X() - projectedStart.X();
    const float wireSpan = projectedExtrapolated.Z() - projectedStart.Z();

    return Grid(projectedStart, driftSpan, wireSpan, m_dimensions, pandoraView, true);
}

/////////////////////////////////////////////////////////////

void GridManager::FillViewGrid(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    GridManager::Grid &grid) const
{
    // Get 2D hits for Pandora view..
    const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);

    PandoraView pandoraView = grid.GetPandoraView();

    for (const art::Ptr<recob::Hit> hit : pfpHits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(hit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;
                     
        const PandoraView thisPandoraView = GetPandora2DView(hit);

        if (thisPandoraView != pandoraView)
            continue;

        // Get 2D hit position
        const TVector3 pandoraHitPosition = ObtainPandoraHitPosition(evt, hit, pandoraView);

        if (!grid.IsInsideGrid(pandoraHitPosition))
            continue;

        // TODO - these are just fillers
        const float energy = 1.f;
        const float weight = 1.f;

        grid.AddToGrid(pandoraHitPosition, energy, weight);
    }

    grid.NormaliseGrid();
}

/////////////////////////////////////////////////////////////

const TVector3 GridManager::ObtainPandoraHitPosition(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
    const GridManager::PandoraView hitType) const
{
    art::ServiceHandle<geo::Geometry const> theGeometry;
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const geo::WireID hitWireID = hit->WireID();
    const geo::CryostatID cryostatID(hitWireID.Cryostat);
    const double hitTime = hit->PeakTime();
    const double xCoord = detProp.ConvertTicksToX(hitTime, hitWireID.Plane, hitWireID.TPC, hitWireID.Cryostat);

    // Get hit Y and Z coordinates, based on central position of wire
    geo::Point_t xyz = theGeometry->Cryostat(cryostatID).TPC(hitWireID.TPC).Plane(hitWireID.Plane).Wire(hitWireID.Wire).GetCenter();

    return TVector3(xCoord, 0.f, hitType == TPC_VIEW_U ? YZToU(xyz.Y(), xyz.Z()) : hitType == TPC_VIEW_V ? YZToV(xyz.Y(), xyz.Z()) : YZToW(xyz.Y(), xyz.Z()));
}

/////////////////////////////////////////////////////////////

const TVector3 GridManager::ProjectIntoPandoraView(const TVector3 &inputPosition3D, const GridManager::PandoraView pandoraView) const
{
    const float xCoord = inputPosition3D.X();
    const float yCoord = inputPosition3D.Y();
    const float zCoord = inputPosition3D.Z();

    return TVector3(xCoord, 0.f, pandoraView == TPC_VIEW_U ? YZToU(yCoord, zCoord) : pandoraView == TPC_VIEW_V ? YZToV(yCoord, zCoord) : YZToW(yCoord, zCoord));
}

/////////////////////////////////////////////////////////////

float GridManager::YZToU(const float yCoord, const float zCoord) const
{
    return (zCoord * std::cos(m_uWireAngle)) - (yCoord * std::sin(m_uWireAngle));
}

/////////////////////////////////////////////////////////////

float GridManager::YZToV(const float yCoord, const float zCoord) const
{
    return (zCoord * std::cos(m_vWireAngle)) - (yCoord * std::sin(m_vWireAngle));
}

/////////////////////////////////////////////////////////////

float GridManager::YZToW(const float yCoord, const float zCoord) const
{
    return (zCoord * std::cos(m_wWireAngle)) - (yCoord * std::sin(m_wWireAngle));
}

/////////////////////////////////////////////////////////////

const GridManager::PandoraView GridManager::GetPandora2DView(const art::Ptr<recob::Hit> &hit) const
{
    const geo::WireID hitWireID(hit->WireID());
    const geo::View_t hitView(hit->View());
    const geo::View_t thisPandoraView(lar_pandora::LArPandoraGeometry::GetGlobalView(hitWireID.Cryostat, hitWireID.TPC, hitView));

    if (thisPandoraView == geo::kW || thisPandoraView == geo::kY)
        return TPC_VIEW_W;
    else if (thisPandoraView == geo::kU)
        return TPC_VIEW_U;
    else if (thisPandoraView == geo::kV)
        return TPC_VIEW_V;
    else
        throw cet::exception("ivysaur::GridManager") << "wire view not recognised";
}

/////////////////////////////////////////////////////////////

}
