#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaSpacePointUtils.h"
#include "dunereco/FDSelections/HierarchyUtils.h"

namespace HierarchyUtils
{

/////////////////////////////////////////////////////////////

double GetTrackScore(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel)
{
    double trackScore = DEFAULT_DOUBLE;

    try
    {
        const art::Ptr<larpandoraobj::PFParticleMetadata> metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfp, evt, recoModuleLabel);

        if (metadata->GetPropertiesMap().find("TrackScore") != metadata->GetPropertiesMap().end())
            trackScore = metadata->GetPropertiesMap().at("TrackScore");
    }
    catch (...) {}

    return trackScore;
}

/////////////////////////////////////////////////////////////

double GetNuVertexSeparation(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, 
    const std::string recoModuleLabel)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, recoModuleLabel))
        return DEFAULT_DOUBLE;

    const art::Ptr<recob::PFParticle> &nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, recoModuleLabel);

    try
    {
        const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, recoModuleLabel);
        const art::Ptr<recob::Vertex> &pfpVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, recoModuleLabel);

        const float dx = std::fabs(nuVertex->position().X() - pfpVertex->position().X());
        const float dy = std::fabs(nuVertex->position().Y() - pfpVertex->position().Y());
        const float dz = std::fabs(nuVertex->position().Z() - pfpVertex->position().Z());
        const float sep = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));

        return sep;
    }
    catch (...)
    {
        return DEFAULT_DOUBLE;
    }
}

/////////////////////////////////////////////////////////////

double GetBraggVariable()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

int GetEndRegionNHits(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string trackModuleLabel,  const std::string recoModuleLabel, 
    const double separationThreshold)
{
    // If shower-like then return -1
    const double trackScore = HierarchyUtils::GetTrackScore(evt, pfp, recoModuleLabel);

    if (trackScore < 0.5)
        return -1;

    // Get endpoint
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, recoModuleLabel, trackModuleLabel))
        return -1;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, recoModuleLabel, trackModuleLabel);
    const TVector3 trackEndpoint = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    // Get pfps from event
    const std::vector<art::Ptr<recob::PFParticle>> &eventPFPs = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, recoModuleLabel);

    int count = 0;

    for (const art::Ptr<recob::PFParticle> &eventPFP : eventPFPs)
    {
        if (eventPFP == pfp)
            continue;

        const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(eventPFP, evt, recoModuleLabel);

        for (const art::Ptr<recob::SpacePoint> &pfpSpacepoint : pfpSpacepoints)
        {
            const TVector3 pos = TVector3(pfpSpacepoint->XYZ()[0], pfpSpacepoint->XYZ()[1], pfpSpacepoint->XYZ()[2]);
            const double sep = (pos - trackEndpoint).Mag();

            if (sep < separationThreshold)
                ++count;
        }
    }

    return count;
}

/////////////////////////////////////////////////////////////

int GetEndRegionNParticles(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string trackModuleLabel,  const std::string recoModuleLabel, 
    const double separationThreshold)
{
    // If shower-like then return 0
    const double trackScore = HierarchyUtils::GetTrackScore(evt, pfp, recoModuleLabel);

    if (trackScore < 0.5)
        return 0;

    // Get endpoint
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, recoModuleLabel, trackModuleLabel))
        return 0;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, recoModuleLabel, trackModuleLabel);
    const TVector3 trackEndpoint = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    // Get pfps from event
    const std::vector<art::Ptr<recob::PFParticle>> &eventPFPs = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, recoModuleLabel);

    int count = 0;

    for (const art::Ptr<recob::PFParticle> &eventPFP : eventPFPs)
    {
        if (eventPFP == pfp)
            continue;

        const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(eventPFP, evt, recoModuleLabel);

        for (const art::Ptr<recob::SpacePoint> &pfpSpacepoint : pfpSpacepoints)
        {
            const TVector3 pos = TVector3(pfpSpacepoint->XYZ()[0], pfpSpacepoint->XYZ()[1], pfpSpacepoint->XYZ()[2]);
            const double sep = (pos - trackEndpoint).Mag();

            if (sep < separationThreshold)
            {
                ++count;
                break;
            }
        }
    }

    return count;
}

/////////////////////////////////////////////////////////////

double GetEndRegionRToWall(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp,  
    const std::string recoModuleLabel, const std::string trackModuleLabel)
{
    // If shower-like then return -1
    const double trackScore = HierarchyUtils::GetTrackScore(evt, pfp, recoModuleLabel);

    if (trackScore < 0.5)
        return -1;

    // Get endpoint
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, recoModuleLabel, trackModuleLabel))
        return -1;

    const art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, recoModuleLabel, trackModuleLabel);
    const TVector3 trackEndpoint = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    // Detector boundaries
    std::vector<std::vector<float>> detectorBoundaries = {
        std::vector<float>({-360.0, 360.0}),
        std::vector<float>({-600.0, 600.0}),
        std::vector<float>({0.0, 1394.0})
    };

    // Get closest distance
    float closestDistance = std::numeric_limits<float>::max();

    for (int iView : {0, 1, 2})
    {
        for (float iBoundary : {0, 1})
        {
            const float sep = std::fabs(trackEndpoint(iView) - detectorBoundaries.at(iView).at(iBoundary));

            closestDistance = std::min(closestDistance, sep);
        }
    }

    return closestDistance;
}

/////////////////////////////////////////////////////////////

double GetVertexSeparation(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    double sep = DEFAULT_DOUBLE;

    try
    {
        const art::Ptr<recob::Vertex> parentVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(parentPFP, evt, recoModuleLabel);
        const art::Ptr<recob::Vertex> childVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(childPFP, evt, recoModuleLabel);

        const TVector3 parentPos = TVector3(parentVertex->position().X(), parentVertex->position().Y(), parentVertex->position().Z());
        const TVector3 childPos = TVector3(childVertex->position().X(), childVertex->position().Y(), childVertex->position().Z());

        sep = (parentPos - childPos).Mag();
    }
    catch (...)
    {
    }

    return sep;
}

/////////////////////////////////////////////////////////////

double GetSeparation3D(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    const std::vector<art::Ptr<recob::SpacePoint>> parentSPs = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(parentPFP, evt, recoModuleLabel);
    const std::vector<art::Ptr<recob::SpacePoint>> childSPs = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(childPFP, evt, recoModuleLabel);

    if (parentSPs.empty() || childSPs.empty())
        return DEFAULT_DOUBLE;

    double closestDistanceSq = std::numeric_limits<double>::max();

    for (const art::Ptr<recob::SpacePoint> &parentSP : parentSPs)
    {
        for (const art::Ptr<recob::SpacePoint> &childSP : childSPs)
        {
            const double dx = parentSP->XYZ()[0] - childSP->XYZ()[0];
            const double dy = parentSP->XYZ()[1] - childSP->XYZ()[1];
            const double dz = parentSP->XYZ()[2] - childSP->XYZ()[2];
            const double separation = (dx * dx) + (dy * dy) + (dz * dz);

            if (separation < closestDistanceSq)
                closestDistanceSq = separation;
        }
    }

    return std::sqrt(closestDistanceSq);
}

/////////////////////////////////////////////////////////////
// ratio = child / parent

double GetChargeRatio(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    const double parentCharge = HierarchyUtils::GetPFPCharge(evt, parentPFP, recoModuleLabel);
    const double childCharge = HierarchyUtils::GetPFPCharge(evt, childPFP, recoModuleLabel);
    const double ratio = (parentCharge < std::numeric_limits<float>::epsilon()) ? 0.0 : (childCharge / parentCharge);

    return ratio;
}

/////////////////////////////////////////////////////////////

double GetPFPCharge(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel)
{
    std::vector<art::Ptr<recob::Hit>> collectionViewHits = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfp, evt, recoModuleLabel, 2);

    double totalCharge = 0.0;

    for (const art::Ptr<recob::Hit> &hit : collectionViewHits)
        totalCharge += hit->Integral();

    return totalCharge;
}

/////////////////////////////////////////////////////////////

int GetPIDLinkType(const int parentPDG, const int childPDG)
{
    const int absParentPDG = std::abs(parentPDG);
    const int absChildPDG = std::abs(childPDG);

    if ((absParentPDG < 0) || (absChildPDG < 0))
        return - 1;
    else if ((absParentPDG == 13) && (absChildPDG == 13))     // mu - mu
        return 0;
    else if ((absParentPDG == 13) && (absChildPDG == 2212))   // mu - p
        return 1;
    else if ((absParentPDG == 13) && (absChildPDG == 211))    // mu - pi
        return 2;
    else if ((absParentPDG == 13) && (absChildPDG == 11))     // mu - e
        return 3;
    else if ((absParentPDG == 13) && (absChildPDG == 22))     // mu - gamma
        return 4;
    else if ((absParentPDG == 2212) && (absChildPDG == 13))   // p - mu
        return 5;
    else if ((absParentPDG == 2212) && (absChildPDG == 2212)) // p - p
        return 6;
    else if ((absParentPDG == 2212) && (absChildPDG == 211))  // p - pi
        return 7;
    else if ((absParentPDG == 2212) && (absChildPDG == 11))   // p - e
        return 8;
    else if ((absParentPDG == 2212) && (absChildPDG == 22))   // p - gamma
        return 9;
    else if ((absParentPDG == 211) && (absChildPDG == 13))    // pi - mu
        return 10;
    else if ((absParentPDG == 211) && (absChildPDG == 2212))  // pi - p
        return 11;
    else if ((absParentPDG == 211) && (absChildPDG == 211))   // pi - pi
        return 12;
    else if ((absParentPDG == 211) && (absChildPDG == 11))    // pi - e
        return 13;
    else if ((absParentPDG == 211) && (absChildPDG == 22))    // pi - gamma
        return 14;
    else if ((absParentPDG == 11) && (absChildPDG == 13))     // e - mu
        return 15;
    else if ((absParentPDG == 11) && (absChildPDG == 2212))   // e - p
        return 16;
    else if ((absParentPDG == 11) && (absChildPDG == 211))    // e - pi
        return 17;
    else if ((absParentPDG == 11) && (absChildPDG == 11))     // e - e
        return 18;
    else if ((absParentPDG == 11) && (absChildPDG == 22))     // e - gamma
        return 19;
    else if ((absParentPDG == 22) && (absChildPDG == 13))     // gamma - mu
        return 20;
    else if ((absParentPDG == 22) && (absChildPDG == 2212))   // gamma - p
        return 21;
    else if ((absParentPDG == 22) && (absChildPDG == 211))    // gamma - pi
        return 22;
    else if ((absParentPDG == 22) && (absChildPDG == 11))     // gamma - e
        return 23;
    else if ((absParentPDG == 22) && (absChildPDG == 22))     // gamma - gamma
        return 24;
    else
        return 25;
}

/////////////////////////////////////////////////////////////

double GetOpeningAngle(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    TVector3 parentDirection(0.f, 0.f, 0.f);
    TVector3 childDirection(0.f, 0.f, 0.f);

    // I know the 25cm is hard coded... shhhh...
    if (!HierarchyUtils::GetParticleDirection(evt, parentPFP, recoModuleLabel, 25.0, parentDirection))
        return DEFAULT_DOUBLE;

    // I know the 25cm is hard coded... shhhh...
    if (!HierarchyUtils::GetParticleDirection(evt, childPFP, recoModuleLabel, 25.0, childDirection))
        return DEFAULT_DOUBLE;

    const double openingAngle = parentDirection.Angle(childDirection) * 180.0 / M_PI;

    return openingAngle;
}

/////////////////////////////////////////////////////////////

bool GetParticleDirection(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfp, const std::string recoModuleLabel, 
    const float searchRegion, TVector3 &pfpDirection) 
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, recoModuleLabel); 

    if (spacepoints.empty())
        return false;

    try
    {
        const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, recoModuleLabel);
        const TVector3 vertexPos = TVector3(vertex->position().X(), vertex->position().Y(), vertex->position().Z());

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
            const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z());
            const TVector3 displacement = spacepointPos - vertexPos;
            const float mag = sqrt((displacement.X() * displacement.X()) + (displacement.Y() * displacement.Y()) + (displacement.Z() * displacement.Z()));

            if (mag > searchRegion)
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

            const std::vector<art::Ptr<recob::Hit>> assocHits = dune_ana::DUNEAnaSpacePointUtils::GetHits(spacepoint, evt, recoModuleLabel);

            if (assocHits.empty())
                continue;

            const int bin0YZ = std::floor(theta0YZ / binWidth);
            const int bin0XZ = std::floor(theta0XZ / binWidth);

            spatialDist[bin0YZ][bin0XZ] += 1;
            energyDist[bin0YZ][bin0XZ] += assocHits.front()->Integral(); // very basic tie-breaker..

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

        pfpDirection = TVector3(std::fabs(std::sin(bestTheta0YZ) * std::sin(bestTheta0XZ)), std::fabs(std::cos(bestTheta0YZ)), 
            std::fabs(std::sin(bestTheta0YZ) * std::cos(bestTheta0XZ)));

        if (bestTheta0XZ > M_PI)
            pfpDirection.SetX(pfpDirection.X() * -1.f);

        if (bestTheta0YZ > M_PI)
            pfpDirection.SetZ(pfpDirection.Z() * -1.f);

        if ((bestTheta0YZ > (M_PI / 2.f)) && (bestTheta0YZ < (M_PI * 3.f / 2.f)))
            pfpDirection.SetY(pfpDirection.Y() * -1.f);
    }
    catch (...) { return false; }

    return true;
}

/////////////////////////////////////////////////////////////

int GetTrackShowerLinkType(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    const double parentTrackScore = HierarchyUtils::GetTrackScore(evt, parentPFP, recoModuleLabel);
    const double childTrackScore = HierarchyUtils::GetTrackScore(evt, childPFP, recoModuleLabel);

    if ((parentTrackScore < 0.0) || (childTrackScore < 0.0))
        return -1;

    const bool isParentTrack = parentTrackScore > 0.5;
    const bool isChildTrack = parentTrackScore > 0.5;

    if (isParentTrack && isChildTrack)
        return 0;
    else if (isParentTrack && !isChildTrack)
        return 1;
    else if (!isParentTrack && isChildTrack)
        return 2;
    else 
        return 3;
}

/////////////////////////////////////////////////////////////











/*
void EdgeVars::SetSeparation(art::Event const & evt)
{
    std::vector<TVector3> parentPositions;

    // If neutrino...
    if ((std::abs(m_parentPFP->PdgCode()) == 12) || (std::abs(m_parentPFP->PdgCode()) == 14) || (std::abs(m_parentPFP->PdgCode()) == 16))
    {
        try
        {
            const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(m_parentPFP, evt, m_recoModuleLabel);
            parentPositions.push_back(TVector3(nuVertex->position().X(), nuVertex->position().Y(), nuVertex->position().Z()));
        }
        catch (...) { return; }
    }
    else
    {
        const std::vector<art::Ptr<recob::SpacePoint>> parentSPs = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(m_parentPFP, evt, m_recoModuleLabel);

        for (const art::Ptr<recob::SpacePoint> &parentSP : parentSPs)
            parentPositions.push_back(TVector3(parentSP->XYZ()[0], parentSP->XYZ()[1], parentSP->XYZ()[2]));
    }

    const std::vector<art::Ptr<recob::SpacePoint>> childSPs = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(m_childPFP, evt, m_recoModuleLabel);

    if (parentPositions.empty() || childSPs.empty())
        return;

    double closestDistanceSq = std::numeric_limits<double>::max();

    for (const TVector3 &parentPosition : parentPositions)
    {
        for (const art::Ptr<recob::SpacePoint> &childSP : childSPs)
        {
            const double dx = parentPosition.X() - childSP->XYZ()[0];
            const double dy = parentPosition.Y() - childSP->XYZ()[1];
            const double dz = parentPosition.Z() - childSP->XYZ()[2];
            const double separation = (dx * dx) + (dy * dy) + (dz * dz);

            if (separation < closestDistanceSq)
                closestDistanceSq = separation;
        }
    }

    m_separation = std::sqrt(closestDistanceSq);
}
*/
/////////////////////////////////////////////////////////////
/*
void EdgeVars::SetOpeningAngle(art::Event const & evt)
{

}
*/
/////////////////////////////////////////////////////////////
/*

*/
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

}




