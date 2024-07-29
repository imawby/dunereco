#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

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

double GetBraggVariable()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

int GetEndRegionNHits()
{
    return DEFAULT_INT;
}

/////////////////////////////////////////////////////////////

int GetEndRegionNParticles()
{
    return DEFAULT_INT;
}

/////////////////////////////////////////////////////////////

double GetEndRegionRToWall()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

double GetVertexSeparation()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

double GetSeparation3D()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

double GetEnergyRatio()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

int GetPIDLinkType()
{
    return DEFAULT_INT;
}

/////////////////////////////////////////////////////////////

double GetOpeningAngle()
{
    return DEFAULT_DOUBLE;
}

/////////////////////////////////////////////////////////////

int GetTrackShowerLinkType()
{
    return DEFAULT_INT;
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
    TVector3 parentDirection(0.f, 0.f, 0.f);
    TVector3 childDirection(0.f, 0.f, 0.f);

    // I know the 25cm is hard coded... shhhh...
    if (!this->GetParticleDirection(evt, m_parentPFP, 25.0, parentDirection))
        return;

    // I know the 25cm is hard coded... shhhh...
    if (!this->GetParticleDirection(evt, m_childPFP, 25.0, childDirection))
        return;

    m_openingAngle = parentDirection.Angle(childDirection);
}
*/
/////////////////////////////////////////////////////////////
/*
bool EdgeVars::GetParticleDirection(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfp, const float searchRegion, 
    TVector3 &pfpDirection)
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, m_recoModuleLabel); 

    if (spacepoints.empty())
        return false;

    try
    {
        const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfp, evt, m_recoModuleLabel);
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

            const std::vector<art::Ptr<recob::Hit>> assocHits = dune_ana::DUNEAnaSpacePointUtils::GetHits(spacepoint, evt, m_recoModuleLabel);

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
*/
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

}




