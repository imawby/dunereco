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

#include "TPrincipal.h"

namespace HierarchyUtils
{

/////////////////////////////////////////////////////////////

void GetLinkInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const TVector3 &trueParentEndpoint, const TVector3 &trueChildStartpoint,  
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &linkVars)
{
    ////////////////////////////////////////////////
    // Set everything to default values
    ////////////////////////////////////////////////
    // General vars
    linkVars["ParentNuVertexSeparation"] = DEFAULT_DOUBLE;
    linkVars["ChildNuVertexSeparation"] = DEFAULT_DOUBLE;    
    linkVars["ParentBraggVariable"] = DEFAULT_DOUBLE;
    linkVars["OpeningAngle"] = DEFAULT_DOUBLE;
    linkVars["VertexSeparation"] = DEFAULT_DOUBLE;
    // End region vars
    linkVars["ParentEndRegionNHits"] = DEFAULT_DOUBLE;
    linkVars["ParentEndRegionNParticles"] = DEFAULT_DOUBLE;
    linkVars["ParentEndRegionRToWall"] = DEFAULT_DOUBLE;    
    // Parent/child directions
    linkVars["IsParentSet"] = false;
    linkVars["ReverseParent"] = false;    
    linkVars["ParentStartX"] = DEFAULT_DOUBLE; linkVars["ParentStartY"] = DEFAULT_DOUBLE; linkVars["ParentStartZ"] = DEFAULT_DOUBLE;    
    linkVars["ParentEndX"] = DEFAULT_DOUBLE; linkVars["ParentEndY"] = DEFAULT_DOUBLE; linkVars["ParentEndZ"] = DEFAULT_DOUBLE;    
    linkVars["ParentStartDX"] = DEFAULT_DOUBLE; linkVars["ParentStartDY"] = DEFAULT_DOUBLE; linkVars["ParentStartDZ"] = DEFAULT_DOUBLE;
    linkVars["ParentEndDX"] = DEFAULT_DOUBLE; linkVars["ParentEndDY"] = DEFAULT_DOUBLE; linkVars["ParentEndDZ"] = DEFAULT_DOUBLE;
    linkVars["IsChildSet"] = false;    
    linkVars["ReverseChild"] = false;    
    linkVars["ChildStartX"] = DEFAULT_DOUBLE; linkVars["ChildStartY"] = DEFAULT_DOUBLE; linkVars["ChildStartZ"] = DEFAULT_DOUBLE;    
    linkVars["ChildStartDX"] = DEFAULT_DOUBLE; linkVars["ChildStartDY"] = DEFAULT_DOUBLE; linkVars["ChildStartDZ"] = DEFAULT_DOUBLE;
    // Does/where child connects vars
    linkVars["OvershootStartDCA"] = DEFAULT_DOUBLE;
    linkVars["OvershootStartL"] = DEFAULT_DOUBLE;
    linkVars["OvershootEndDCA"] = DEFAULT_DOUBLE;
    linkVars["OvershootEndL"] = DEFAULT_DOUBLE;
    linkVars["DoesChildConnect"] = false;
    linkVars["ChildConnectionX"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionY"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionZ"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionDX"] = DEFAULT_DOUBLE; // I have named this var stupidly, it's actually the direction of the parent. 
    linkVars["ChildConnectionDY"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionDZ"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionDCA"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionExtrapDistance"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionL"] = DEFAULT_DOUBLE;
    linkVars["TrackLength"] = DEFAULT_DOUBLE;
    linkVars["ChildConnectionLRatio"] = DEFAULT_DOUBLE;
    // Splitting parent vars
    linkVars["ParentConnectionNUpstreamHits"] = DEFAULT_DOUBLE;
    linkVars["ParentConnectionNDownstreamHits"] = DEFAULT_DOUBLE;
    linkVars["ParentConnectionNHitRatio"] = DEFAULT_DOUBLE;
    linkVars["ParentConnectionEigenValueRatio"] = DEFAULT_DOUBLE;
    linkVars["ParentConnectionOpeningAngle"] = DEFAULT_DOUBLE;

    // If we can't find the connection pair then abort!
    if (!HierarchyUtils::CheatGetParentEndpointAndDirection(evt, parentPFP, trueParentEndpoint, recoModuleLabel, trackModuleLabel, linkVars))
        return;

    if (!HierarchyUtils::CheatGetChildStartpointAndDirection(evt, childPFP, trueChildStartpoint, recoModuleLabel, trackModuleLabel, linkVars))
        return;
    
    if (!linkVars["IsParentSet"] || !linkVars["IsChildSet"])
        return;

    // Set general vars
    const TVector3 parentStartpoint(linkVars["ParentStartX"], linkVars["ParentStartY"], linkVars["ParentStartZ"]);    
    const TVector3 parentEndpoint(linkVars["ParentEndX"], linkVars["ParentEndY"], linkVars["ParentEndZ"]);    
    const TVector3 childStartpoint(linkVars["ChildStartX"], linkVars["ChildStartY"], linkVars["ChildStartZ"]);    

    linkVars["ParentNuVertexSeparation"] = HierarchyUtils::GetNuVertexSeparation(evt, parentStartpoint, recoModuleLabel);
    linkVars["ChildNuVertexSeparation"] = HierarchyUtils::GetNuVertexSeparation(evt, childStartpoint, recoModuleLabel);
    linkVars["ParentBraggVariable"] = DEFAULT_DOUBLE; // TODO
    linkVars["VertexSeparation"] = (parentEndpoint - childStartpoint).Mag();

    // Set end region vars
    HierarchyUtils::GetEndRegionVars(evt, parentPFP, recoModuleLabel, linkVars);
    // Find if and where child would connect on parent
    HierarchyUtils::GetConnectionVars(evt, parentPFP, childPFP, recoModuleLabel, trackModuleLabel, linkVars);
    // Splitting parent vars
    HierarchyUtils::GetParentConnectionPointVars(evt, parentPFP, recoModuleLabel, 10.0, linkVars);    
}

/////////////////////////////////////////////////////////////

bool CheatGetParentEndpointAndDirection(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const TVector3 &trueParentEndpoint,
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &linkVars)
{
    if (!HierarchyUtils::IsPandoraApprovedTrack(evt, parentPFP, recoModuleLabel, trackModuleLabel))
        return false;

    linkVars["IsParentSet"] = true;

    const art::Ptr<recob::Track> &parentTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(parentPFP, evt, recoModuleLabel, trackModuleLabel);
    const float separationSq1 = (trueParentEndpoint - TVector3(parentTrack->Start().X(), parentTrack->Start().Y(), parentTrack->Start().Z())).Mag2();
    const float separationSq2 = (trueParentEndpoint - TVector3(parentTrack->End().X(), parentTrack->End().Y(), parentTrack->End().Z())).Mag2();

    TVector3 parentStartpoint = TVector3(0.f, 0.f, 0.f);
    TVector3 parentEndpoint = TVector3(0.f, 0.f, 0.f);
    TVector3 parentStartDirection = TVector3(0.f, 0.f, 0.f);
    TVector3 parentEndDirection = TVector3(0.f, 0.f, 0.f);

    if (separationSq1 < separationSq2)
    {
        linkVars["ReverseParent"] = true;
        parentStartpoint = TVector3(parentTrack->End().X(), parentTrack->End().Y(), parentTrack->End().Z());
        parentEndpoint = TVector3(parentTrack->Start().X(), parentTrack->Start().Y(), parentTrack->Start().Z());
        parentEndDirection = TVector3(parentTrack->StartDirection().X(), parentTrack->StartDirection().Y(), parentTrack->StartDirection().Z()) * (-1.0); // want direction to point out
        parentStartDirection = TVector3(parentTrack->EndDirection().X(), parentTrack->EndDirection().Y(), parentTrack->EndDirection().Z()) * (-1.0); // want direction to point along track
    }
    else
    {
        linkVars["ReverseParent"] = false;
        parentStartpoint = TVector3(parentTrack->Start().X(), parentTrack->Start().Y(), parentTrack->Start().Z());
        parentEndpoint = TVector3(parentTrack->End().X(), parentTrack->End().Y(), parentTrack->End().Z());
        parentEndDirection = TVector3(parentTrack->EndDirection().X(), parentTrack->EndDirection().Y(), parentTrack->EndDirection().Z());
        parentStartDirection = TVector3(parentTrack->StartDirection().X(), parentTrack->StartDirection().Y(), parentTrack->StartDirection().Z());
    }

    linkVars["ParentStartX"] = parentStartpoint.X();
    linkVars["ParentStartY"] = parentStartpoint.Y();
    linkVars["ParentStartZ"] = parentStartpoint.Z();
    linkVars["ParentStartDX"] = parentStartDirection.X();
    linkVars["ParentStartDY"] = parentStartDirection.Y();
    linkVars["ParentStartDZ"] = parentStartDirection.Z();
    linkVars["ParentEndX"] = parentEndpoint.X();
    linkVars["ParentEndY"] = parentEndpoint.Y();
    linkVars["ParentEndZ"] = parentEndpoint.Z();
    linkVars["ParentEndDX"] = parentEndDirection.X();
    linkVars["ParentEndDY"] = parentEndDirection.Y();
    linkVars["ParentEndDZ"] = parentEndDirection.Z();

    return true;
}

/////////////////////////////////////////////////////////////

bool CheatGetChildStartpointAndDirection(art::Event const & evt, const art::Ptr<recob::PFParticle> childPFP, const TVector3 &trueChildStartpoint,
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &linkVars)
{
    TVector3 childStartpoint = TVector3(0.f, 0.f, 0.f);
    TVector3 childDirection = TVector3(0.f, 0.f, 0.f);

    if (HierarchyUtils::IsPandoraApprovedTrack(evt, childPFP, recoModuleLabel, trackModuleLabel))
    {
        linkVars["IsChildSet"] = true;

        const art::Ptr<recob::Track> &childTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, evt, recoModuleLabel, trackModuleLabel);

        const float separationSq1 = (trueChildStartpoint - TVector3(childTrack->Start().X(), childTrack->Start().Y(), childTrack->Start().Z())).Mag2();
        const float separationSq2 = (trueChildStartpoint - TVector3(childTrack->End().X(), childTrack->End().Y(), childTrack->End().Z())).Mag2();

        if (separationSq1 < separationSq2)
        {
            linkVars["ReverseChild"] = false;
            childStartpoint = TVector3(childTrack->Start().X(), childTrack->Start().Y(), childTrack->Start().Z());
            childDirection = TVector3(childTrack->StartDirection().X(), childTrack->StartDirection().Y(), childTrack->StartDirection().Z());
        }
        else
        {
            linkVars["ReverseChild"] = true;
            childStartpoint = TVector3(childTrack->End().X(), childTrack->End().Y(), childTrack->End().Z());
            childDirection = TVector3(childTrack->EndDirection().X(), childTrack->EndDirection().Y(), childTrack->EndDirection().Z()) * (-1.0); // want dir to point along track
        }
    }
    else
    {
        try
        {
            const art::Ptr<recob::Vertex> &childRecobVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(childPFP, evt, recoModuleLabel);
            childStartpoint = TVector3(childRecobVertex->position().X(), childRecobVertex->position().Y(), childRecobVertex->position().Z());

            // I know the 25cm is hard coded... shhhh...
            if (!HierarchyUtils::GetParticleDirection(evt, childPFP, childStartpoint, recoModuleLabel, 25.0, childDirection))
                return false;

            linkVars["IsChildSet"] = true;
            linkVars["ReverseChild"] = false;
        }
        catch (...)
        {
            return false;
        }
    }

    linkVars["ChildStartX"] = childStartpoint.X();
    linkVars["ChildStartY"] = childStartpoint.Y();
    linkVars["ChildStartZ"] = childStartpoint.Z();
    linkVars["ChildStartDX"] = childDirection.X();
    linkVars["ChildStartDY"] = childDirection.Y();
    linkVars["ChildStartDZ"] = childDirection.Z();

    return true;
}
    
/////////////////////////////////////////////////////////////

double GetNuVertexSeparation(art::Event const & evt, const TVector3 &pfpVertex, const std::string recoModuleLabel)
{
    if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, recoModuleLabel))
        return DEFAULT_DOUBLE;

    const art::Ptr<recob::PFParticle> &nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, recoModuleLabel);

    try
    {
        const art::Ptr<recob::Vertex> &nuVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(nuPFP, evt, recoModuleLabel);

        const float dx = std::fabs(nuVertex->position().X() - pfpVertex.X());
        const float dy = std::fabs(nuVertex->position().Y() - pfpVertex.Y());
        const float dz = std::fabs(nuVertex->position().Z() - pfpVertex.Z());
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

void GetEndRegionVars(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP,  
    const std::string recoModuleLabel, std::map<std::string, double> &linkVars)
{
    const double separationThreshold = 5.0;
    
    HierarchyUtils::GetEndRegionNParticlesAndHits(evt, parentPFP, recoModuleLabel, separationThreshold, linkVars);
    HierarchyUtils::GetEndRegionRToWall(linkVars);
}

/////////////////////////////////////////////////////////////

void GetEndRegionNParticlesAndHits(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const std::string recoModuleLabel, 
    const double separationThreshold, std::map<std::string, double> &linkVars)
{
    const TVector3 parentEndpoint(linkVars["ParentEndX"], linkVars["ParentEndY"], linkVars["ParentEndZ"]);
    const std::vector<art::Ptr<recob::PFParticle>> &eventPFPs = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, recoModuleLabel);

    // Count hits near parent endpoint
    int hitCount = 0;
    int particleCount = 0;    
    double separationThresholdSq = separationThreshold * separationThreshold;
    
    for (const art::Ptr<recob::PFParticle> &eventPFP : eventPFPs)
    {
        if (eventPFP == parentPFP)
            continue;

        bool isClose = false;
        const std::vector<art::Ptr<recob::SpacePoint>> &pfpSpacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(eventPFP, evt, recoModuleLabel);

        for (const art::Ptr<recob::SpacePoint> &pfpSpacepoint : pfpSpacepoints)
        {
            const TVector3 pos = TVector3(pfpSpacepoint->XYZ()[0], pfpSpacepoint->XYZ()[1], pfpSpacepoint->XYZ()[2]);
            const double sepSq = (pos - parentEndpoint).Mag2();

            if (sepSq < separationThresholdSq)
            {
                isClose = true;
                ++hitCount;
            }
        }

        if (isClose)
            ++particleCount;
    }

    linkVars["ParentEndRegionNHits"] = hitCount;
    linkVars["ParentEndRegionNParticles"] = particleCount;    
}

/////////////////////////////////////////////////////////////

void GetEndRegionRToWall(std::map<std::string, double> &linkVars)
{
    const TVector3 parentEndpoint(linkVars["ParentEndX"], linkVars["ParentEndY"], linkVars["ParentEndZ"]);
    
    // Detector boundaries
    const std::vector<std::vector<float>> detectorBoundaries = {
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
            const float sep = std::fabs(parentEndpoint(iView) - detectorBoundaries.at(iView).at(iBoundary));

            closestDistance = std::min(closestDistance, sep);
        }
    }

    linkVars["ParentEndRegionRToWall"] = closestDistance;
}

/////////////////////////////////////////////////////////////

void GetConnectionVars(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &linkVars)
{
    const art::Ptr<recob::Track> &parentTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(parentPFP, evt, recoModuleLabel, trackModuleLabel);    
    const TVector3 parentStartpoint(linkVars["ParentStartX"], linkVars["ParentStartY"], linkVars["ParentStartZ"]);    
    const TVector3 parentEndpoint(linkVars["ParentEndX"], linkVars["ParentEndY"], linkVars["ParentEndZ"]);    
    const TVector3 parentStartDirection(linkVars["ParentStartDX"], linkVars["ParentStartDY"], linkVars["ParentStartDZ"]);
    const TVector3 parentEndDirection(linkVars["ParentEndDX"], linkVars["ParentEndDY"], linkVars["ParentEndDZ"]);
    const TVector3 childVertex(linkVars["ChildStartX"], linkVars["ChildStartY"], linkVars["ChildStartZ"]);    
    const TVector3 childDirection(linkVars["ChildStartDX"], linkVars["ChildStartDY"], linkVars["ChildStartDZ"]);

    ////////////////////////////////////////////////    
    // Look for the child connection point variables
    ////////////////////////////////////////////////    
    const int nTrajPoints = parentTrack->NumberTrajectoryPoints();

    if (nTrajPoints < 1)
        return;
    
    double cumulativeL = 0.0;
    double minDCA = std::numeric_limits<double>::max();

    // Loop through trajectory points - will loop from track start to track end (so we may have to invert it at the end)
    for (int i = 0; i < (nTrajPoints - 1); ++i)
    {
        const bool arePointsValid = (parentTrack->HasValidPoint(i)) && (parentTrack->HasValidPoint(i+1));

        if (!arePointsValid)
            continue;

        const TVector3 firstPoint = TVector3(parentTrack->TrajectoryPoint(i).position.X(), parentTrack->TrajectoryPoint(i).position.Y(), parentTrack->TrajectoryPoint(i).position.Z());
        const TVector3 secondPoint = TVector3(parentTrack->TrajectoryPoint(i+1).position.X(), parentTrack->TrajectoryPoint(i+1).position.Y(), parentTrack->TrajectoryPoint(i+1).position.Z());
        const TVector3 midPoint = (secondPoint + firstPoint) * 0.5;
        
        // Extrapolate the child to the parent
        TVector3 extrapolationPoint(0.f, 0.f, 0.f);
        const bool extrapolateBackwards = HierarchyUtils::ExtrapolateChildToParent(midPoint, childVertex, childDirection, extrapolationPoint);
        const float thisDCA = (midPoint - extrapolationPoint).Mag();

        if (thisDCA < minDCA)
        {
            minDCA = thisDCA;

            // Does this connect - be very loose?
            const float buffer = 50.0;
            if (HierarchyUtils::DoesConnect(firstPoint, secondPoint, extrapolationPoint, buffer))
            {
                linkVars["DoesChildConnect"] = true;
                linkVars["ChildConnectionDCA"] = thisDCA;
                linkVars["ChildConnectionL"] = cumulativeL + (midPoint - firstPoint).Mag();
                linkVars["ChildConnectionX"] = midPoint.X();
                linkVars["ChildConnectionY"] = midPoint.Y();
                linkVars["ChildConnectionZ"] = midPoint.Z();
                linkVars["ChildConnectionExtrapDistance"] = (extrapolationPoint - childVertex).Mag() * (extrapolateBackwards ? 1.0 : (-1.0));

                const TVector3 midPointDir = (secondPoint - firstPoint).Unit();

                linkVars["ChildConnectionDX"] = midPointDir.X() * (linkVars["ReverseParent"] ? (-1.0) : 1.0);
                linkVars["ChildConnectionDY"] = midPointDir.Y() * (linkVars["ReverseParent"] ? (-1.0) : 1.0);
                linkVars["ChildConnectionDZ"] = midPointDir.Z() * (linkVars["ReverseParent"] ? (-1.0) : 1.0);
            }
        }

        cumulativeL += (secondPoint - firstPoint).Mag();
    }

    linkVars["TrackLength"] = cumulativeL;

    if (linkVars["DoesChildConnect"])
    {
        if (linkVars["ReverseParent"] == true)
            linkVars["ChildConnectionL"] = linkVars["TrackLength"] - linkVars["ChildConnectionL"];

        linkVars["ChildConnectionLRatio"] = (linkVars["TrackLength"]  < std::numeric_limits<double>::epsilon()) ? 
            DEFAULT_DOUBLE : linkVars["ChildConnectionL"] / linkVars["TrackLength"];

        linkVars["OpeningAngle"] = TVector3(linkVars["ChildConnectionDX"], linkVars["ChildConnectionDY"], linkVars["ChildConnectionDZ"]).Angle(childDirection) * 180.0 / M_PI;
    }

    // If it didn't connect, fill out 'Overshoot' info
    if (!linkVars["DoesChildConnect"])
    {
        TVector3 extrapolationPoint_start = TVector3(0.0, 0.0, 0.0);
        const bool extrapolateBackwards_start = HierarchyUtils::ExtrapolateChildToParent(parentStartpoint, childVertex, childDirection, extrapolationPoint_start);

        linkVars["OvershootStartDCA"] = (extrapolationPoint_start - parentStartpoint).Mag() * (extrapolateBackwards_start ? 1.0 : (-1.0));
        linkVars["OvershootStartL"] = std::fabs((childVertex - parentStartpoint).Dot(parentStartDirection));
 
        TVector3 extrapolationPoint_end = TVector3(0.0, 0.0, 0.0);
        const bool extrapolateBackwards_end = HierarchyUtils::ExtrapolateChildToParent(parentEndpoint, childVertex, childDirection, extrapolationPoint_end);

        linkVars["OvershootEndDCA"] = (extrapolationPoint_end - parentEndpoint).Mag() * (extrapolateBackwards_end ? 1.0 : (-1.0));
        linkVars["OvershootEndL"] = std::fabs((childVertex - parentEndpoint).Dot(parentEndDirection));
    }
}


/////////////////////////////////////////////////////////////

bool IsPandoraApprovedTrack(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const std::string recoModuleLabel, 
    const std::string trackModuleLabel)
{
    const bool isPandoraTrack = (HierarchyUtils::GetTrackScore(evt, pfp, recoModuleLabel) > 0.5);

    if (!isPandoraTrack)
        return false;

    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, recoModuleLabel, trackModuleLabel))
        return false;

    return true;
}
    
/////////////////////////////////////////////////////////////

bool ExtrapolateChildToParent(const TVector3 &parentPosition, const TVector3 &childVertex, const TVector3 &childDirection,
    TVector3 &extrapolationPoint)
{                          
    const double extrapDistance((parentPosition - childVertex).Dot(childDirection));
    extrapolationPoint = childVertex + (extrapDistance * childDirection);

    // make sure we extrapolate the child backwards
    if (extrapDistance > 0.0)
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

bool DoesConnect(const TVector3 &boundary1, const TVector3 &boundary2, const TVector3 &testPoint, const float buffer)
{
    const TVector3 displacement = boundary2 - boundary1;
    const double segmentLength = displacement.Mag();
    const TVector3 unitVector = displacement * (1.0 / segmentLength);

    const double l = unitVector.Dot(testPoint - boundary1);

    if ((l < 0.0) || (l > segmentLength))
        return false;

    const double t = unitVector.Cross(testPoint - boundary1).Mag2();

    if (t > (buffer * buffer))
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

void GetParentConnectionPointVars(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP,
    const std::string recoModuleLabel, const double searchRegion, std::map<std::string, double> &linkVars)
{
    // Does it connect?
    if (!linkVars["DoesChildConnect"])
        return;

    // Find the connection point on the parent
    const TVector3 parentConnectionPoint = TVector3(linkVars["ChildConnectionX"], 
        linkVars["ChildConnectionY"], linkVars["ChildConnectionZ"]);
    const TVector3 parentConnectionDir = TVector3(linkVars["ChildConnectionDX"], 
        linkVars["ChildConnectionDY"], linkVars["ChildConnectionDZ"]);

    // Use direction to current direction to split spacepoints into two groups
    const std::vector<art::Ptr<recob::SpacePoint>> parentSPs = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(parentPFP, evt, recoModuleLabel); 

    if (parentSPs.empty())
        return;

    std::vector<art::Ptr<recob::SpacePoint>> upstreamGroup;
    std::vector<art::Ptr<recob::SpacePoint>> downstreamGroup;

    for (const art::Ptr<recob::SpacePoint> &parentSP : parentSPs)
    {
        const TVector3 parentSPPos = TVector3(parentSP->position().X(), parentSP->position().Y(), parentSP->position().Z());

        const float dx = std::fabs(parentSPPos.X() - parentConnectionPoint.X());
        const float dy = std::fabs(parentSPPos.Y() - parentConnectionPoint.Y());
        const float dz = std::fabs(parentSPPos.Z() - parentConnectionPoint.Z());
        const float sepSq = (dx * dx) + (dy * dy) + (dz * dz);

        if (sepSq > (searchRegion * searchRegion))
            continue;

        // Assign to group
        const double thisL = parentConnectionDir.Dot(parentSPPos - parentConnectionPoint);

        std::vector<art::Ptr<recob::SpacePoint>> &group = (thisL > 0) ? downstreamGroup : upstreamGroup;
        group.push_back(parentSP);
    }

    linkVars["ParentConnectionNUpstreamHits"] = upstreamGroup.size();
    linkVars["ParentConnectionNDownstreamHits"] = downstreamGroup.size();

    if ((upstreamGroup.size() == 0) || (downstreamGroup.size() == 0))
        return;

    linkVars["ParentConnectionNHitRatio"] = linkVars["ParentConnectionNDownstreamHits"] / linkVars["ParentConnectionNUpstreamHits"];

    // Now do PCA
    // Demand a reasonable number of hits
    if ((upstreamGroup.size() < 3) || (downstreamGroup.size() < 3))
        return;

    std::vector<double> upstreamEigenvalues;
    std::vector<TVector3> upstreamEigenvectors;
    HierarchyUtils::RunPCA(upstreamGroup, upstreamEigenvalues, upstreamEigenvectors);

    // Make sure the PCA fit is sensible
    for (double upstreamEigenvalue : upstreamEigenvalues)
        if (std::isnan(upstreamEigenvalue))
            return;

    std::vector<double> downstreamEigenvalues;
    std::vector<TVector3> downstreamEigenvectors;
    HierarchyUtils::RunPCA(downstreamGroup, downstreamEigenvalues, downstreamEigenvectors);

    // Make sure the PCA fit is sensible
    for (double downstreamEigenvalue : downstreamEigenvalues)
        if (std::isnan(downstreamEigenvalue))
            return;

    // Get opening angle from first eigenvectors (this is the longitudinal one) - straight would be around 180
    linkVars["ParentConnectionOpeningAngle"] = upstreamEigenvectors.at(0).Angle(downstreamEigenvectors.at(0)) * 180.0 / M_PI;

    // Get average transverse eigenvalues, get ratio
    const double upstreamAvTransverseE = (upstreamEigenvalues.at(1) + upstreamEigenvalues.at(2)) / 2.0;
    const double downstreamAvTransverseE = (downstreamEigenvalues.at(1) + downstreamEigenvalues.at(2)) / 2.0;

    linkVars["ParentConnectionEigenValueRatio"] = upstreamAvTransverseE < std::numeric_limits<double>::epsilon() ? 0.0 : (downstreamAvTransverseE / upstreamAvTransverseE);
}

/////////////////////////////////////////////////////////////

void RunPCA(const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints, std::vector<double> &eVals, std::vector<TVector3> &eVecs)
{
    TPrincipal* principal = new TPrincipal(3, "D");

    for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
        principal->AddRow(spacepoint->XYZ());

    // PERFORM PCA
    principal->MakePrincipals();
    // GET EIGENVALUES AND EIGENVECTORS
    for (unsigned int i = 0; i < 3; ++i)
         eVals.push_back(principal->GetEigenValues()->GetMatrixArray()[i]);

    for (int i : {0, 3, 6})
    {
        const double eVec_x = principal->GetEigenVectors()->GetMatrixArray()[i];
        const double eVec_y = principal->GetEigenVectors()->GetMatrixArray()[i + 1];
        const double eVec_z = principal->GetEigenVectors()->GetMatrixArray()[i + 2];

        eVecs.push_back(TVector3(eVec_x, eVec_y, eVec_z));
    }
}

/////////////////////////////////////////////////////////////

bool GetParticleDirection(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfp, const TVector3 &pfpVertex, const std::string recoModuleLabel, 
    const float searchRegion, TVector3 &pfpDirection) 
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, recoModuleLabel); 

    if (spacepoints.empty())
        return false;

    int nBins = 180;
    float angleMin = 0.f, angleMax = 2.f * M_PI;
    float binWidth = (angleMax - angleMin) / static_cast<float>(nBins);

    std::vector<std::vector<int>> spatialDist(nBins, std::vector<int>(nBins, 0));
    std::vector<std::vector<float>> energyDist(nBins, std::vector<float>(nBins, 0.f));

    // theta0YZ then theta0XZ
    int highestSP = 0;
    float highestEnergy = 0.f;
    int bestTheta0YZBin = -1;
    int bestTheta0XZBin = -1; 

    for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
    {
        const TVector3 spacepointPos = TVector3(spacepoint->position().X(), spacepoint->position().Y(), spacepoint->position().Z());
        const TVector3 displacement = spacepointPos - pfpVertex;
        const float mag = displacement.Mag();

        if (mag > searchRegion)
            continue;

        const float magXZ = sqrt((displacement.X() * displacement.X()) + (displacement.Z() * displacement.Z()));

        float theta0YZ = (mag < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.Y() / mag) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f : 
             std::acos(displacement.Y() / mag);

        float theta0XZ = (magXZ < std::numeric_limits<float>::epsilon()) ? 0.f : 
            (std::fabs(std::fabs(displacement.X() / magXZ) - 1.f) < std::numeric_limits<float>::epsilon()) ? 0.f :
             std::acos(displacement.X() / magXZ);

        // try do signed-ness
        if (displacement.Z() < 0.f)
           theta0XZ = (2.0 * M_PI) - theta0XZ;

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

    pfpDirection = TVector3(std::sin(bestTheta0YZ) * std::cos(bestTheta0XZ), std::cos(bestTheta0YZ), std::sin(bestTheta0YZ) * std::sin(bestTheta0XZ));

    return true;
}

/////////////////////////////////////////////////////////////
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

double GetSeparation3D(art::Event const & evt, const art::Ptr<recob::PFParticle> pfp, const TVector3 &position,
    const std::string recoModuleLabel)
{
    const std::vector<art::Ptr<recob::SpacePoint>> sps = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfp, evt, recoModuleLabel);

    if (sps.empty())
        return DEFAULT_DOUBLE;

    double closestDistanceSq = std::numeric_limits<double>::max();

    for (const art::Ptr<recob::SpacePoint> &sp : sps)
    {
        const double dx = sp->XYZ()[0] - position.X();
        const double dy = sp->XYZ()[1] - position.Y();
        const double dz = sp->XYZ()[2] - position.Z();
        const double separation = (dx * dx) + (dy * dy) + (dz * dz);

        if (separation < closestDistanceSq)
            closestDistanceSq = separation;
    }

    return std::sqrt(closestDistanceSq);
}


/////////////////////////////////////////////////////////////

double GetChargeRatio(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    const double parentCharge = HierarchyUtils::GetPFPCharge(evt, parentPFP, recoModuleLabel);
    const double childCharge = HierarchyUtils::GetPFPCharge(evt, childPFP, recoModuleLabel);
    // ratio = child / parent
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

int GetPIDLinkTypeWithIvysaurus(const int parentParticleType, const int childParticleType)
{
    int absParentPDG = DEFAULT_DOUBLE;
    int absChildPDG = DEFAULT_DOUBLE;

    std::vector<int> pdg = {13, 2212, 211, 11, 22};

    if ((parentParticleType >= 0) && (parentParticleType <= 4))
        absParentPDG = pdg.at(parentParticleType);

    if ((childParticleType >= 0) && (childParticleType <= 4))
        absChildPDG = pdg.at(childParticleType);

    return HierarchyUtils::GetPIDLinkTypeWithPDG(absParentPDG, absChildPDG);
}

/////////////////////////////////////////////////////////////

int GetPIDLinkTypeWithPDG(const int parentPDG, const int childPDG)
{
    const int absParentPDG = std::abs(parentPDG);
    const int absChildPDG = std::abs(childPDG);

    if ((parentPDG == DEFAULT_DOUBLE) || (childPDG == DEFAULT_DOUBLE))
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

int GetTrackShowerLinkType(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel)
{
    const double parentTrackScore = HierarchyUtils::GetTrackScore(evt, parentPFP, recoModuleLabel);
    const double childTrackScore = HierarchyUtils::GetTrackScore(evt, childPFP, recoModuleLabel);

    if ((parentTrackScore < 0.0) || (childTrackScore < 0.0))
        return -1;

    const bool isParentTrack = parentTrackScore > 0.5;
    const bool isChildTrack = childTrackScore > 0.5;

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

}



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////// IDEAS!! /////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/*
bool GetParentChildVerticesAndDirections(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &linkVars)
{
    //////////////////////////
    // Positions
    //////////////////////////
    std::vector<TVector3> parentPositions;
    std::vector<TVector3> childPositions;    
    std::vector<TVector3> childPositions_temp;

    // Parent
    if (!HierarchyUtils::IsPandoraApprovedTrack(evt, parentPFP, recoModuleLabel, trackModuleLabel))
        return false;

    const art::Ptr<recob::Track> &parentTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(parentPFP, evt, recoModuleLabel, trackModuleLabel);

    parentPositions.push_back(TVector3(parentTrack->Start().X(), parentTrack->Start().Y(), parentTrack->Start().Z()));
    parentPositions.push_back(TVector3(parentTrack->End().X(), parentTrack->End().Y(), parentTrack->End().Z()));    

    // Child
    if (HierarchyUtils::IsPandoraApprovedTrack(evt, childPFP, recoModuleLabel, trackModuleLabel))
    {
        const art::Ptr<recob::Track> &childTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, evt, recoModuleLabel, trackModuleLabel);

        childPositions_temp.push_back(TVector3(childTrack->Start().X(), childTrack->Start().Y(), childTrack->Start().Z()));
        childPositions_temp.push_back(TVector3(childTrack->End().X(), childTrack->End().Y(), childTrack->End().Z()));    
    }
    else
    {
        try
        {
            const art::Ptr<recob::Vertex> &childRecobVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(childPFP, evt, recoModuleLabel);
            childPositions_temp.push_back(TVector3(childRecobVertex->position().X(), childRecobVertex->position().Y(), childRecobVertex->position().Z()));
        }
        catch (...)
        {
            return false;
        }
    }

    //////////////////////////
    // Directions
    //////////////////////////
    std::vector<TVector3> parentDirections;
    std::vector<TVector3> childDirections;    

    // Parent
    // make sure start points out    
    parentDirections.push_back(TVector3(parentTrack->StartDirection().X(), parentTrack->StartDirection().Y(), parentTrack->StartDirection().Z()) * -1.0); 
    parentDirections.push_back(TVector3(parentTrack->EndDirection().X(), parentTrack->EndDirection().Y(), parentTrack->EndDirection().Z()));    

    // Child
    for (const TVector3 &childPosition : childPositions_temp)
    {
        TVector3 childDirection(0.f, 0.f, 0.f);

        // I know the 25cm is hard coded... shhhh...
        if (!HierarchyUtils::GetParticleDirection(evt, childPFP, childPosition, recoModuleLabel, 25.0, childDirection))
            continue;

        childPositions.push_back(childPosition);
        childDirections.push_back(childDirection);
    }

    // Did we actually find any child directions?
    if (childDirections.empty())
        return false;

    //////////////////////////
    // Work out the best pair
    //////////////////////////
    double closestSep = std::numeric_limits<double>::max();
    int closestParentIndex = -1;
    int closestChildIndex = -1;    
    
    for (unsigned int iParent = 0; iParent < parentPositions.size(); ++iParent)
    {
        const TVector3 &parentPosition = parentPositions.at(iParent);
        
        for (unsigned int iChild = 0; iChild < childPositions.size(); ++iChild)
        {
            std::cout << "iParent: " << iParent << std::endl;
            std::cout << "iChild: " << iChild << std::endl;

            const TVector3 &childPosition = childPositions.at(iChild);
            const TVector3 &childDirection = childDirections.at(iChild);            

            TVector3 extrapolationPoint(0.f, 0.f, 0.f);
            if (!ExtrapolateChildToParent(parentPosition, childPosition, childDirection, extrapolationPoint))
            {
                std::cout << "CANNOT EXTRAOPLATE!" << std::endl;
                continue;
            }

            std::cout << "CAN EXTRAPOLATE!" << std::endl;

            const double sepSq = (extrapolationPoint - parentPosition).Mag2();

            if (sepSq < closestSep)
            {
                closestSep = sepSq;
                closestParentIndex = iParent;
                closestChildIndex = iChild;
            }
        }
    }

    //////////////////////////
    // Set best pair
    //////////////////////////
    if ((closestParentIndex == -1) || (closestChildIndex == -1))
        return false;

    linkVars["IsParentSet"] = true;
    linkVars["ReverseParent"] = closestParentIndex == 0; // Is the start actually the endpoint?
    //linkVars["ParentX"] = parentPositions.at(closestParentIndex).X();
    //linkVars["ParentY"] = parentPositions.at(closestParentIndex).Y();
    //linkVars["ParentZ"] = parentPositions.at(closestParentIndex).Z();
    //linkVars["ParentDX"] = parentDirections.at(closestParentIndex).X(); 
    //linkVars["ParentDY"] = parentDirections.at(closestParentIndex).Y();
    //linkVars["ParentDZ"] = parentDirections.at(closestParentIndex).Z(); // need to change this
    linkVars["IsChildSet"] = true;
    linkVars["ChildStartX"] = childPositions.at(closestChildIndex).X();
    linkVars["ChildStartY"] = childPositions.at(closestChildIndex).Y();
    linkVars["ChildStartZ"] = childPositions.at(closestChildIndex).Z();
    linkVars["ChildStartDX"] = childDirections.at(closestChildIndex).X(); 
    linkVars["ChildStartDY"] = childDirections.at(closestChildIndex).Y();
    linkVars["ChildStartDZ"] = childDirections.at(closestChildIndex).Z();

    return true;    
}
*/
