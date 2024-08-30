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

// Implement reverse by working out which end is closer?

void GetLinkConnectionInfo(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP, const art::Ptr<recob::PFParticle> childPFP, 
    const std::string recoModuleLabel, const std::string trackModuleLabel, std::map<std::string, double> &connectionVars)
{
    connectionVars["UnderOvershootDCA"] = (-1.0) * DEFAULT_DOUBLE; // needs to be positive for below function to work
    connectionVars["UnderOvershootL"] = (-1.0) * DEFAULT_DOUBLE; // needs to be positive for below function to work
    connectionVars["DoesChildConnect"] = false;
    connectionVars["ChildConnectionX"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionY"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionZ"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionDX"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionDY"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionDZ"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionDCA"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionExtrapDistance"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionL"] = DEFAULT_DOUBLE;
    connectionVars["TrackLength"] = DEFAULT_DOUBLE;
    connectionVars["ChildConnectionLRatio"] = DEFAULT_DOUBLE;

    // First work out the parent and child vertices





    
    TVector3 childVertex = TVector3(0.0, 0.0, 0.0);

    // Get child direction
    try
    {
        const art::Ptr<recob::Vertex> &childRecobVertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(childPFP, evt, recoModuleLabel);
        childVertex = TVector3(childRecobVertex->position().X(), childRecobVertex->position().Y(), childRecobVertex->position().Z());
    }
    catch (...)
    {
        return;
    }

    // Get parent track, and check that it's good before we do the heavy lifting
    if (!HierarchyUtils::IsPandoraApprovedTrack(evt, parentPFP, recoModuleLabel, trackModuleLabel))
        return;

    const art::Ptr<recob::Track> &parentTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(parentPFP, evt, recoModuleLabel, trackModuleLabel);

    const int nTrajPoints = parentTrack->NumberTrajectoryPoints();

    if (nTrajPoints < 1)
        return;

    // Alright heavy lifting, get child direction
    TVector3 childDirection(0.f, 0.f, 0.f);

    // I know the 25cm is hard coded... shhhh...
    if (!HierarchyUtils::GetParticleDirection(evt, childPFP, recoModuleLabel, 25.0, childDirection))
        return;

    // Now loop through trajectory points - we can work out several things here...
    double cumulativeL = 0.0;
    double minDCA = std::numeric_limits<double>::max();

    // Loop through trajectory points
    for (int i = 0; i < (nTrajPoints - 1); ++i)
    {
        const bool arePointsValid = (parentTrack->HasValidPoint(i)) && (parentTrack->HasValidPoint(i+1));

        if (!arePointsValid)
            continue;

        const TVector3 firstPoint = TVector3(parentTrack->TrajectoryPoint(i).position.X(), parentTrack->TrajectoryPoint(i).position.Y(), parentTrack->TrajectoryPoint(i).position.Z());
        const TVector3 secondPoint = TVector3(parentTrack->TrajectoryPoint(i+1).position.X(), parentTrack->TrajectoryPoint(i+1).position.Y(), parentTrack->TrajectoryPoint(i+1).position.Z());
        const TVector3 midPoint = (secondPoint + firstPoint) * 0.5;

        // Extrapolate the child to the parent
        const double extrapFactor((midPoint - childVertex).Dot(childDirection));
        const TVector3 extrapPoint = childVertex + (extrapFactor * childDirection);

        // Work out DCA
         const float thisDCA = (midPoint - extrapPoint).Mag();

        if (thisDCA < minDCA)
        {
            // Extrap point should move closer to parent track
            if (extrapFactor > 0.0)
                continue;

            minDCA = thisDCA;

            // Does this connect - be very loose?
            const float buffer = 5.0;
            if (HierarchyUtils::IsInBoundingBox(firstPoint, secondPoint, extrapPoint, buffer))
            {
                connectionVars["DoesChildConnect"] = true;
                connectionVars["ChildConnectionDCA"] = thisDCA;
                connectionVars["ChildConnectionL"] = cumulativeL + (midPoint - firstPoint).Mag();
                connectionVars["ChildConnectionX"] = midPoint.X();
                connectionVars["ChildConnectionY"] = midPoint.Y();
                connectionVars["ChildConnectionZ"] = midPoint.Z();
                connectionVars["ChildConnectionExtrapDistance"] = extrapFactor;

                const TVector3 midPointDir = (secondPoint - firstPoint).Unit();

                connectionVars["ChildConnectionDX"] = midPointDir.X();
                connectionVars["ChildConnectionDY"] = midPointDir.Y();
                connectionVars["ChildConnectionDZ"] = midPointDir.Z();
            }
        }

        // Whilst we're here, we might as well work out the track length
        cumulativeL += (secondPoint - firstPoint).Mag();
    }

    connectionVars["TrackLength"] = cumulativeL;

    if (connectionVars["DoesChildConnect"])
        connectionVars["ChildConnectionLRatio"] = (connectionVars["TrackLength"]  < std::numeric_limits<double>::epsilon()) ? 
            DEFAULT_DOUBLE : connectionVars["ChildConnectionL"] / connectionVars["TrackLength"];

    if (!connectionVars["DoesChildConnect"])
    {
        const TVector3 parentTrackStart = TVector3(parentTrack->Start().X(), parentTrack->Start().Y(), parentTrack->Start().Z());
        const TVector3 parentTrackEnd = TVector3(parentTrack->End().X(), parentTrack->End().Y(), parentTrack->End().Z());
        const double gamma_start((childVertex - parentTrackStart).Dot(childDirection));
        const double gamma_end((childVertex - parentTrackEnd).Dot(childDirection));

        // Child should be coming out of the track
        if (gamma_start < 0.0)
        {
            const TVector3 extrap_start = childVertex + (gamma_start * childDirection);
            const double dca_start = (extrap_start - parentTrackStart).Mag();
            connectionVars["UnderOvershootDCA"] = std::min(dca_start, connectionVars["UnderOvershootDCA"]);

            // Get L on parent track
            const TVector3 parentTrackStartDirection = TVector3(parentTrack->StartDirection().X(), parentTrack->StartDirection().Y(), parentTrack->StartDirection().Z());
            const double parentL_start = std::fabs((childVertex - parentTrackStart).Dot(parentTrackStartDirection));
            connectionVars["UnderOvershootL"] = std::min(parentL_start, connectionVars["UnderOvershootL"]);
        }

        if (gamma_end < 0.0)
        {
            const TVector3 extrap_end = childVertex + (gamma_end * childDirection);
            const double dca_end = (extrap_end - parentTrackEnd).Mag();
            connectionVars["UnderOvershootDCA"] = std::min(dca_end, connectionVars["UnderOvershootDCA"]);

            // Get L on parent track
            const TVector3 parentTrackEndDirection = TVector3(parentTrack->EndDirection().X(), parentTrack->EndDirection().Y(), parentTrack->EndDirection().Z());
            const double parentL_end = std::fabs((childVertex - parentTrackEnd).Dot(parentTrackEndDirection));
            connectionVars["UnderOvershootL"] = std::min(parentL_end, connectionVars["UnderOvershootL"]);
        }
    }

    // I know the 10cm is hard coded... shhhh...
    HierarchyUtils::GetParentConnectionPointVars(evt, parentPFP, recoModuleLabel, 10.0, connectionVars);
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

bool IsInBoundingBox(const TVector3 &boundary1, const TVector3 &boundary2, const TVector3 &testPoint, const float buffer)
{
    /*
    const double minX = std::min(boundary1.X(), boundary2.X()) - buffer;
    const double maxX = std::max(boundary1.X(), boundary2.X()) + buffer;
    const double minY = std::min(boundary1.Y(), boundary2.Y()) - buffer;
    const double maxY = std::max(boundary1.Y(), boundary2.Y()) + buffer;
    const double minZ = std::min(boundary1.Z(), boundary2.Z()) - buffer;
    const double maxZ = std::max(boundary1.Z(), boundary2.Z()) + buffer;

    if ((testPoint.X() < minX) || (testPoint.X() > maxX))
        return false;

    if ((testPoint.Y() < minY) || (testPoint.Y() > maxY))
        return false;

    if ((testPoint.Z() < minZ) || (testPoint.Z() > maxZ))
        return false;
    */

    const TVector3 displacement = boundary2 - boundary1;
    const double segmentLength = displacement.Mag();
    const TVector3 unitVector = displacement * (1.0 / segmentLength);

    const double l = unitVector.Dot(testPoint - boundary1);

    if ((l < 0.0) || (l > segmentLength))
        return false;

    const double t = unitVector.Cross(testPoint - boundary1).Mag();

    if (t > buffer)
        return false;

    return true;
}

/////////////////////////////////////////////////////////////

void GetParentConnectionPointVars(art::Event const & evt, const art::Ptr<recob::PFParticle> parentPFP,
    const std::string recoModuleLabel, const double searchRegion, std::map<std::string, double> &connectionVars)
{
    // Ratios taken downstream/upstream
    connectionVars["ParentConnectionNUpstreamHits"] = DEFAULT_DOUBLE;
    connectionVars["ParentConnectionNDownstreamHits"] = DEFAULT_DOUBLE;
    connectionVars["ParentConnectionNHitRatio"] = DEFAULT_DOUBLE;
    connectionVars["ParentConnectionEigenValueRatio"] = DEFAULT_DOUBLE;
    connectionVars["ParentConnectionOpeningAngle"] = DEFAULT_DOUBLE;

    // Does it connect?
    if (!connectionVars["DoesChildConnect"])
        return;

    // Find the connection point on the parent
    const TVector3 parentConnectionPoint = TVector3(connectionVars["ChildConnectionX"], 
        connectionVars["ChildConnectionY"], connectionVars["ChildConnectionZ"]);
    const TVector3 parentConnectionDir = TVector3(connectionVars["ChildConnectionDX"], 
        connectionVars["ChildConnectionDY"], connectionVars["ChildConnectionDZ"]);

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

    connectionVars["ParentConnectionNUpstreamHits"] = upstreamGroup.size();
    connectionVars["ParentConnectionNDownstreamHits"] = downstreamGroup.size();

    if ((upstreamGroup.size() == 0) || (downstreamGroup.size() == 0))
        return;

    connectionVars["ParentConnectionNHitRatio"] = connectionVars["ParentConnectionNDownstreamHits"] / connectionVars["ParentConnectionNUpstreamHits"];

    // Now do PCA
    std::vector<double> upstreamEigenvalues;
    std::vector<TVector3> upstreamEigenvectors;
    HierarchyUtils::RunPCA(upstreamGroup, upstreamEigenvalues, upstreamEigenvectors);

    std::vector<double> downstreamEigenvalues;
    std::vector<TVector3> downstreamEigenvectors;
    HierarchyUtils::RunPCA(downstreamGroup, downstreamEigenvalues, downstreamEigenvectors);

    // Get opening angle from first eigenvectors (this is the longitudinal one)
    connectionVars["ParentConnectionOpeningAngle"] = upstreamEigenvectors.at(0).Angle(downstreamEigenvectors.at(0)) * 180.0 / M_PI;

    // Get average transverse eigenvalues, get ratio
    const double upstreamAvTransverseE = (upstreamEigenvalues.at(1) + upstreamEigenvalues.at(2)) / 2.0;
    const double downstreamAvTransverseE = (downstreamEigenvalues.at(1) + downstreamEigenvalues.at(2)) / 2.0;

    connectionVars["ParentConnectionEigenValueRatio"] = downstreamAvTransverseE / upstreamAvTransverseE;
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

        /*
        if (bestTheta0XZ > M_PI)
            pfpDirection.SetX(pfpDirection.X() * -1.f);

        if (bestTheta0YZ > M_PI)
            pfpDirection.SetZ(pfpDirection.Z() * -1.f);

        if ((bestTheta0YZ > (M_PI / 2.f)) && (bestTheta0YZ < (M_PI * 3.f / 2.f)))
            pfpDirection.SetY(pfpDirection.Y() * -1.f);
        */
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




