////////////////////////////////////////////////////////////////////////
/// \file    ProfileManager.cxx
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

#include "dunereco/Ivysaurus/Managers/ProfileManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaSpacePointUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

namespace ivysaurus
{

/////////////////////////////////////////////////////////////

ProfileManager::EnergyProfiles::EnergyProfiles(const std::vector<float> &longitudinal, const std::vector<float> &transverse, 
    const std::vector<float> &energy) :
        m_longitudinal(longitudinal),
        m_transverse(transverse),
        m_energy(energy)
{
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProfileManager::ProfileManager(const fhicl::ParameterSet& pset) :
    m_hitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_trackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    m_showerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    m_recombFactor(pset.get<float>("RecombFactor")),
    m_calorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
}

/////////////////////////////////////////////////////////////

ProfileManager::~ProfileManager()
{
}

/////////////////////////////////////////////////////////////

void ProfileManager::GetProfiles(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    const int directionMode, std::vector<float> &longitudinalProfile, std::vector<float> &transverseProfile, 
    std::vector<float> &energies, std::vector<std::vector<float>> &fitPositions) const
{
    // Just incase
    longitudinalProfile.clear();
    transverseProfile.clear();
    energies.clear();
    fitPositions.clear();

    // Method 0 - true particles
    if (directionMode == 0)
    {
        TVector3 truePosition = TVector3(0.f, 0.f, 0.f);
        TVector3 trueDirection = TVector3(0.f, 0.f, 0.f);

        if (!GetTrueDirectionAndPosition(evt, pfparticle, trueDirection, truePosition))
            return;

        EnergyProfiles energyProfiles = CreateProfiles(evt, pfparticle, trueDirection, truePosition);
        longitudinalProfile = energyProfiles.m_longitudinal;
        transverseProfile = energyProfiles.m_transverse;
        energies = energyProfiles.m_energy;

        FillFitPositions(evt, pfparticle, trueDirection, truePosition, fitPositions);
    }

    // Method 1 - reco start direction
    else if (directionMode == 1)
    {
        TVector3 recoPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 recoStartDirection = TVector3(0.f, 0.f, 0.f);

        if (!GetRecoStartDirectionAndPosition(evt, pfparticle, recoStartDirection, recoPosition))
            return;

        EnergyProfiles energyProfiles = CreateProfiles(evt, pfparticle, recoStartDirection, recoPosition);
        longitudinalProfile = energyProfiles.m_longitudinal;
        transverseProfile = energyProfiles.m_transverse;
        energies = energyProfiles.m_energy;

        FillFitPositions(evt, pfparticle, recoStartDirection, recoPosition, fitPositions);
    }

    // Method 2 - PCA
    else if (directionMode == 2)
    {
        TVector3 pcaPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 pcaDirection = TVector3(0.f, 0.f, 0.f);

        if (!GetPCADirectionAndPosition(evt, pfparticle, pcaDirection, pcaPosition))
            return;

        EnergyProfiles energyProfiles = CreateProfiles(evt, pfparticle, pcaDirection, pcaPosition);
        longitudinalProfile = energyProfiles.m_longitudinal;
        transverseProfile = energyProfiles.m_transverse;
        energies = energyProfiles.m_energy;

        FillFitPositions(evt, pfparticle, pcaDirection, pcaPosition, fitPositions);
    }

    // Method 3 - Initial track
    else if (directionMode == 3)
    {
        TVector3 trackStubPosition = TVector3(0.f, 0.f, 0.f);
        TVector3 trackStubDirection = TVector3(0.f, 0.f, 0.f);

        if (!GetTrackStubDirectionAndPosition(evt, pfparticle, trackStubDirection, trackStubPosition))
            return;

        EnergyProfiles energyProfiles = CreateProfiles(evt, pfparticle, trackStubDirection, trackStubPosition);
        longitudinalProfile = energyProfiles.m_longitudinal;
        transverseProfile = energyProfiles.m_transverse;
        energies = energyProfiles.m_energy;

        FillFitPositions(evt, pfparticle, trackStubDirection, trackStubPosition, fitPositions);

    }

    // Method 4 - Sliding linear fit (essentially a track fit?)
    else if (directionMode == 4)
    {
        EnergyProfiles energyProfiles = CreateProfilesFromTrack(evt, pfparticle);
        longitudinalProfile = energyProfiles.m_longitudinal;
        transverseProfile = energyProfiles.m_transverse;
        energies = energyProfiles.m_energy;

        FillFitPositionsFromTrack(evt, pfparticle, fitPositions);
    }
}

/////////////////////////////////////////////////////////////

bool ProfileManager::GetTrueDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    TVector3 &trueDirection, TVector3 &truePosition) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    const std::vector<art::Ptr<recob::Hit>> pfpHits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfparticle, evt, m_recoModuleLabel);
    int g4ID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfpHits, 1);

    if (!TruthMatchUtils::Valid(g4ID))
        return false;

    art::ServiceHandle<cheat::ParticleInventoryService> piServ;
    const simb::MCParticle* matchedMC = piServ->ParticleList().at(g4ID);

    if (!matchedMC)
        return false;

    trueDirection = TVector3(matchedMC->Momentum().X(), matchedMC->Momentum().Y(), matchedMC->Momentum().Z());

    if (trueDirection.Mag() < std::numeric_limits<float>::epsilon())
        return false;

    trueDirection = trueDirection.Unit();

    //////////////////////
    truePosition = TVector3(matchedMC->Vx(), matchedMC->Vy(), matchedMC->Vz());

    // If the true particle is a photon, we need to be smarter.
    // Select the first traj point where tne photon loses energy
    if (abs(matchedMC->PdgCode()) == 22)
    {
        const double initialEnergy = matchedMC->E();
        const unsigned int nTrajPoints = matchedMC->NumberTrajectoryPoints();

        for (unsigned int trajPoint = 0; trajPoint < nTrajPoints; trajPoint++)
        {
            if (matchedMC->E(trajPoint) < initialEnergy)
            {
                truePosition = matchedMC->Position(trajPoint).Vect();
                break;
            }
        }
    }

    return true;
}

/////////////////////////////////////////////////////////////

bool ProfileManager::GetRecoStartDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    TVector3 &recoStartDirection, TVector3 &recoPosition) const
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel); 
    const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);
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

        // ISOBEL THIS IS A PLACEHOLDER
        if (mag > 25.0)
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

    recoStartDirection = TVector3(std::fabs(std::sin(bestTheta0YZ) * std::sin(bestTheta0XZ)), std::fabs(std::cos(bestTheta0YZ)), 
        std::fabs(std::sin(bestTheta0YZ) * std::cos(bestTheta0XZ)));

    if (bestTheta0XZ > M_PI)
        recoStartDirection.SetX(recoStartDirection.X() * -1.f);

    if (bestTheta0YZ > M_PI)
        recoStartDirection.SetZ(recoStartDirection.Z() * -1.f);

    if ((bestTheta0YZ > (M_PI / 2.f)) && (bestTheta0YZ < (M_PI * 3.f / 2.f)))
        recoStartDirection.SetY(recoStartDirection.Y() * -1.f);

    if (recoStartDirection.Mag() < std::numeric_limits<float>::epsilon())
        return false;

    recoPosition = vertexPos;
    recoStartDirection = recoStartDirection.Unit();

    return true;
}

/////////////////////////////////////////////////////////////

bool ProfileManager::GetPCADirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    TVector3 &pcaDirection, TVector3 &pcaPosition) const
{
    const std::vector<art::Ptr<recob::SpacePoint>> spacepoints = dune_ana::DUNEAnaPFParticleUtils::GetSpacePoints(pfparticle, evt, m_recoModuleLabel); 
    const art::Ptr<recob::Vertex> vertex = dune_ana::DUNEAnaPFParticleUtils::GetVertex(pfparticle, evt, m_recoModuleLabel);

    // To use Pandora's PCA methods, we need CartesianVectors
    const pandora::CartesianVector vertexPosition(vertex->position().X(), vertex->position().Y(), vertex->position().Z());
    pandora::CartesianPointVector positionVector;
    for (const art::Ptr<recob::SpacePoint> &spacepoint : spacepoints)
        positionVector.emplace_back(pandora::CartesianVector(spacepoint->XYZ()[0], spacepoint->XYZ()[1], spacepoint->XYZ()[2]));

    // Now do PCA
    try {
        // Access centroid of shower via this method
        const lar_content::LArShowerPCA initialLArShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(positionVector, vertexPosition));

        // Ensure successful creation of all structures before placing results in output containers, remaking LArShowerPCA with updated vertex
        const pandora::CartesianVector& centroid(initialLArShowerPCA.GetCentroid());
        const pandora::CartesianVector& primaryAxis(initialLArShowerPCA.GetPrimaryAxis());

        // Project the PFParticle vertex onto the PCA axis
        const pandora::CartesianVector projectedVertex(centroid - primaryAxis.GetUnitVector() * (centroid - vertexPosition).GetDotProduct(primaryAxis));

        // By convention, principal axis should always point away from vertex
        const float testProjection(primaryAxis.GetDotProduct(projectedVertex - centroid));
        const float directionSF((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

        pcaPosition = TVector3(projectedVertex.GetX(), projectedVertex.GetY(), projectedVertex.GetZ());
        pcaDirection = TVector3(primaryAxis.GetX() * directionSF, primaryAxis.GetY() * directionSF, primaryAxis.GetZ() * directionSF);
    }
    catch (...) { return false; }

    return true;
}


/////////////////////////////////////////////////////////////

bool ProfileManager::GetTrackStubDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    TVector3 &trackStubDirection, TVector3 &trackStubPosition) const
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

    trackStubPosition = TVector3(initialTrackStub->Vertex().X(), initialTrackStub->Vertex().Y(), initialTrackStub->Vertex().Z());
    trackStubDirection = TVector3(initialTrackStub->StartDirection().X(), initialTrackStub->StartDirection().Y(), initialTrackStub->StartDirection().Z());

    return true;
}

/////////////////////////////////////////////////////////////

ProfileManager::EnergyProfiles ProfileManager::CreateProfiles(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    const TVector3 &fittingAxisDirection, const TVector3 &fittingAxisStart) const
{
    const std::vector<art::Ptr<recob::Hit>> &collectionViewHits = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfparticle, evt, m_recoModuleLabel, 2);

    std::vector<float> longitudinal, transverse, energy;

    for (const art::Ptr<recob::Hit> &collectionHit : collectionViewHits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(collectionHit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;

        const TVector3 spacepointPositionTV = TVector3(spacePoints.at(0)->XYZ()[0], spacePoints.at(0)->XYZ()[1], spacePoints.at(0)->XYZ()[2]);
        const TVector3 displacementTV = spacepointPositionTV - fittingAxisStart;

        // Longitudinal
        const float l = fittingAxisDirection.Dot(displacementTV);

        // Transverse
        const float t = fittingAxisDirection.Cross(displacementTV).Mag();

        // Energy
        const float hitEnergy = ObtainHitEnergy(evt, collectionHit);

        longitudinal.push_back(l);
        transverse.push_back(t);
        energy.push_back(hitEnergy);
    }

    return EnergyProfiles(longitudinal, transverse, energy);
}

/////////////////////////////////////////////////////////////

ProfileManager::EnergyProfiles ProfileManager::CreateProfilesFromTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return EnergyProfiles(std::vector<float>(), std::vector<float>(), std::vector<float>());

    art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    // Work out valid points...

    std::vector<recob::tracking::TrajectoryPoint_t> validTrajectoryPoints;

    for (unsigned int i = 0; i < track->NumberTrajectoryPoints(); ++i)
    {
        const float validDeltaX(std::fabs(track->TrajectoryPoint(i).position.X() - (-999.f)));
        const float validDeltaY(std::fabs(track->TrajectoryPoint(i).position.Y() - (-999.f)));
        const float validDeltaZ(std::fabs(track->TrajectoryPoint(i).position.Z() - (-999.f)));

        if ((validDeltaX < std::numeric_limits<float>::epsilon()) && (validDeltaY < std::numeric_limits<float>::epsilon()) && (validDeltaZ < std::numeric_limits<float>::epsilon()))
            continue;

        validTrajectoryPoints.push_back(track->TrajectoryPoint(i));
    }

    if (validTrajectoryPoints.empty())
        return EnergyProfiles(std::vector<float>(), std::vector<float>(), std::vector<float>());

    // Work out cumulativeL between trajpoints
    float cumulativeL = 0.f;
    std::vector<float> cumulativeLVector = std::vector<float>({0.f});

    for (unsigned int i = 0; i < (validTrajectoryPoints.size() - 1); ++i)
    {
        TVector3 firstPoint = TVector3(validTrajectoryPoints.at(i).position.X(), validTrajectoryPoints.at(i).position.Y(), validTrajectoryPoints.at(i).position.Z());
        TVector3 secondPoint = TVector3(validTrajectoryPoints.at(i + 1).position.X(), validTrajectoryPoints.at(i + 1).position.Y(), validTrajectoryPoints.at(i + 1).position.Z());
        cumulativeL += (secondPoint - firstPoint).Mag();
        cumulativeLVector.push_back(cumulativeL);
    }

    // This is going to be so long...
    std::vector<float> longitudinal, transverse, energy;
    const std::vector<art::Ptr<recob::Hit>> &collectionViewHits = dune_ana::DUNEAnaPFParticleUtils::GetViewHits(pfparticle, evt, m_recoModuleLabel, 2);

    for (const art::Ptr<recob::Hit> &collectionHit : collectionViewHits)
    {
        // Make sure 2D hit has an associated space point
        std::vector<art::Ptr<recob::SpacePoint>> spacePoints = dune_ana::DUNEAnaHitUtils::GetSpacePoints(collectionHit, evt, m_hitModuleLabel, m_recoModuleLabel);

        if (spacePoints.empty())
            continue;

        const art::Ptr<recob::SpacePoint> spacePoint = spacePoints.at(0);
        const TVector3 spacepointPosition = TVector3(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]);

        for (unsigned int i = 0; i < (validTrajectoryPoints.size() - 1); ++i)
        {
            TVector3 firstPoint = TVector3(validTrajectoryPoints.at(i).position.X(), validTrajectoryPoints.at(i).position.Y(), validTrajectoryPoints.at(i).position.Z());
            TVector3 secondPoint = TVector3(validTrajectoryPoints.at(i + 1).position.X(), validTrajectoryPoints.at(i + 1).position.Y(), validTrajectoryPoints.at(i + 1).position.Z());

            if ((secondPoint - firstPoint).Mag() < std::numeric_limits<float>::epsilon())
                continue;

            TVector3 localAxis = (secondPoint - firstPoint).Unit();
            TVector3 displacement = (spacepointPosition - firstPoint);
            float localL = localAxis.Dot(displacement);
            float segmentLength = cumulativeLVector.at(i + 1) - cumulativeLVector.at(i);

            if ((i + 1) != (validTrajectoryPoints.size() - 1))
            {
                if ((localL < 0.f) || (localL > segmentLength))
                    continue;
            }

            // Found segment!
            float globalL = localL + cumulativeLVector.at(i);
            float globalT = localAxis.Cross(displacement).Mag();
            float hitEnergy = ObtainHitEnergy(evt, collectionHit);
            longitudinal.push_back(globalL);
            transverse.push_back(globalT);
            energy.push_back(hitEnergy);

            break;
        }
    }

    return EnergyProfiles(longitudinal, transverse, energy);
}

/////////////////////////////////////////////////////////////

void ProfileManager::FillFitPositions(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    const TVector3 &fittingAxisDirection, const TVector3 &fittingAxisStart, std::vector<std::vector<float>> &fitPositions) const
{
    // Track score
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return;

    const float trackScore = metaMap.at("TrackScore");

    // Now find appropriate length
    float length = -1.f;

    if (dune_ana::DUNEAnaPFParticleUtils::IsShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel))
    {
        art::Ptr<recob::Shower> shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(pfparticle, evt, m_recoModuleLabel, m_showerModuleLabel);

        const float showerLength = shower->Length();

        if (showerLength > length)
            length = showerLength;
    }

    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
    {
        art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

        const float trackLength = track->Length();

        if ((trackScore > 0.5f) || (length < 0.f))
            length = trackLength;
    }

    const TVector3 endPosition = fittingAxisStart + (length * fittingAxisDirection);
    fitPositions.push_back(std::vector<float>({static_cast<float>(fittingAxisStart.X()), static_cast<float>(fittingAxisStart.Y()), static_cast<float>(fittingAxisStart.Z())}));
    fitPositions.push_back(std::vector<float>({static_cast<float>(endPosition.X()), static_cast<float>(endPosition.Y()), static_cast<float>(endPosition.Z())}));
}

/////////////////////////////////////////////////////////////

void ProfileManager::FillFitPositionsFromTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
    std::vector<std::vector<float>> &fitPositions) const
{
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel))
        return;

    art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfparticle, evt, m_recoModuleLabel, m_trackModuleLabel);

    for (unsigned int i = 0; i < track->NumberTrajectoryPoints(); ++i)
    {
        const float validDeltaX(std::fabs(track->TrajectoryPoint(i).position.X() - (-999.f)));
        const float validDeltaY(std::fabs(track->TrajectoryPoint(i).position.Y() - (-999.f)));
        const float validDeltaZ(std::fabs(track->TrajectoryPoint(i).position.Z() - (-999.f)));

        if ((validDeltaX < std::numeric_limits<float>::epsilon()) && (validDeltaY < std::numeric_limits<float>::epsilon()) && (validDeltaZ < std::numeric_limits<float>::epsilon()))
            continue;

        fitPositions.push_back(std::vector<float>({static_cast<float>(track->TrajectoryPoint(i).position.X()), 
            static_cast<float>(track->TrajectoryPoint(i).position.Y()), static_cast<float>(track->TrajectoryPoint(i).position.Z())}));
    }
}

/////////////////////////////////////////////////////////////

float ProfileManager::ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    const double charge = dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, {hit});
    const double nElectrons = m_calorimetryAlg.ElectronsFromADCArea(charge, hit->WireID().Plane);
    const double hitEnergy = nElectrons / m_recombFactor / util::kGeVToElectrons;

    return hitEnergy;
}

/////////////////////////////////////////////////////////////

}
