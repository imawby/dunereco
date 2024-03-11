////////////////////////////////////////////////////////////////////////
/// \file    ProfileManager.cxx
/// \brief   A class to manage the Ivysaurus 2D grid input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef PROFILEMANAGER_H
#define PROFILEMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace ivysaurus
{

class ProfileManager
{
  public:
    class EnergyProfiles
    {
        public:
            EnergyProfiles(const std::vector<float> &longitudinal, const std::vector<float> &transverse, const std::vector<float> &energy);

            //private:
            std::vector<float> m_longitudinal;
            std::vector<float> m_transverse;
            std::vector<float> m_energy;
    };

    ProfileManager(const fhicl::ParameterSet& pset);
    ~ProfileManager();

    void GetProfiles(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, const int directionMode, 
        std::vector<float> &longitudinalProfile, std::vector<float> &transverseProfile, std::vector<float> &energies, 
        std::vector<std::vector<float>> &fitPositions) const;

    bool GetTrueDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        TVector3 &trueDirection, TVector3 &truePosition) const;

    bool GetRecoStartDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        TVector3 &recoStartDirection, TVector3 &recoPosition) const;

    bool GetPCADirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        TVector3 &pcaDirection, TVector3 &pcaPosition) const;

    bool GetTrackStubDirectionAndPosition(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        TVector3 &trackStubDirection, TVector3 &trackStubPosition) const;

    EnergyProfiles CreateProfiles(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        const TVector3 &fittingAxis, const TVector3 &fittingAxisStart) const;

    EnergyProfiles CreateProfilesFromTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle) const;

    void FillFitPositions(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        const TVector3 &fittingAxisDirection, const TVector3 &fittingAxisStart, std::vector<std::vector<float>> &fitPositions) const;

    void FillFitPositionsFromTrack(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle,
        std::vector<std::vector<float>> &fitPositions) const;

    float ObtainHitEnergy(const art::Event &evt, const art::Ptr<recob::Hit> &hit) const;

    std::string m_hitModuleLabel;
    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;
    std::string m_showerModuleLabel;

    float m_recombFactor;
    calo::CalorimetryAlg m_calorimetryAlg;
};

}

#endif  // PROFILEMANAGER_H
