#include "HighestDeepPanScoreRecoVertexTrackSelector.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/TrackPID/algorithms/CTPHelper.h"
#include "dunereco/TrackPID/products/CTPResult.h"

FDSelectionTools::HighestDeepPanScoreRecoVertexTrackSelector::HighestDeepPanScoreRecoVertexTrackSelector(fhicl::ParameterSet const& pset) :
    fTrackModuleLabel(pset.get<std::string>("ModuleLabels.TrackModuleLabel")),
    fPFParticleModuleLabel(pset.get<std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fDeepPanHelper(pset.get<fhicl::ParameterSet>("DeepPanHelper"))
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Track> FDSelectionTools::HighestDeepPanScoreRecoVertexTrackSelector::SelectTrack(art::Event const & evt)
{
  art::Ptr<recob::Track> selTrack;

  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return selTrack;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  // Loop over primaries, get the associated track and then find that with the highest pandizzle score
  double highestDeepPanScore(std::numeric_limits<double>::lowest());

  for (art::Ptr<recob::PFParticle> childPFP : nuChildren) 
  {
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel))
      continue;

    art::Ptr<recob::Track> childTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel);

    // DeepPan variables
    ctp::CTPResult deepPanResult = fDeepPanHelper.RunConvolutionalTrackPID(childPFP, evt);

    if (deepPanResult.IsValid())
    {
        std::cout << "WE HAVE A VALID DEEPPAN RESULT!" << std::endl;

      float muonScore = deepPanResult.GetMuonScore();
      float pionScore = deepPanResult.GetPionScore();
      float protonScore = deepPanResult.GetProtonScore();

      std::cout << "muonScore: " << muonScore << std::endl;
      std::cout << "pionScore: " << pionScore << std::endl;
      std::cout << "protonScore: " << protonScore << std::endl;

      if (muonScore > highestDeepPanScore)
      {
        highestDeepPanScore = muonScore;
        selTrack = childTrack;
      }
    }
    else { std::cout << "WE DO NOT HAVE A VALID DEEPPAN RESULT!" << std::endl; }

  }

  return selTrack;
}
