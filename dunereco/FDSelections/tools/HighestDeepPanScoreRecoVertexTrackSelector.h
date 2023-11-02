#ifndef HIGHESTDEEPPANSCOREVERTEXTRACKSELECTOR_H_SEEN
#define HIGHESTDEEPPANSCOREVERTEXTRACKSELECTOR_H_SEEN

//STL
#include <iostream>
#include <limits>

//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"

//LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"

//DUNE
#include "dunereco/TrackPID/algorithms/CTPHelper.h"
#include "dunereco/TrackPID/products/CTPResult.h"

//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class HighestDeepPanScoreRecoVertexTrackSelector : public RecoTrackSelector{
    public:
      explicit HighestDeepPanScoreRecoVertexTrackSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;

      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;

      ctp::CTPHelper fDeepPanHelper;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestDeepPanScoreRecoVertexTrackSelector)
#endif
