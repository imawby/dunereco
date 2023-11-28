////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       IvysaurusGraph
//// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
///               P.Plonski,                      from DUNE, WUT, Sept. 2017
////              S. Alonso Monsalve,             from DUNE, CERN, Aug. 2018
//// Iterface to run Tensorflow graph saved to a file. First attempts, almost functional.
////
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef IVYSAURUSGRAPH_H
#define IVYSAURUSGRAPH_H

#include <memory>
#include <vector>
#include <string>

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h" 

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "dunereco/Ivysaurus/Managers/GridManager.h"
#include "dunereco/Ivysaurus/Managers/TrackVarManager.h"

namespace ivysaurus
{

class IvysaurusGraph
{
  struct IvysaurusScores
  {
      float m_muonScore;
      float m_protonScore;
      float m_pionScore;
      float m_electronScore;
      float m_photonScore;
      float m_otherScore;
      int m_particleType;

      IvysaurusScores() : m_muonScore(-1.f), m_protonScore(-1.f), m_pionScore(-1.f), m_electronScore(-1.f), m_photonScore(-1.f), m_otherScore(-1.f), m_particleType(-1) {};
  };

public:
    IvysaurusGraph(fhicl::ParameterSet const &pset);

    void IvysaurusUseEvaluate(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

private:
    tensorflow::Tensor ObtainInputGridTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
        const bool isStart, const IvysaurusUtils::PandoraView &pandoraView);

    tensorflow::Tensor ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    double GetTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle);

    std::string m_networkDirectory;
    GridManager m_gridManager;
    TrackVarManager m_trackVarManager;
    tensorflow::SavedModelBundleLite m_savedModelBundle;
    std::string m_recoModuleLabel;
    int m_nTrackVars;

};

} // namespace ivysaurus

#endif // IVYSAURUSGRAPH_H
