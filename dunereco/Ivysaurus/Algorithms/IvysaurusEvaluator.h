////////////////////////////////////////////////////////////////////////
/// \file    IvysaurusEvaluator.cxx
/// \brief   Functions to use the Ivysaurus network
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef IVYSAURUSEVALUATOR_H
#define IVYSAURUSEVALUATOR_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"

#include "dunereco/Ivysaurus/TensorFlow/IvysaurusGraph.h"

namespace ivysaurus
{

class IvysaurusEvaluator
{
  public:
      IvysaurusEvaluator(const fhicl::ParameterSet& pset);

  private:
      void IvysaurusUseEvaluate(const art::Ptr<recob::PFParticle> pfparticle, const art::Event &evt) const;

      std::string m_networkDir;
      std::string m_networkName;
};
}

#endif  // IVYSAURUSEVALUATOR_H
