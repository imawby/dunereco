#include "IvysaurusUtils.h"

namespace IvysaurusUtils
{

double CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits, 
  const std::vector<art::Ptr<recob::Hit>> &eventHits, const int trackID)
{
    int nMatchesInSelHits = 0;
    int nMatchesInAllHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1); 
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);;
        if (matchedID == trackID) nMatchesInAllHits++;
    }

    const double completeness = (nMatchesInAllHits > 0) ? static_cast<double>(nMatchesInSelHits) / static_cast<double>(nMatchesInAllHits) : 0.0;

    return completeness;
}

/////////////////////////////////////////////////////////////////////////

double HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits,
    const int trackID)
{
    int nMatchesInSelHits = 0;

    for (art::Ptr<recob::Hit> hit : selectedHits)
    {
        const int matchedID = TruthMatchUtils::TrueParticleID(clockData, hit, 1);
        if (matchedID == trackID) nMatchesInSelHits++;
    }

    const double purity = (selectedHits.size() > 0) ? static_cast<double>(nMatchesInSelHits) / selectedHits.size() : 0.0;

    return purity;
}
}




