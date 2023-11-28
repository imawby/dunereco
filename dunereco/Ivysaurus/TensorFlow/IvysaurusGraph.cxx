////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       IvysaurusGraph
// Authors:     R.Sulej (Robert.Sulej@cern.ch), from DUNE, FNAL/NCBJ, Sept. 2017
//              P.Plonski,                      from DUNE, WUT, Sept. 2017
//              S.Alonso-Monsalve,              from DUNE, CERN, Aug. 2018
// Iterface to run Tensorflow graph saved to a file. First attempts, quite functional.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"

#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/core/framework/logging.h" 

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

#include "dunereco/Ivysaurus/Managers/GridManager.h"
#include "dunereco/Ivysaurus/Managers/TrackVarManager.h"
#include "dunereco/Ivysaurus/Utils/IvysaurusUtils.h"

#include "IvysaurusGraph.h"

/////////////////////////////////////////////////////////////

ivysaurus::IvysaurusGraph::IvysaurusGraph(fhicl::ParameterSet const &pset) :
    m_networkDirectory(pset.get<std::string>("NetworkDirectory")), 
    m_gridManager(pset.get<fhicl::ParameterSet>("GridManager")),
    m_trackVarManager(pset.get<fhicl::ParameterSet>("TrackVarManager")),
    m_recoModuleLabel(pset.get<std::string>("RecoModuleLabel")),
    m_nTrackVars(pset.get<int>("NTrackVars"))
{
    std::cout << "Starting to build the graph" << std::endl;

    // Create dummy options.
    tensorflow::SessionOptions sessionOptions;
    tensorflow::RunOptions runOptions;

    std::cout << "Session started... reading network architecture" << std::endl;

    // Load the model bundle. (this returns a status code)
    const auto loadResult = tensorflow::LoadSavedModel(sessionOptions, runOptions, m_networkDirectory, 
        {tensorflow::kSavedModelTagServe}, &m_savedModelBundle);

    // Check if loading was okay.
    TF_CHECK_OK(loadResult);
}

/////////////////////////////////////////////////////////////

void ivysaurus::IvysaurusGraph::IvysaurusUseEvaluate(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    // Obtain the input grid tensors
    std::cout << "Making the grid tensors..." << std::endl;
    tensorflow::Tensor startTensorU = ObtainInputGridTensor(evt, pfparticle, true, IvysaurusUtils::PandoraView::TPC_VIEW_U);
    tensorflow::Tensor endTensorU = ObtainInputGridTensor(evt, pfparticle, false, IvysaurusUtils::PandoraView::TPC_VIEW_U);

    tensorflow::Tensor startTensorV = ObtainInputGridTensor(evt, pfparticle, true, IvysaurusUtils::PandoraView::TPC_VIEW_V);
    tensorflow::Tensor endTensorV = ObtainInputGridTensor(evt, pfparticle, false, IvysaurusUtils::PandoraView::TPC_VIEW_V);

    tensorflow::Tensor startTensorW = ObtainInputGridTensor(evt, pfparticle, true, IvysaurusUtils::PandoraView::TPC_VIEW_W);
    tensorflow::Tensor endTensorW = ObtainInputGridTensor(evt, pfparticle, false, IvysaurusUtils::PandoraView::TPC_VIEW_W);

    // Obtain the input track variable tensors
    std::cout << "Making the track variable tensors..." << std::endl;
    tensorflow::Tensor trackVarTensor = ObtainInputTrackTensor(evt, pfparticle);

    auto trackVarTensorMap = trackVarTensor.tensor<float, 2>();
    for (int i =0; i < 10; ++i)
        std::cout << "trackVarTensorMap(0, i): " << trackVarTensorMap(0, i) << std::endl;


    // Link the data with some tags so tensorflow know where to put those data entries.
    std::cout << "Hooking up input/output layers..." << std::endl;
    std::vector<std::pair<std::string, tensorflow::Tensor>> feedInputs = {{"serving_default_input_1:0", startTensorU}, {"serving_default_input_2:0", endTensorU},
                                                                          {"serving_default_input_3:0", startTensorV}, {"serving_default_input_4:0", endTensorV},  
                                                                          {"serving_default_input_5:0", startTensorW}, {"serving_default_input_6:0", endTensorW},
                                                                          {"serving_default_input_7:0", trackVarTensor}};
    std::vector<std::string> fetches = { "StatefulPartitionedCall:0" };

    // We need to store the results somewhere.
    std::vector<tensorflow::Tensor> outputs;

    // Let's run the model...
    std::cout << "Running the model..." << std::endl;
    auto status = m_savedModelBundle.GetSession()->Run(feedInputs, fetches, {}, &outputs);

    // ... and print out it's predictions.
    std::cout << "What does the model predict?..." << std::endl;
    auto output_map = outputs[0].tensor<float, 2>();

    IvysaurusScores ivysaurusScores;

    ivysaurusScores.m_muonScore = output_map(0, 0);
    ivysaurusScores.m_protonScore = output_map(0, 1);
    ivysaurusScores.m_pionScore = output_map(0, 2);
    ivysaurusScores.m_electronScore = output_map(0, 3);
    ivysaurusScores.m_photonScore = output_map(0, 4);
    ivysaurusScores.m_otherScore = output_map(0, 5);

    float highestScore = -std::numeric_limits<float>::max();
    int count = 0;

    for (float score : {ivysaurusScores.m_muonScore, ivysaurusScores.m_protonScore, ivysaurusScores.m_pionScore, 
        ivysaurusScores.m_electronScore, ivysaurusScores.m_photonScore, ivysaurusScores.m_otherScore})
    {
        if (score > highestScore)
        {
            highestScore = score;
            ivysaurusScores.m_particleType = count;
        }

        ++count;
    }

    std::cout << "muonScore: " << ivysaurusScores.m_muonScore << std::endl;
    std::cout << "protonScore: " << ivysaurusScores.m_protonScore << std::endl;
    std::cout << "pionScore: " << ivysaurusScores.m_pionScore << std::endl;
    std::cout << "electronScore: " << ivysaurusScores.m_electronScore << std::endl;
    std::cout << "photonScore: " << ivysaurusScores.m_photonScore << std::endl;
    std::cout << "otherScore: " << ivysaurusScores.m_otherScore << std::endl;
    std::cout << "particleType: " << ivysaurusScores.m_particleType << std::endl;

    TF_CHECK_OK(status);
}

/////////////////////////////////////////////////////////////

tensorflow::Tensor ivysaurus::IvysaurusGraph::ObtainInputGridTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, 
    const bool isStart, const IvysaurusUtils::PandoraView &pandoraView)
{
    GridManager::Grid grid = m_gridManager.ObtainViewGrid(evt, pfparticle, pandoraView, isStart);
    m_gridManager.FillViewGrid(evt, pfparticle, grid);

    tensorflow::Tensor gridTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, grid.GetAxisDimensions(), grid.GetAxisDimensions(), 1 }));

    // Fill the tensors
    auto tensorMap = gridTensor.tensor<float, 4>();

    for (unsigned int driftIndex = 0; driftIndex < grid.GetAxisDimensions(); ++driftIndex)
        for (unsigned int wireIndex = 0; wireIndex < grid.GetAxisDimensions(); ++wireIndex)
            tensorMap(0, driftIndex, wireIndex, 0) = grid.GetGridValues().at(driftIndex).at(wireIndex);

    return gridTensor;
}

/////////////////////////////////////////////////////////////

tensorflow::Tensor ivysaurus::IvysaurusGraph::ObtainInputTrackTensor(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    // Order should be: 
    // nTrackChildren, nShowerChildren, nGrandChildren, nChildHits, childEnergy, childTrackScore, trackLength, trackWobble, trackScore, momComparison  

    tensorflow::Tensor trackVarTensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, m_nTrackVars}));
    auto trackVarsTensorMap = trackVarTensor.tensor<float, 2>();

    TrackVarManager::TrackVars trackVars;

    m_trackVarManager.EvaluateTrackVars(evt, pfparticle, trackVars);
    m_trackVarManager.NormaliseTrackVars(trackVars);

    trackVarsTensorMap(0, 0) = trackVars.GetNTrackChildren();
    trackVarsTensorMap(0, 1) = trackVars.GetNShowerChildren();
    trackVarsTensorMap(0, 2) = trackVars.GetNGrandChildren();
    trackVarsTensorMap(0, 3) = trackVars.GetNChildHits();
    trackVarsTensorMap(0, 4) = trackVars.GetChildEnergy();
    trackVarsTensorMap(0, 5) = trackVars.GetChildTrackScore();
    trackVarsTensorMap(0, 6) = trackVars.GetTrackLength();
    trackVarsTensorMap(0, 7) = trackVars.GetWobble();
    trackVarsTensorMap(0, 8) = GetTrackShowerScore(evt, pfparticle);
    trackVarsTensorMap(0, 9) = trackVars.GetMomentumComparison();

    return trackVarTensor;
}

/////////////////////////////////////////////////////////////

double ivysaurus::IvysaurusGraph::GetTrackShowerScore(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata = dune_ana::DUNEAnaPFParticleUtils::GetMetadata(pfparticle, evt, m_recoModuleLabel);
    const auto metaMap = metadata->GetPropertiesMap();

    if (metaMap.find("TrackScore") == metaMap.end())
        return -1.0;

    const double trackScore = metaMap.at("TrackScore");

    return trackScore;
}

/////////////////////////////////////////////////////////////

// -------------------------------------------------------------------

/*
    // Trying to work out the input/output layers
    const auto signatures = m_savedModelBundle.GetSignatures();
    for (auto &entry : signatures)
    {
        std::cout << "entry.first: " << entry.first << std::endl;

        auto inputsMap = entry.second.inputs();

        std::cout << "inputsMap.size(): " << inputsMap.size() << std::endl;

        for (auto &inputEntry : inputsMap)
        {
            std::cout << "inputEntry.first: " << inputEntry.first << std::endl;
            std::cout << "inputEntry.second.getName(): " << inputEntry.second.name() << std::endl;
        }

        auto outputsMap = entry.second.outputs();

        std::cout << "outputsMap.size(): " << outputsMap.size() << std::endl;

        for (auto &outputEntry : outputsMap)
        {
            std::cout << "outputEntry.first: " << outputEntry.first << std::endl;
            std::cout << "outputEntry.second.getName(): " << outputEntry.second.name() << std::endl;
        }
    }
*/
