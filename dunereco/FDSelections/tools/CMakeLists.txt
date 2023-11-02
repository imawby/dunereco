simple_plugin(LongestRecoVertexTrackSelector "tool"
                      dunereco_AnaUtils
                      lardataobj::RecoBase
                      larpandora::LArPandoraInterface
                      art::Persistency_Common 
                      art::Persistency_Provenance 

                       
             )

simple_plugin(HighestPandrizzleScoreRecoVertexShowerSelector "tool"
                      lardataobj::RecoBase
                      larpandora::LArPandoraInterface
                      art::Persistency_Common 
                      art::Persistency_Provenance 

                       

                      PandrizzleAlg    
                      #dune_FDSelections_pandizzle
             )


simple_plugin(HighestPandizzleScoreRecoVertexTrackSelector "tool"
                      lardataobj::RecoBase
                      larpandora::LArPandoraInterface
                      art::Persistency_Common 
                      art::Persistency_Provenance 

                       
                      PandizzleAlg
             )

simple_plugin(HighestDeepPanScoreRecoVertexTrackSelector "tool"
    	  	        lardataobj_RecoBase
    	  	        larpandora_LArPandoraInterface
                ${ROOT_BASIC_LIB_LIST}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_CORE}
		       	   ${ART_FRAMEWORK_SERVICES_REGISTRY}
 	  	           art_Persistency_Common canvas
 	  	           art_Persistency_Provenance canvas
 	  	           art_Utilities canvas
                   dunereco_TrackPID_algorithms
             )

simple_plugin(HighestEnergyRecoVertexShowerSelector "tool"
                      lardataobj::RecoBase
                      larcorealg::Geometry
                      larcore::Geometry_Geometry_service
                      larreco::RecoAlg
                      larpandora::LArPandoraInterface
                      dunereco_AnaUtils
                      art::Persistency_Common 
                      art::Persistency_Provenance 

                       
             )
