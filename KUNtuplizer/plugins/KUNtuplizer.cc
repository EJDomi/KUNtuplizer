// -*- C++ -*-
//
// Package:    Framework/KUNtuplizer
// Class:      KUNtuplizer
// 
/**\class KUNtuplizer KUNtuplizer.cc Framework/KUNtuplizer/plugins/KUNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sadia Khalil
//         Created:  Fri, 21 Sep 2018 22:01:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

//#include "Geometry/Records/interface/MuonGeometryRecord.h"
//#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
//#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "Framework/KUNtuplizer/interface/EventInfoTree.h"
#include "Framework/KUNtuplizer/interface/GenInfoTree.h"
//#include "Framework/KUNtuplizer/interface/METTree.h"
//#include "Framework/KUNtuplizer/interface/ElectronTree.h"
//#include "Framework/KUNtuplizer/interface/MuonTree.h"
//#include "Framework/KUNtuplizer/interface/JetTree.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
//
// class declaration
//

class KUNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit KUNtuplizer(const edm::ParameterSet&);
      ~KUNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginRun(edm::Run const&, edm::EventSetup const&); 
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&); 
      virtual void endJob() override;

      // ----------member data ---------------------------  
   edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
   edm::EDGetTokenT<GenEventInfoProduct>         genToken_;
   edm::EDGetTokenT<LHEEventProduct>             genlheToken_;
   edm::InputTag                                 puInfo_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>   vtxToken_;

   edm::Service<TFileService> fs_;
   TTree* tree_;    
   GenInfoTree   genevt_; 
   EventInfoTree evt_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KUNtuplizer::KUNtuplizer(const edm::ParameterSet& iConfig):
   genPartsToken_  (consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genToken_       (consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   genlheToken_    (consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
   puInfo_         (iConfig.getParameter<edm::InputTag>("puInfo")),
   vtxToken_ (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices")))
   
{
   consumes<std::vector<PileupSummaryInfo>>(puInfo_);
   //now do what ever initialization is needed
   usesResource("TFileService");
   
}


KUNtuplizer::~KUNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
KUNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm; 
   genevt_.clearTreeVectors(); 
   evt_.clearTreeVectors(); 

   // Basic event info
   evt_.runno = iEvent.eventAuxiliary().run();
   evt_.lumisec = iEvent.eventAuxiliary().luminosityBlock();
   evt_.evtno = iEvent.eventAuxiliary().event();
   evt_.countEvents = 1.;//test

   //Pileup
   edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
   iEvent.getByLabel(puInfo_, puInfo);
      
   edm::Handle<std::vector<pat::PackedGenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   std::cout << "???" << std::endl; 
   tree_->Fill();
}

// ------------ method called once each run ----------------
void
KUNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
   //edm::ESHandle<ME0Geometry> hGeom;
   //iSetup.get<MuonGeometryRecord>().get(hGeom);
   //ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
void
KUNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just before starting event loop  ------------
void 
KUNtuplizer::beginJob()
{
   tree_ = fs_->make<TTree>("kutree", "kutree") ;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KUNtuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KUNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KUNtuplizer);
