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
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/VertexTypes.h"
#include "DataFormats/BTauReco/interface/ParticleMasses.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

//#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include "Framework/KUNtuplizer/interface/EventInfoTree.h"
#include "Framework/KUNtuplizer/interface/GenInfoTree.h"
#include "Framework/KUNtuplizer/interface/ElectronTree.h"
//#include "Framework/KUNtuplizer/interface/MuonTree.h"
//#include "Framework/KUNtuplizer/interface/METTree.h"
//#include "Framework/KUNtuplizer/interface/JetTree.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include "Framework/KUNtuplizer/interface/SVDiscrTool.hh"
#include "Framework/KUNtuplizer/src/SVDiscrTool.cc"
//
// class declaration
//

class KUNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit KUNtuplizer(const edm::ParameterSet&);
      ~KUNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   protected:
      float catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound);
      const float& catchInfs (const float& in,const float& replace_value);     

   private:
      virtual void beginRun(edm::Run const&, edm::EventSetup const&); 
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&); 
      virtual void endJob() override;

      // ----------member data ---------------------------  
   //edm::EDGetTokenT<std::vector<pat::PackedGenParticle>>         genPartsToken_;
   edm::EDGetTokenT<std::vector<reco::GenParticle>>              genPartsToken_;
   edm::EDGetTokenT<GenEventInfoProduct>                         genToken_;
   edm::EDGetTokenT<LHEEventProduct>                             genlheToken_;
   edm::InputTag                                                 puInfo_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>                   vtxToken_;
   edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> svToken_;
   edm::EDGetTokenT<std::vector<pat::Electron>>                  elecsToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection>              packedPFCands_;
   edm::EDGetTokenT<edm::View<pat::Jet> >                        ak4jetsToken_;
   edm::EDGetTokenT< reco::ConversionCollection >                conv_;
   edm::EDGetTokenT< reco::BeamSpot >                            beamSpot_;
   edm::EDGetTokenT< double >                                    rho_;
   edm::Service<TFileService>                                    fs_;
   std::map<std::string, TH1D*>                                  h1_; 
   TTree*        tree_;    
   GenInfoTree   genevt_; 
   EventInfoTree evt_;
   ElectronTree  ele_;
   int           ndofPV_;
   double        zPV_;
   double        rhoPV_;
   double        ptMaxSV_;
   double        dlenSigSV_;
   double        d2dSV_;
   double        cosPdotVSV_;
   int           ndauSV_;
   double        jetdRSV_;
   //double        ak4ptmax_;
   double        genPartdRSV_;
   std::string   sv_model_file_;
   SVDiscrTool   m_SVDiscrTool;
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
   //genPartsToken_  (consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genPartsToken_  (consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genToken_       (consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   genlheToken_    (consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
   puInfo_         (iConfig.getParameter<edm::InputTag>("puInfo")),
   vtxToken_       (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   svToken_        (consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
   elecsToken_     (consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
   packedPFCands_  (consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCands"))), 
   ak4jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
   conv_           (consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversion"))),
   beamSpot_       (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
   rho_            (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
   ndofPV_         (iConfig.getParameter<int>("ndofPV")),
   zPV_            (iConfig.getParameter<double>("zPV")),
   rhoPV_          (iConfig.getParameter<double>("rhoPV")),
   ptMaxSV_        (iConfig.getParameter<double>("ptMaxSV")), 
   dlenSigSV_      (iConfig.getParameter<double>("dlenSigSV")),
   d2dSV_          (iConfig.getParameter<double>("d2dSV")),
   cosPdotVSV_     (iConfig.getParameter<double>("cosPdotVSV")),
   ndauSV_         (iConfig.getParameter<int>("ndauSV")),
   jetdRSV_        (iConfig.getParameter<double>("jetdRSV")),
   //ak4ptmax_       (iConfig.getParameter<double>("ak4ptmax")),
   genPartdRSV_    (iConfig.getParameter<double>("genPartdRSV")),
   sv_model_file_  (iConfig.getParameter<std::string>("sv_model_file"))
                    
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

   //Generator weights
   float genwt = 1.0;
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken( genToken_,genInfo);
   
   if (genInfo.isValid()){
      genwt = genInfo->weight();
      genwt /= std::abs(genwt);   
   }
   genevt_.genWt = genwt;
   
   //LHE weight
   edm::Handle<LHEEventProduct> lheInfo;
   iEvent.getByToken(genlheToken_, lheInfo);
   
   if(lheInfo.isValid()) {
      int size = lheInfo->weights().size();
      double asdd = 1.;
      if (size != 0){
         //asdd = lheInfo->weights()[0].wgt;
         asdd = lheInfo->originalXWGTUP();
      }
      
      genevt_.lheWtIDs.reserve(size);
      genevt_.lheWts.reserve(size);
      for(unsigned int i=0; i<lheInfo->weights().size(); ++i) {         
         double asdde =lheInfo->weights()[i].wgt;
         int asddeID = std::stoi(lheInfo->weights()[i].id);
         genevt_.lheWtIDs.push_back(asddeID);
         genevt_.lheWts.push_back(asdde/asdd);
      }
   }

   // Basic event info
   evt_.runno = iEvent.eventAuxiliary().run();
   evt_.lumisec = iEvent.eventAuxiliary().luminosityBlock();
   evt_.evtno = iEvent.eventAuxiliary().event();
   //evt_.countEvents = 1.;//test

   //Pileup
   edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
   iEvent.getByLabel(puInfo_, puInfo);
   std::vector<PileupSummaryInfo>::const_iterator pvi;
   for(pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
      evt_.npuTrue = pvi->getTrueNumInteractions(); 
      evt_.npuInt = pvi->getBunchCrossing(); 
      evt_.puBX = pvi->getPU_NumInteractions();
   }
   
   // Vertex
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);

   std::vector<reco::Vertex> goodvtx;
   goodvtx.clear();

   if (vertices->empty()) return;
   evt_.nGoodVtx = 0;
   evt_.npv = vertices->size();
   
   int prVtx = -1;
   // store the vertex related variables after basic PV selection
   for (size_t i = 0; i < vertices->size(); i++) {
      if (vertices->at(i).isFake()) continue;
      if (vertices->at(i).ndof() <= ndofPV_) continue;
      if (fabs(vertices->at(i).position().z()) > zPV_) continue;
      if (vertices->at(i).position().rho() > rhoPV_) continue;    
      if (prVtx < 0) {prVtx = i; }
      evt_.ndofVtx.push_back(vertices->at(i).ndof());
      evt_.chi2Vtx.push_back(vertices->at(i).chi2());
      evt_.zVtx.push_back(vertices->at(i).position().z());
      evt_.rhoVtx.push_back(vertices->at(i).position().rho());
      evt_.etaVtx.push_back(vertices->at(i).position().eta());
      evt_.phiVtx.push_back(vertices->at(i).position().phi());
      goodvtx.push_back(vertices->at(i));
      evt_.nGoodVtx++;
   }
   
   if (evt_.nGoodVtx==0) return;
   //auto primaryVertex=vertices->at(prVtx);   
   const auto & pv = goodvtx[0]; 
  
   // Jets
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByToken(ak4jetsToken_, ak4jets);

   // GenParticles
   //edm::Handle<std::vector<pat::PackedGenParticle>> genParts;
   edm::Handle<std::vector<reco::GenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   
   // Secondary Vertices
   edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
   iEvent.getByToken(svToken_, secVertices);
   evt_.nSV = 0; evt_.nGoodSV = 0;
   evt_.nSV = secVertices->size();
    
   VertexDistance3D vdist;  
   VertexDistanceXY vdistXY;
  
   edm::ESHandle<TransientTrackBuilder> builder_;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);


 
   for (const auto & sv : *secVertices) {

 
      //sum pt of all tracks in SV
      if(sv.pt() > ptMaxSV_) continue;//20

     //SV direction vector, SV flight direction w.r.t PV 
      math::XYZVector svDirU = sv.momentum().Unit();
      TVector3 svDir3(svDirU.x(), svDirU.y(), svDirU.z());      
      GlobalVector svDir(sv.px(), sv.py(), sv.pz());
      GlobalVector svDirPV(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(), sv.vertex().z() - pv.z()); 


      //store the branches after gen particle matching
      
      float genSVdR = -1000.0; float matchGenSVdR = -1000.0;
      int matchGenMomID = 0; bool matchGenIsHardProcess = false; bool matchGenIsLastCopyBeforeFSR = false;     
      float matchGenPt = -1000.0; float matchGenEta = -1000.0; float matchGenMass = -1000.0; 
      float matchGenPhi = -1000.0; float matchGenEnergy = -1000.0;

      int nbs = 0.0, ngs = 0.0, ncs = 0.0, nuds = 0.0, nss = 0.0, nmatchothers = 0, nothers = 0;
      bool isgb = false, isb = false, isc = false;
      bool isg = false, iss = false, isud = false, ismatchother = false, isother = false;

      bool isbMeson = false; bool isbBaryon = false; 
      bool iscMeson = false; bool iscBaryon = false;
      bool issMeson = false; bool issBaryon = false;
      bool isudMeson = false; bool isudBaryon = false;

      bool fill1stHeavyGen = false;
      bool fill1stNotHeavyGen = false;
      float mindR = 100.0; float maxdR = 0.0;

      //float matchGenPt = 0.0, matchGenEta = 0.0, matchGenMass = 0.0, matchGenPhi = 0.0;

      //for (const pat::PackedGenParticle & gp : *genParts) {
      for (const reco::GenParticle & gp : *genParts) {
         genSVdR = reco::deltaR(gp.p4(), sv.p4());
         if (genSVdR <= mindR) { mindR = genSVdR;}
         if (genSVdR > maxdR)  { maxdR = genSVdR;}

         if(genSVdR <= genPartdRSV_){//0.4
            
            int pdgid = abs(gp.pdgId());

            isbMeson = pdgid%1000 > 500 && pdgid%1000 < 600;
            isbBaryon = pdgid > 5000 &&  pdgid < 6000;
            
            iscMeson = pdgid%1000 > 400 && pdgid%1000 < 500;
            iscBaryon = pdgid > 4000 &&  pdgid < 5000;

            issMeson = pdgid%1000 > 300 && pdgid%1000 < 400;
            issBaryon = pdgid > 3000 &&  pdgid < 4000;
 
            isudMeson = pdgid%1000 > 100 && pdgid%1000 < 300;
            isudBaryon = pdgid > 1000 &&  pdgid < 3000;

            if      (pdgid == 5 || isbMeson || isbBaryon){nbs++;}
            else if (pdgid == 4 || iscMeson || iscBaryon){ncs++;}
            else if (pdgid == 21){ngs++;}
            else if (pdgid == 3 || pdgid == 130 || issMeson || issBaryon){nss++;}
            else if (pdgid == 2 || pdgid == 1 || isudMeson || isudBaryon){nuds++;}
            else {nmatchothers++;}

            if (nbs+ncs == 0 && fill1stNotHeavyGen == false){                           
               matchGenSVdR   = genSVdR; 
               matchGenPt     = gp.p4().pt();               
               matchGenEta    = gp.p4().eta();
               matchGenPhi    = gp.p4().phi();
               matchGenMass   = gp.p4().mass();
               matchGenEnergy = gp.p4().energy();
               matchGenMomID               = gp.mother()->pdgId();               
               matchGenIsHardProcess       = gp.isHardProcess();              
               matchGenIsLastCopyBeforeFSR = gp.statusFlags().isLastCopyBeforeFSR();               
               fill1stNotHeavyGen = true;               
            }// if a b or c quark is found, then overwrite the P4 to the leading one
            else if (nbs+ncs != 0 && fill1stHeavyGen == false) {
               matchGenSVdR   = genSVdR; 
               matchGenPt     = gp.p4().pt();
               matchGenEta    = gp.p4().eta();
               matchGenPhi    = gp.p4().phi();
               matchGenMass   = gp.p4().mass();
               matchGenEnergy = gp.p4().energy();
               matchGenMomID               = gp.mother()->pdgId();
               matchGenIsHardProcess       = gp.isHardProcess();
               matchGenIsLastCopyBeforeFSR = gp.statusFlags().isLastCopyBeforeFSR();
               fill1stHeavyGen = true;
         }
        }
        else {nothers++;}
      }

      if (nbs > 0 && ngs > 0) isgb = true;
      else if (nbs > 0) isb = true;
      else if (ncs > 0) isc = true;
      else if (ngs > 0) isg = true;
      else if (nss > 0) iss = true;
      else if (nuds > 0) isud = true;
      else if (nmatchothers > 0) ismatchother = true;
      else isother = true;
       
      //else continue;


      //SV track separation from jets     
      bool isInJetDir = false;      
      for (const pat::Jet & jet : *ak4jets ){
         GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(), sv.vertex().z() - pv.z());
         GlobalVector jetDir(jet.px(),jet.py(),jet.pz());
         if(reco::deltaR(flightDir, jetDir ) <=  jetdRSV_) {//0.4
            isInJetDir = true; break;}
      }
      
      //distance in the transverse plane between the SV and PV
      Measurement1D d2d= vdistXY.distance(pv, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error()))); 
      
      //pointing angle (i.e. the angle between the sum of the momentum of the 
      //tracks in the SV and the flight direction betwen PV and SV)  
      double dx = (sv.vx() - pv.x()), dy = (sv.vy() - pv.y()), dz = (sv.vz() - pv.z());    
      double cosPdotV = (dx * sv.px() + dy*sv.py() + dz*sv.pz())/(sv.p()*std::sqrt(dx * dx + dy * dy + dz * dz));
      
      //Number of tracks
      int ndaus  = sv.numberOfDaughters();
      
      //significance of distance in 3d space point between the SV and PV
      Measurement1D dl= vdist.distance(pv, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
     

      //calculate sv b-tag discriminator
      std::map<std::string, double> sv_map;

      sv_map["Evt_pt"] = sv.pt();
      sv_map["Evt_eta"] = sv.eta();
      sv_map["Evt_mass"] = sv.mass();
      sv_map["Evt_d3d"] = dl.value();
      sv_map["Evt_d3dsig"] = dl.significance();
      sv_map["Evt_dxy"] = d2d.value();
      sv_map["Evt_costhetaSvPv"] = cosPdotV;
      sv_map["Evt_ndof"] = sv.vertexNdof();

      std::map<std::string, double> probs = m_SVDiscrTool.PROB(sv_map); 

      double probb, probc;
      
      probb = probs["prob_isB"];
      probc = probs["prob_isC"];


     //Add track based info
      double trackDeltaR(0.);   double trackPtRel(0.);   double trackPtRatio(0.);
      double trackEta(0.);      double trackEtaRel(0.);  double trackPPar(0.);     double trackPParRatio(0.);
      double trackSip2dVal(0.); double trackSip2dSig(0.);double trackSip3dVal(0.); double trackSip3dSig(0.);
      double trackSVDistVal(0.);double trackSVDistSig(0.);

      for (unsigned int i = 0; i < sv.numberOfDaughters(); i++) {
        auto const* cand = sv.daughter(i);
        auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
        auto pf_cand = dynamic_cast<const reco::PFCandidate*>(cand);

        // deal with PAT/AOD polymorphism to get the track
        const reco::Track *track_ptr = nullptr;

        if (pf_cand){
           track_ptr = pf_cand->bestTrack();  // trackRef was sometimes null
        }
        else if (packed_cand && packed_cand->hasTrackDetails()){ // if no track info is available, this will give an exception
           track_ptr = &(packed_cand->pseudoTrack());
        }

        math::XYZVector trackMom = track_ptr->momentum();
        double trackMag          = std::sqrt(trackMom.Mag2());
        TVector3 trackMom3(trackMom.x(), trackMom.y(), trackMom.z()); 

      //Fill the branches
        trackDeltaR    = reco::deltaR(trackMom, svDirU);
        trackPtRel     = trackMom3.Perp(svDir3); 
        trackEta       = trackMom.Eta();
        trackEtaRel    = reco::btau::etaRel(svDirU, trackMom);
        trackPPar      = svDirU.Dot(trackMom);
        trackPtRatio   = trackMom3.Perp(svDir3) / trackMag;
        trackPParRatio = svDirU.Dot(trackMom) / trackMag;

        // Use transient track that has a pointer to the magnetic field, and geometry to estimate  
        // track parameters along the track trajectory 
        reco::TransientTrack transientTrack;
        transientTrack = builder_->build(*track_ptr);

        //Extrapolate the track to the point of closest approach to the PV in transverse (3D) plane to compute IP and its error
        Measurement1D meas_ip2d = IPTools::signedTransverseImpactParameter(transientTrack, svDir, pv).second;
        Measurement1D meas_ip3d = IPTools::signedImpactParameter3D(transientTrack, svDir, pv).second;
        trackSip2dVal = static_cast<float>(meas_ip2d.value());
        trackSip2dSig = static_cast<float>(meas_ip2d.significance());
        trackSip3dVal = static_cast<float>(meas_ip3d.value());
        trackSip3dSig = static_cast<float>(meas_ip3d.significance());

        //Extrapolate the track to the point of closest approach to the SV axis
        Measurement1D svdist = IPTools::jetTrackDistance(transientTrack, svDir, pv).second;
        trackSVDistVal = static_cast<float>(svdist.value());
        trackSVDistSig = static_cast<float>(svdist.significance());

        // Fill the branches
        evt_.trkDeltaRSV.push_back(trackDeltaR); 
        evt_.trkPtRelSV.push_back(trackPtRel); 
        evt_.trkEtaSV.push_back(trackEta); 
        evt_.trkEtaRelSV.push_back(trackEtaRel); 
        evt_.trkPParSV.push_back(trackPPar);
        evt_.trkPtRatioSV.push_back(trackPtRatio);
        evt_.trkPParRatioSV.push_back(trackPParRatio); 
        evt_.trkSip2dValSV.push_back(trackSip2dVal);
        evt_.trkSip2dSigSV.push_back(trackSip2dSig);
        evt_.trkSip3dValSV.push_back(trackSip3dVal);
        evt_.trkSip3dSigSV.push_back(trackSip3dSig);
        evt_.trkSVDistValSV.push_back(trackSVDistVal);
        evt_.trkSVDistSigSV.push_back(trackSVDistSig); 
      }

      evt_.nGoodSV++;
      evt_.pt.push_back(sv.pt()); 
      evt_.eta.push_back(sv.eta());
      evt_.phi.push_back(sv.phi());
      evt_.energy.push_back(sv.energy());
      evt_.mass.push_back(sv.mass());
      evt_.chi2.push_back(sv.vertexChi2()); 
      evt_.normChi2.push_back(sv.vertexNormalizedChi2());
      evt_.dxy.push_back(d2d.value());
      evt_.dxyerr.push_back(d2d.error());
      evt_.dxysig.push_back(d2d.significance());  
      evt_.d3d.push_back(dl.value());
      evt_.d3derr.push_back(dl.error());
      evt_.d3dsig.push_back(dl.significance());
      evt_.costhetaSvPv.push_back(cosPdotV);
      evt_.ndau.push_back(ndaus);
      evt_.ndof.push_back(sv.vertexNdof());

      evt_.isMatchedToJet.push_back(isInJetDir);
      evt_.isB.push_back(isb);
      evt_.isGB.push_back(isgb);
      evt_.isC.push_back(isc);
      evt_.isG.push_back(isg);
      evt_.isS.push_back(iss);
      evt_.isUD.push_back(isud);
      evt_.isMatchOther.push_back(ismatchother);
      evt_.isOther.push_back(isother);

      evt_.nB.push_back(nbs);
      evt_.nC.push_back(ncs);
      evt_.nG.push_back(ngs);
      evt_.nS.push_back(nss);
      evt_.nUD.push_back(nuds);
      evt_.nMatchOther.push_back(nmatchothers);
      evt_.nOther.push_back(nothers);
      evt_.ProbB.push_back(probb);
      evt_.ProbC.push_back(probc);

      evt_.matchedGenIsHardProcess.push_back(matchGenIsHardProcess);
      evt_.matchedGenIsLastCopyBeforeFSR.push_back(matchGenIsLastCopyBeforeFSR);
      evt_.matchedGenPtSV.push_back(matchGenPt);
      evt_.matchedGenEtaSV.push_back(matchGenEta);
      evt_.matchedGenMassSV.push_back(matchGenMass);
      evt_.matchedGenPhiSV.push_back(matchGenPhi);
      evt_.matchedGenEnergySV.push_back(matchGenEnergy);
      evt_.matchedGenMomSV.push_back(matchGenMomID);
      evt_.matchedGenMindRSV.push_back(mindR);
      evt_.matchedGendRSV.push_back(matchGenSVdR);

   }
   //std::cout << "???" << std::endl; 
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
   tree_ = fs_->make<TTree>("svtree", "svtree") ;
   evt_.RegisterTree(tree_, "Evt") ;
   genevt_.RegisterTree(tree_, "GenEvt") ;
   m_SVDiscrTool.CreateNN(sv_model_file_);
   
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


const float& 
KUNtuplizer::catchInfs(const float& in,const float& replace_value) {
  if(in==in){
    if(std::isinf(in))
      return replace_value;
    else if(in < -1e32 || in > 1e32)
      return replace_value;
    return in;
  }
  return replace_value;
}

float 
KUNtuplizer::catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound) {
  float withoutinfs=catchInfs(in,replace_value);
  if(withoutinfs<lowerbound) return lowerbound;
  if(withoutinfs>upperbound) return upperbound;
  return withoutinfs;
}

//define this as a plug-in
DEFINE_FWK_MODULE(KUNtuplizer);
