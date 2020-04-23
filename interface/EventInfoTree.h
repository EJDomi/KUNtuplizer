#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class EventInfoTree {
  public:
    int    countEvents ;
    int    runno ;
    int    lumisec ;
    int    evtno;
    //int    npv;   
    //int    npuTrue;
    //int    puBX;
    //int    npuInt;
    int    nGoodVtx;
    //std::vector<int>   ndofVtx;
    //std::vector<float> chi2Vtx;
    //std::vector<float> zVtx;
    //std::vector<float> rhoVtx;
    //std::vector<float> ptVtx;

    //int    nSV;
    //int    nGoodSV;
    float pt;
    float eta;
    float phi;
    float energy;
    float mass;
    int   isMatchedToJet; 
    float chi2;
    float normChi2;
    float dxy;
    float dxyerr;
    float dxysig;
    float d3d;
    float d3derr;
    float d3dsig;
    float costhetaSvPv;
    int   ntrk_;
    float   ndau;
    float   ndof;
    int   isB;
    int   isGB;
    int   isC;
    int   isG;
    int   isS;
    int   isUD;
    int   isMatchOther;
    int   isOther;
    int   nB;
    int   nC;
    int   nG;
    int   nS;
    int   nUD;
    int   nMatchOther;
    int   nOther;
    float matchedGenPt;
    float matchedGenEta;
    float matchedGenMass;
    float matchedGenPhi;
    float ProbB;
    float ProbC;

    // track variables
    static constexpr size_t max_trk=10;
    float trkDeltaRSV[max_trk];
    float trkPtRelSV[max_trk];
    float trkEtaSV[max_trk];
    float trkEtaRelSV[max_trk];
    float trkPParSV[max_trk];
    float trkPtRatioSV[max_trk];  
    float trkPParRatioSV[max_trk];   
    float trkSip2dValSV[max_trk];
    float trkSip2dSigSV[max_trk];
    float trkSip3dValSV[max_trk];
    float trkSip3dSigSV[max_trk];
    float trkSVDistValSV[max_trk];
    float trkSVDistSigSV[max_trk];




   void clearTreeVectors() {
      //std::fill(std::begin(trkDeltaRSV), std::end(trkDeltaRSV), 0.);
      //std::fill(std::begin(trkPtRelSV), std::end(trkPtRelSV), 0.);
      //std::fill(std::begin(trkEtaSV), std::end(trkEtaSV), 0.);
      //std::fill(std::begin(trkEtaRelSV), std::end(trkEtaRelSV), 0.);
      //std::fill(std::begin(trkPParSV), std::end(trkPParSV), 0.);
      //std::fill(std::begin(trkPtRatioSV), std::end(trkPtRatioSV), 0.);
      //std::fill(std::begin(trkPParRatioSV), std::end(trkPParRatioSV), 0.);
      //std::fill(std::begin(trkSip2dValSV), std::end(trkSip2dValSV), 0.);
      //std::fill(std::begin(trkSip2dSigSV), std::end(trkSip2dSigSV), 0.);
      //std::fill(std::begin(trkSip3dValSV), std::end(trkSip3dValSV), 0.);
      //std::fill(std::begin(trkSip3dSigSV), std::end(trkSip3dSigSV), 0.);
      //std::fill(std::begin(trkSVDistValSV), std::end(trkSVDistValSV), 0.);
      //std::fill(std::begin(trkSVDistSigSV), std::end(trkSVDistSigSV), 0.);
    //   ndofVtx.clear();
    //   chi2Vtx.clear();
    //   zVtx.clear();
    //   rhoVtx.clear();
    //   ptVtx.clear();
    //   
    //   pt.clear();
    //   eta.clear();
    //   phi.clear();
    //   energy.clear();
    //   mass.clear();
    //   isMatchedToJet.clear();
    //   chi2.clear();
    //   normChi2.clear();
    //   dxy.clear();
    //   dxyerr.clear();
    //   dxysig.clear();
    //   d3d.clear();
    //   d3derr.clear();
    //   d3dsig.clear();       
    //   costhetaSvPv.clear();
    //   ntrk.clear();
    //   ndau.clear();
    //   ndof.clear();
    //   isB.clear();
    //   isC.clear();
    //   isOther.clear();
    //   matchedGenPt.clear();
    //   matchedGenEta.clear();
    //   matchedGenMass.clear();
    //   matchedGenPhi.clear();
    }

    void RegisterTree(TTree* tree, std::string name="SV") {
      //tree->Branch((name+"_countEvents").c_str(), &countEvents, (name+"_countEvents/I").c_str());
      tree->Branch((name+"_runno").c_str(),   &runno,    (name+"_runno/I").c_str());
      tree->Branch((name+"_lumisec").c_str(), &lumisec,  (name+"_lumisec/I").c_str());
      tree->Branch((name+"_evtno").c_str(),   &evtno,    (name+"_evtno/I").c_str());
      //tree->Branch((name+"_npv").c_str(),     &npv,      (name+"_npv/I").c_str());
      //tree->Branch((name+"_npuTrue").c_str(), &npuTrue,  (name+"_npuTrue/I").c_str());
      //tree->Branch((name+"_puBX").c_str(),    &puBX,     (name+"_puBX/I").c_str());
      //tree->Branch((name+"_npuInt").c_str(),  &npuInt,   (name+"_npuInt/I").c_str());      
      tree->Branch((name+"_nGoodVtx").c_str(),&nGoodVtx, (name+"_nGoodVtx/I").c_str());
      //tree->Branch((name+"_ndofVtx").c_str(), &ndofVtx);
      //tree->Branch((name+"_chi2Vtx").c_str(), &chi2Vtx);
      //tree->Branch((name+"_zVtx").c_str(),    &zVtx);
      //tree->Branch((name+"_rhoVtx").c_str(),  &rhoVtx); 
      //tree->Branch((name+"_ptVtx").c_str(),   &ptVtx);
      
      //tree->Branch((name+"_nSV").c_str(),             &nSV,      (name+"_nSV/I").c_str());
      //tree->Branch((name+"_nGoodSV").c_str(),         &nGoodSV,  (name+"_nGoodSV/I").c_str());
      tree->Branch((name+"_pt").c_str(),            &pt, (name+"_pt/F").c_str());
      tree->Branch((name+"_eta").c_str(),           &eta, (name+"_eta/F").c_str());
      tree->Branch((name+"_phi").c_str(),           &phi, (name+"_phi/F").c_str());
      tree->Branch((name+"_energy").c_str(),        &energy, (name+"_energy/F").c_str());
      tree->Branch((name+"_mass").c_str(),          &mass, (name+"_mass/F").c_str());
      tree->Branch((name+"_isMatchedToJet").c_str(),&isMatchedToJet, (name+"_isMatchedToJet/I").c_str());
      tree->Branch((name+"_chi2").c_str(),          &chi2, (name+"_chi2/F").c_str()); 
      tree->Branch((name+"_normChi2").c_str(),      &normChi2, (name+"_normChi2/F").c_str());
      tree->Branch((name+"_dxy").c_str(),           &dxy, (name+"_dxy/F").c_str()); 
      tree->Branch((name+"_dxyerr").c_str(),        &dxyerr, (name+"_dxyerr/F").c_str()); 
      tree->Branch((name+"_dxysig").c_str(),        &dxysig, (name+"_dxysig/F").c_str());
      tree->Branch((name+"_d3d").c_str(),           &d3d, (name+"_d3d/F").c_str()); 
      tree->Branch((name+"_d3derr").c_str(),        &d3derr, (name+"_d3derr/F").c_str()); 
      tree->Branch((name+"_d3dsig").c_str(),        &d3dsig, (name+"_d3dsig/F").c_str());
      tree->Branch((name+"_costhetaSvPv").c_str(),  &costhetaSvPv, (name+"_costhetaSvPv/F").c_str()); 
      tree->Branch((name+"_ntrk_").c_str(),          &ntrk_, (name+"_ntrk_/I").c_str());
      tree->Branch((name+"_ndau").c_str(),          &ndau, (name+"_ndau/F").c_str());
      tree->Branch((name+"_ndof").c_str(),          &ndof, (name+"_ndof/F").c_str());
      tree->Branch((name+"_isB").c_str(),           &isB, (name+"_isB/i").c_str());
      tree->Branch((name+"_isGB").c_str(),           &isGB, (name+"_isGB/i").c_str());
      tree->Branch((name+"_isC").c_str(),           &isC, (name+"_isC/i").c_str());
      tree->Branch((name+"_isG").c_str(),           &isG, (name+"_isG/i").c_str());
      tree->Branch((name+"_isS").c_str(),           &isS, (name+"_isS/i").c_str());
      tree->Branch((name+"_isUD").c_str(),           &isUD, (name+"_isUD/i").c_str());
      tree->Branch((name+"_isMatchOther").c_str(),         &isMatchOther, (name+"_isMatchOther/i").c_str());
      tree->Branch((name+"_isOther").c_str(),         &isOther, (name+"_isOther/i").c_str());
      tree->Branch((name+"_nB").c_str(),           &nB, (name+"_nB/I").c_str());
      tree->Branch((name+"_nC").c_str(),           &nC, (name+"_nC/I").c_str());
      tree->Branch((name+"_nG").c_str(),           &nG, (name+"_nG/I").c_str());
      tree->Branch((name+"_nS").c_str(),           &nS, (name+"_nS/I").c_str());
      tree->Branch((name+"_nUD").c_str(),           &nUD, (name+"_nUD/I").c_str());
      tree->Branch((name+"_nMatchOther").c_str(),         &nMatchOther, (name+"_nMatchOther/I").c_str());
      tree->Branch((name+"_nOther").c_str(),         &nOther, (name+"_nOther/I").c_str());
      tree->Branch((name+"_ProbB").c_str(),         &ProbB, (name+"_ProbB/F").c_str());
      tree->Branch((name+"_ProbC").c_str(),         &ProbB, (name+"_ProbC/F").c_str());
      //tree->Branch((name+"_matchedGenPt").c_str(),  &matchedGenPt, (name+"_matchedGenPt").c_str());
      //tree->Branch((name+"_matchedGenEta").c_str(), &matchedGenEta, (name+"_matchedGenEta").c_str());
      //tree->Branch((name+"_matchedGenMass").c_str(),&matchedGenMass, (name+"_matchedGenMass").c_str());
      //tree->Branch((name+"_matchedGenPhi").c_str(), &matchedGenPhi, (name+"_matchedGenPhi").c_str());
      //
      // track variables
      tree->Branch((name+"_trkDeltaRSV").c_str(), &trkDeltaRSV, (name+"_trkDeltaRSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkPtRelSV").c_str(), &trkPtRelSV, (name+"_trkPtRelSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkEtaSV").c_str(), &trkEtaSV, (name+"_trkEtaSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkEtaRelSV").c_str(), &trkEtaRelSV, (name+"_trkEtaRelSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkPParSV").c_str(), &trkPParSV, (name+"_trkPParSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkPtRatioSV").c_str(), &trkPtRatioSV, (name+"_trkPtRatioSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkPParRatioSV").c_str(), &trkPParRatioSV, (name+"_trkPParRatioSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSip2dValSV").c_str(), &trkSip2dValSV, (name+"_trkSip2dValSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSip2dSigSV").c_str(), &trkSip2dSigSV, (name+"_trkSip2dSigSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSip3dValSV").c_str(), &trkSip3dValSV, (name+"_trkSip3dValSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSip3dSigSV").c_str(), &trkSip3dSigSV, (name+"_trkSip3dSigSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSVDistValSV").c_str(), &trkSVDistValSV, (name+"_trkSVDistValSV["+name+"_ntrk_]/F").c_str());
      tree->Branch((name+"_trkSVDistSigSV").c_str(), &trkSVDistSigSV, (name+"_trkSVDistSigSV["+name+"_ntrk_]/F").c_str());

   
    }
};
