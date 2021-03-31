// -*- C++ -*-
//
// Package:    GenAnalysis/TQGenAnalyzer
// Class:      TQGenAnalyzer
//
/**\class TQGenAnalyzer TQGenAnalyzer.cc GenAnalysis/TQGenAnalyzer/plugins/TQGenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Livia Soffi
//         Created:  Mon, 01 Mar 2021 11:15:18 GMT
//
//


// system include files
#include <memory>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"  

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include <SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h>

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronIDAlgo.h"


#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"



#include "DataFormats/Math/interface/deltaR.h"
// root include files
#include "TTree.h"
#include "TLorentzVector.h"


namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }


struct tree_struc_{
  int run;
  int lumi;
  long unsigned int event;
  int nEle;
  int nMu;
  float TQ_genMass;

  std::vector<float>            genLep_pt;
  std::vector<float>            genLep_eta;
  std::vector<float>            genLep_mass;
  std::vector<float>            genLep_phi;
  std::vector<int>            genLep_pdgId;
  std::vector<float>            genMom_pt;
  std::vector<float>            genMom_eta;
  std::vector<float>            genMom_mass;
  std::vector<float>            genMom_phi;
  std::vector<int>            genMom_pdgId;

  int nEleReco;
  int nMuReco;
  float TQ_recoMass;
  std::vector<float>            recoMu_pt;
  std::vector<float>            recoMu_eta;
  std::vector<float>            recoMu_mass;
  std::vector<float>            recoMu_phi;
  std::vector<int>              recoMu_charge;
  std::vector<float>            recoMu_dR1;
  std::vector<float>            recoMu_dR2;
  std::vector<int>              recoMu_matchid1;
  std::vector<int>              recoMu_matchid2;
  std::vector<int>              recoMu_looseID;
  std::vector<int>              recoMu_softID;

  //below variables are those applied in softID. Cuts can be seen here:https://cmssdt.cern.ch/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc?v=CMSSW_10_6_20
  std::vector<int>              recoMu_nTrkLayers;
  std::vector<int>              recoMu_nPixLayers;
  std::vector<int>              recoMu_isHQ;
  std::vector<float>              recoMu_dxy;
  std::vector<float>              recoMu_dz;

  std::vector<float>            recoEle_pt;
  std::vector<float>            recoEle_eta;
  std::vector<float>            recoEle_mass;
  std::vector<float>            recoEle_phi;
  std::vector<int>              recoEle_charge;

  std::vector<int>              recoEle_vx;
  std::vector<int>              recoEle_vy;
  std::vector<int>              recoEle_vz;

  std::vector<float>            recoEle_dR1;
  std::vector<float>            recoEle_dR2;
  std::vector<int>              recoEle_matchid1;
  std::vector<int>              recoEle_matchid2;
  std::vector<float>recoEle_E_ecal_preReg;
  std::vector<float>recoEle_E_ecal_postReg;
  std::vector<float>recoEle_E_ecaltrk_preReg;
  std::vector<float>recoEle_E_ecaltrk_postReg;

  std::vector<float>recoEle_rawE;
  std::vector<float>recoEle_corrEcalE;
  std::vector<float>recoEle_gsfTrkChi2;
  std::vector<float>recoEle_passConvVeto;
  std::vector<float>recoEle_isPF;
  std::vector<float>recoEle_mvaPFID;
  std::vector<float>recoEle_mvaPFValue;
  std::vector<float>recoEle_isPFoverlap;

  //specific fro low pt electrons
  std::vector<float>recoEle_isLowPt;
  std::vector<float>recoEle_mvaID;
  std::vector<float>recoEle_mvaValue;

  std::vector<float>recoEle_mvaEGammaID;
  std::vector<float>recoEle_mvaEGammaValue;
  
  float vx;
  float vy;
  float vz;

};


class TQGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TQGenAnalyzer(const edm::ParameterSet&);
      ~TQGenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void initTreeStructure();
  void clearVectors(); 

      // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GenEventInfoProduct>genInfoToken_;
  
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;


  const edm::EDGetTokenT< edm::View<pat::Muon> > patMuons_; // MINIAOD
  edm::Handle< edm::View<pat::Muon> > patMuonsH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectrons_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsH_;

  const edm::EDGetTokenT< edm::View<reco::GsfElectron> > patElectronsLowPt_; // MINIAOD
  edm::Handle< edm::View<reco::GsfElectron> > gsfElectronsLowPtH_;

  const edm::EDGetTokenT< std::vector<reco::Vertex> > vtx_; // MINIAOD
  edm::Handle< std::vector<reco::Vertex>> vtxH_;


  const double dr_cleaning_;
  const double dz_cleaning_;

  //setup mva ID for electrons
  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValueEGamma_; 
  edm::Handle< edm::ValueMap<float> > mvaValueEGammaH_;
  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaIdEGamma_; 
  edm::Handle< edm::ValueMap<bool> > mvaIdEGammaH_;


  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValuePF_; 
  edm::Handle< edm::ValueMap<float> > mvaValuePFH_;
  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaIdPF_; 
  edm::Handle< edm::ValueMap<bool> > mvaIdPFH_;


  const edm::EDGetTokenT< edm::ValueMap<float> > mvaValue_; 
  edm::Handle< edm::ValueMap<float> > mvaValueH_;
  const edm::EDGetTokenT< edm::ValueMap<bool> > mvaId_; 
  edm::Handle< edm::ValueMap<bool> > mvaIdH_;

  // setup tree;
  TTree* tree;
  tree_struc_ tree_;

  //setup regression for low pt
  std::unique_ptr<ModifyObjectValueBase> regression_;     // Low pt 
  
};


TQGenAnalyzer::TQGenAnalyzer(const edm::ParameterSet& iConfig)
  :
  prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  genParticlesH_(),
  patMuons_(consumes< edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("patMuons"))),
  patMuonsH_(),
  patElectrons_(consumes< edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("patElectrons"))),
  gsfElectronsH_(),
  patElectronsLowPt_(consumes< edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("patElectronsLowPt"))),
  gsfElectronsLowPtH_(),
  vtx_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vtx"))),
  vtxH_(),
  dr_cleaning_(iConfig.getParameter<double>("drForCleaning")),
  dz_cleaning_(iConfig.getParameter<double>("dzForCleaning")),
  mvaValueEGamma_(consumes< edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValueEGamma"))),
  mvaValueEGammaH_(),
  mvaIdEGamma_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaIdEGamma"))),
  mvaIdEGammaH_(),
  mvaValuePF_(consumes< edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuePF"))),
  mvaValuePFH_(),
  mvaIdPF_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaIdPF"))),
  mvaIdPFH_(),
  mvaValue_(consumes< edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValue"))),
  mvaValueH_(),
  mvaId_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("mvaId"))),
  mvaIdH_()



{
  
  // Regression stuff - lowPtElectrons
  if( iConfig.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
    const edm::ParameterSet& iconf = iConfig.getParameterSet("gsfRegressionConfig");
    const std::string& mname = iconf.getParameter<std::string>("modifierName");
    //    std::cout<<mname<<std::endl;
    ModifyObjectValueBase* plugin = ModifyObjectValueFactory::get()->create(mname,iconf);
    regression_.reset(plugin);
    edm::ConsumesCollector sumes = consumesCollector();
    regression_->setConsumes(sumes);
    
  } else {
    regression_.reset(nullptr);
  }

}


TQGenAnalyzer::~TQGenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TQGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    iEvent.getByToken(prunedGenParticles_, genParticlesH_);
    iEvent.getByToken(patMuons_,patMuonsH_);
    iEvent.getByToken(patElectrons_,gsfElectronsH_);
    if(gsfElectronsLowPtH_.isValid()) iEvent.getByToken(patElectronsLowPt_,gsfElectronsLowPtH_);
    iEvent.getByToken(vtx_,vtxH_);
    iEvent.getByToken(mvaIdEGamma_, mvaIdEGammaH_); 
    iEvent.getByToken(mvaValueEGamma_, mvaValueEGammaH_); 
    iEvent.getByToken(mvaIdPF_, mvaIdPFH_); 
    iEvent.getByToken(mvaValuePF_, mvaValuePFH_); 
    iEvent.getByToken(mvaId_, mvaIdH_); 
    iEvent.getByToken(mvaValue_, mvaValueH_); 

    using namespace edm;
    regression_->setEvent(iEvent);
    regression_->setEventContent(iSetup);
  

    std::vector<float>genLep_pt;
    std::vector<float>genLep_mass;
    std::vector<float>genLep_eta;
    std::vector<float>genLep_phi;
    std::vector<int>genLep_pdgId;

    std::vector<float>genMom_pt;
    std::vector<float>genMom_mass;
    std::vector<float>genMom_eta;
    std::vector<float>genMom_phi;
    std::vector<int>genMom_pdgId;

    std::vector<float>recoMu_pt;
    std::vector<float>recoMu_mass;
    std::vector<float>recoMu_eta;
    std::vector<float>recoMu_phi;
    std::vector<float>recoMu_dR1;
    std::vector<float>recoMu_dR2;
    std::vector<int>              recoMu_charge;
    std::vector<int>              recoMu_matchid1;
    std::vector<int>              recoMu_matchid2;
    std::vector<int>              recoMu_looseID;
    std::vector<int>              recoMu_softID;
    std::vector<int>              recoMu_nTrkLayers;
    std::vector<int>              recoMu_nPixLayers;
    std::vector<int>              recoMu_isHQ;
    std::vector<float>              recoMu_dxy;
    std::vector<float>              recoMu_dz;



    std::vector<float>recoEle_pt;
    std::vector<float>recoEle_mass;
    std::vector<float>recoEle_eta;
    std::vector<float>recoEle_phi;
    std::vector<float>recoEle_dR1;
    std::vector<float>recoEle_dR2;
    std::vector<int>              recoEle_charge;

    std::vector<int>              recoEle_vx;
    std::vector<int>              recoEle_vy;
    std::vector<int>              recoEle_vz;


    std::vector<int>              recoEle_matchid1;
    std::vector<int>              recoEle_matchid2;
    std::vector<float>recoEle_E_ecal_preReg;
    std::vector<float>recoEle_E_ecal_postReg;
    std::vector<float>recoEle_E_ecaltrk_preReg;
    std::vector<float>recoEle_E_ecaltrk_postReg;
    std::vector<float>recoEle_rawE;
    std::vector<float>recoEle_corrEcalE;
    std::vector<float>recoEle_gsfTrkChi2;
    std::vector<float>recoEle_passConvVeto;
    std::vector<float>recoEle_isPF;
    std::vector<float>recoEle_mvaPFID;
    std::vector<float>recoEle_mvaPFValue;

    std::vector<float>recoEle_isPFoverlap;

    std::vector<float>recoEle_isLowPt;
    std::vector<float>recoEle_mvaID;
    std::vector<float>recoEle_mvaValue;

    std::vector<float>recoEle_mvaEGammaID;
    std::vector<float>recoEle_mvaEGammaValue;

    // --- general event info 
    unsigned long int event = iEvent.id().event();   
    int run                 = iEvent.id().run();
    int lumi                = iEvent.luminosityBlock();
    int nEle=0;
    int nMu=0;
    int nMuReco=0;
    int nEleReco=0;

    TLorentzVector* lep1=new TLorentzVector();
    TLorentzVector* lep2=new TLorentzVector();
    TLorentzVector* lep3=new TLorentzVector();
    TLorentzVector* lep0=new TLorentzVector();
    
    float vx;
    float vy;
    float vz;
    //    std::cout<<"-----------------------------------------------------"<<std::endl;

    // save gen particles (only leptons)
    
    for (const auto & genpar_iter : *genParticlesH_){ // loop over genparticles
      if(abs(genpar_iter.pdgId())!=13 && abs(genpar_iter.pdgId())!=11)continue;

      if(genpar_iter.mother(0)->pdgId()<20 || genpar_iter.mother(0)->pdgId()>50) continue;
      
      //      std::cout<<"status: "<<genpar_iter.status()<<" pdgid: "<<genpar_iter.pdgId()<<" pt: "<<genpar_iter.pt()<<" mass: "<<genpar_iter.mass()<<" eta: "<<genpar_iter.eta()<<" phi: "<<genpar_iter.phi()<<" mom pdgId: "<<genpar_iter.mother(0)->pdgId()<<" mom mass: "<<genpar_iter.mother(0)->mass()<<std::endl;

      
      genLep_pt.push_back(genpar_iter.pt());
      genLep_mass.push_back(genpar_iter.mass());
      genLep_eta.push_back(genpar_iter.eta());
      genLep_phi.push_back(genpar_iter.phi());
      genLep_pdgId.push_back(genpar_iter.pdgId());
      

      genMom_pt.push_back(genpar_iter.mother(0)->pt());
      genMom_mass.push_back(genpar_iter.mother(0)->mass());
      genMom_eta.push_back(genpar_iter.mother(0)->eta());
      genMom_phi.push_back(genpar_iter.mother(0)->phi());
      genMom_pdgId.push_back(genpar_iter.mother(0)->pdgId());
      
      if(abs(genpar_iter.pdgId())==11)nEle++;
      if(abs(genpar_iter.pdgId())==13)nMu++;
      
    }


    //saving info of primary vertex
    const reco::Vertex & pv = vtxH_->front(); 
    vx=pv.x();
    vy=pv.y();
    vz=pv.z();


    //looping over muons
    for (const auto & mu_iter : *patMuonsH_){
      nMuReco++;

      float dR=999.;
      float dRMin1=999.;
      float dRMin2=999.;
      int recomatchid1=999;
      int recomatchid2=999;
      int recolooseID=999;
      int recosoftID=999;     
      int recontrklayers=999;
      int reconpixlayers=999;
      int recoishq=999;
      float recodxy=999.;
      float recodz=999.;


      float recopt;
      float recoeta;
      float recophi;
      float recomass;
      int recocharge=999;


      float eta=mu_iter.eta();
      float phi=mu_iter.phi();
      
      //store kinematics
      recopt=mu_iter.pt();
      recoeta=mu_iter.eta();
      recophi=mu_iter.phi();
      recomass=mu_iter.mass();
      recocharge=mu_iter.charge();
      

      //store ID variables
      recolooseID=muon::isLooseMuon(mu_iter);
      recosoftID=muon::isSoftMuon(mu_iter, pv);
      //      recontrklayers=mu_iter.innerTrack()->hitPattern().trackerLayersWithMeasurement();
      //reconpixlayers=mu_iter.innerTrack()->hitPattern().pixelLayersWithMeasurement();
      //recoishq=1;//mu_iter.innerTrack()->quality(reco::Track::highPurity); 
      //      recodxy=mu_iter.innerTrack()->dxy(pv.position());
      //recodz=mu_iter.innerTrack()->dz(pv.position());




      //compute dR matching w/ each gen muon      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=13)continue;
	dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin1){
	  dRMin1=dR;
	  recomatchid1=i;
	}
      }
      
      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=13)continue;
	if(i==(unsigned int)recomatchid1)continue;
	dR= deltaR(eta,phi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin2){
	  dRMin2=dR;
	  recomatchid2=i;
	}
      }
      

      
      recoMu_pt.push_back(recopt);
      recoMu_mass.push_back(recomass);
      recoMu_eta.push_back(recoeta);
      recoMu_phi.push_back(recophi);
      recoMu_charge.push_back(recocharge);

      recoMu_dR1.push_back(dRMin1);
      recoMu_dR2.push_back(dRMin2);
      recoMu_matchid1.push_back(recomatchid1);
      recoMu_matchid2.push_back(recomatchid2);
      recoMu_looseID.push_back(recolooseID);
      recoMu_softID.push_back(recosoftID);
      recoMu_nTrkLayers.push_back(recontrklayers);
      recoMu_nPixLayers.push_back(reconpixlayers);
      recoMu_isHQ.push_back(recoishq);
      recoMu_dxy.push_back(recodxy);
      recoMu_dz.push_back(recodz);
    }




    //here we loop on PF electrons only
    for( size_t electronlooper = 0; electronlooper < gsfElectronsH_->size(); electronlooper++ ) {

      // filter candidates: there must be a trk with pT>0.5
      // Low pT electrons
      const reco::GsfElectronPtr ele(gsfElectronsH_, electronlooper);
      
      float dR=999.;
      float dRMin1=999.;
      float dRMin2=999.;
      int recomatchid1=999;
      int recomatchid2=999;

      float recopt=ele->pt();
      float recoeta=ele->eta();
      float recophi=ele->phi();
      float recomass=ele->mass();
      int recocharge=ele->charge();

      float recovx=ele->vx();
      float recovy=ele->vy();
      float recovz=ele->vz();

      
      float recorawe=ele->superCluster()->rawEnergy();
      float recocorrecale=ele->correctedEcalEnergy();
      float recotrkchi2=ele->gsfTrack()->normalizedChi2();
      int recoconvveto=999.;//ele->passConversionVeto();
      

      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=11)continue;
	dR= deltaR(recoeta,recophi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin1){
	  dRMin1=dR;
	  recomatchid1=i;
	}
      }
      
      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=11)continue;
	if(i==(unsigned int)recomatchid1)continue;
	dR= deltaR(recoeta,recophi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin2){
	  dRMin2=dR;
	  recomatchid2=i;
	}
      }
      
      int recoispf=1;
      int recoislowpt=0;
      int recoispfoverlap=0;


      //ID stuff
      

      //save EGM mva ID no Iso
      float mvaEGamma_value = -999.;
      int mvaEGamma_id = -999;

      //      std::cout<<" mva: "<<mvaValueEGammaH_->size()<<"  gsf: "<<gsfElectronsH_->size()<<std::endl;
      if ( mvaValueEGammaH_.isValid() && 
	   mvaValueEGammaH_->size() == gsfElectronsH_->size() ) {
	mvaEGamma_value = mvaValueEGammaH_->get( ele.key() );
	mvaEGamma_id = mvaIdEGammaH_->get( ele.key() );
      } else {
	std::cout << "ERROR! Issue matching EGamma MVA output to GsfElectrons!" << std::endl;
      }
 
      //save retrained EGamma ID for low pt PF
     float mvaPF_value = -999.;
      int mvaPF_id = -999;

      //      std::cout<<" mva: "<<mvaValuePFH_->size()<<"  gsf: "<<gsfElectronsH_->size()<<std::endl;
      if ( mvaValuePFH_.isValid() && 
	   mvaValuePFH_->size() == gsfElectronsH_->size() ) {
	mvaPF_value = mvaValuePFH_->get( ele.key() );
	mvaPF_id = mvaIdPFH_->get( ele.key() );
      } else {
	std::cout << "ERROR! Issue matching PF MVA output to GsfElectrons!" << std::endl;
      }


      float mva_value = -999.;
      int mva_id = -999;
      
      
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_charge.push_back(recocharge);

      recoEle_vx.push_back(recovx);
      recoEle_vy.push_back(recovy);
      recoEle_vz.push_back(recovz);

      recoEle_dR1.push_back(dRMin1);
      recoEle_dR2.push_back(dRMin2);
      recoEle_matchid1.push_back(recomatchid1);
      recoEle_matchid2.push_back(recomatchid2);

      recoEle_rawE.push_back(recorawe);
      recoEle_corrEcalE.push_back(recocorrecale);
      recoEle_gsfTrkChi2.push_back(recotrkchi2);
      recoEle_passConvVeto.push_back(recoconvveto);

      recoEle_isPF.push_back(recoispf);
      recoEle_isLowPt.push_back(recoislowpt);
      recoEle_isPFoverlap.push_back(recoispfoverlap);

      recoEle_mvaEGammaID.push_back(mvaEGamma_id);
      recoEle_mvaEGammaValue.push_back(mvaEGamma_value);

      recoEle_mvaPFID.push_back(mvaPF_id);
      recoEle_mvaPFValue.push_back(mvaPF_value);

      recoEle_mvaID.push_back(mva_id);
      recoEle_mvaValue.push_back(mva_value);

    // Here I apply the regression and i save also the postregression quantitiess
    TVector3 eleTV3(0,0,0);
    eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());
    // Regression stuff
    float recopre_ecal     = -999.;
    float recopre_ecaltrk  = -999.;
    float recopost_ecal    = -999.;
    float recopost_ecaltrk = -999.;

    reco::GsfElectron newElectron(*ele);
    recopre_ecal = newElectron.correctedEcalEnergy();
    recopre_ecaltrk = newElectron.energy();
    
    //std::cout<<" pre ecal: "<<recopre_ecal<<" pre ecaltrk: "<<recopre_ecaltrk<<std::endl;

    //    regression_->modifyObject(newElectron);
    recopost_ecal = newElectron.correctedEcalEnergy();
    recopost_ecaltrk = newElectron.energy();

    //std::cout<<" post ecal: "<<recopost_ecal<<" post ecaltrk: "<<recopost_ecaltrk<<std::endl;

    recoEle_E_ecal_preReg.push_back(recopre_ecal);
    recoEle_E_ecal_postReg.push_back(recopost_ecal);
    recoEle_E_ecaltrk_preReg.push_back(recopre_ecaltrk);
    recoEle_E_ecaltrk_postReg.push_back(recopost_ecaltrk);

}

    
    unsigned int nPFEle=recoEle_pt.size();

    if(gsfElectronsLowPtH_.isValid()){

    //here we loop on LowPt electrons only
    for( size_t electronlooper = 0; electronlooper < gsfElectronsLowPtH_->size(); electronlooper++ ) {

      // filter candidates: there must be a trk with pT>0.5
      // Low pT electrons
      const reco::GsfElectronPtr ele(gsfElectronsLowPtH_, electronlooper);
      
      float dR=999.;
      float dRMin1=999.;
      float dRMin2=999.;
      int recomatchid1=999;
      int recomatchid2=999;

      float recopt=ele->pt();
      float recoeta=ele->eta();
      float recophi=ele->phi();
      float recomass=ele->mass();
      int recocharge=ele->charge();

      float recovx=ele->vx();
      float recovy=ele->vy();
      float recovz=ele->vz();

      int recoispfoverlap=0;

      //remove overlap with PF
      //pf cleaning    
      bool clean_out = false;
      for(unsigned int iEle=0; iEle<nPFEle; ++iEle) {

	clean_out |= (
		      fabs(recoEle_vz[iEle] - recovz) < dz_cleaning_ &&
		      reco::deltaR(recoeta, recophi, recoEle_eta[iEle], recoEle_phi[iEle]) < dr_cleaning_   );

      }

            
      int recoispf=0;
      int recoislowpt=1;
      if(clean_out) recoispfoverlap=1;
      else recoispfoverlap=0;
 

      
      float recorawe=ele->superCluster()->rawEnergy();
      float recocorrecale=ele->correctedEcalEnergy();
      float recotrkchi2=ele->gsfTrack()->normalizedChi2();
      int recoconvveto=999.;//ele->passConversionVeto();
      

      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=11)continue;
	dR= deltaR(recoeta,recophi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin1){
	  dRMin1=dR;
	  recomatchid1=i;
	}
      }
      
      
      for(unsigned int i=0;i<genLep_pdgId.size();i++){
	if(abs(genLep_pdgId[i])!=11)continue;
	if(i==(unsigned int)recomatchid1)continue;
	dR= deltaR(recoeta,recophi,genLep_eta[i],genLep_phi[i]);
	if(dR<dRMin2){
	  dRMin2=dR;
	  recomatchid2=i;
	}
      }



      //ID stuff
      


      //save EGM mva ID no Iso
      float mvaEGamma_value = -999.;
      int mvaEGamma_id = -999;
      //save retrained EGamma ID for low pt PF
      float mvaPF_value = -999.;
      int mvaPF_id = -999;



     float mva_value = -999.;
      int mva_id = -999;

      //      std::cout<<" mva: "<<mvaValueH_->size()<<"  gsf: "<<gsfElectronsH_->size()<<std::endl;
      if ( mvaValueH_.isValid() && 
	   mvaValueH_->size() == gsfElectronsLowPtH_->size() ) {
	mva_value = mvaValueH_->get( ele.key() );
	//	mva_id = mvaIdH_->get( ele.key() );
      } else {
	std::cout << "ERROR! Issue matching Low Pt MVA output to GsfElectrons!" << std::endl;
      }
      
      
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_charge.push_back(recocharge);

      recoEle_vx.push_back(recovx);
      recoEle_vy.push_back(recovy);
      recoEle_vz.push_back(recovz);


      recoEle_dR1.push_back(dRMin1);
      recoEle_dR2.push_back(dRMin2);
      recoEle_matchid1.push_back(recomatchid1);
      recoEle_matchid2.push_back(recomatchid2);

      recoEle_rawE.push_back(recorawe);
      recoEle_corrEcalE.push_back(recocorrecale);
      recoEle_gsfTrkChi2.push_back(recotrkchi2);
      recoEle_passConvVeto.push_back(recoconvveto);

      recoEle_isPF.push_back(recoispf);
      recoEle_isLowPt.push_back(recoislowpt);
      recoEle_isPFoverlap.push_back(recoispfoverlap);

      recoEle_mvaEGammaID.push_back(mvaEGamma_id);
      recoEle_mvaEGammaValue.push_back(mvaEGamma_value);

      recoEle_mvaPFID.push_back(mvaPF_id);
      recoEle_mvaPFValue.push_back(mvaPF_value);

      recoEle_mvaID.push_back(mva_id);
      recoEle_mvaValue.push_back(mva_value);

    // Here I apply the regression and i save also the postregression quantitiess
    TVector3 eleTV3(0,0,0);
    eleTV3.SetPtEtaPhi(ele->pt(), ele->eta(), ele->phi());
    // Regression stuff
    float recopre_ecal     = -999.;
    float recopre_ecaltrk  = -999.;
    float recopost_ecal    = -999.;
    float recopost_ecaltrk = -999.;

    reco::GsfElectron newElectron(*ele);
    recopre_ecal = newElectron.correctedEcalEnergy();
    recopre_ecaltrk = newElectron.energy();
    
    //std::cout<<" pre ecal: "<<recopre_ecal<<" pre ecaltrk: "<<recopre_ecaltrk<<std::endl;

    //    regression_->modifyObject(newElectron);
    recopost_ecal = newElectron.correctedEcalEnergy();
    recopost_ecaltrk = newElectron.energy();

    //std::cout<<" post ecal: "<<recopost_ecal<<" post ecaltrk: "<<recopost_ecaltrk<<std::endl;

    recoEle_E_ecal_preReg.push_back(recopre_ecal);
    recoEle_E_ecal_postReg.push_back(recopost_ecal);
    recoEle_E_ecaltrk_preReg.push_back(recopre_ecaltrk);
    recoEle_E_ecaltrk_postReg.push_back(recopost_ecaltrk);

}
 

    }//close if low pt collection is valid

    // --- setup tree values
    initTreeStructure();

    // --- clear all vectors
    clearVectors();
  
    // --- fill trees
    tree_.run=run;
    tree_.event=event;
    tree_.lumi=lumi;
    tree_.nEle=nEle;
    tree_.nMu=nMu;
    tree_.vx=vx;
    tree_.vy=vy;
    tree_.vz=vz;

    for(unsigned int i=0;i<genLep_pt.size();i++){
      tree_.genLep_pt.push_back(genLep_pt[i]);
      tree_.genLep_mass.push_back(genLep_mass[i]);
      tree_.genLep_eta.push_back(genLep_eta[i]);
      tree_.genLep_phi.push_back(genLep_phi[i]);
      tree_.genLep_pdgId.push_back(genLep_pdgId[i]);

      tree_.genMom_pt.push_back(genMom_pt[i]);
      tree_.genMom_mass.push_back(genMom_mass[i]);
      tree_.genMom_eta.push_back(genMom_eta[i]);
      tree_.genMom_phi.push_back(genMom_phi[i]);
      tree_.genMom_pdgId.push_back(genMom_pdgId[i]);
    

      if(i==0)lep0->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==1)lep1->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==2)lep2->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
      if(i==3)lep3->SetPtEtaPhiM(genLep_pt[i],genLep_eta[i],genLep_phi[i],genLep_mass[i]);
    }

    TLorentzVector TQ(*lep0+*lep1+*lep2+*lep3);
    tree_.TQ_genMass=TQ.M();

    tree_.nMuReco=nMuReco;
    for(unsigned int i=0;i<recoMu_pt.size();i++){
      tree_.recoMu_pt.push_back(recoMu_pt[i]);
      tree_.recoMu_mass.push_back(recoMu_mass[i]);
      tree_.recoMu_eta.push_back(recoMu_eta[i]);
      tree_.recoMu_phi.push_back(recoMu_phi[i]);   
      tree_.recoMu_charge.push_back(recoMu_charge[i]);   

      tree_.recoMu_dR1.push_back(recoMu_dR1[i]);   
      tree_.recoMu_matchid1.push_back(recoMu_matchid1[i]);   
      tree_.recoMu_dR2.push_back(recoMu_dR2[i]);   
      tree_.recoMu_matchid2.push_back(recoMu_matchid2[i]);   
      tree_.recoMu_looseID.push_back(recoMu_looseID[i]);   
      tree_.recoMu_softID.push_back(recoMu_softID[i]);   
      tree_.recoMu_nTrkLayers.push_back(recoMu_nTrkLayers[i]);   
      tree_.recoMu_nPixLayers.push_back(recoMu_nPixLayers[i]);   
      tree_.recoMu_isHQ.push_back(recoMu_isHQ[i]);   
      tree_.recoMu_dxy.push_back(recoMu_dxy[i]);   
      tree_.recoMu_dz.push_back(recoMu_dz[i]);   
      
      //      if(i==0)recolep0->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);
      // if(i==1)recolep1->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);

    }


    nEleReco=recoEle_pt.size();
    tree_.nEleReco=nEleReco;
    for(unsigned int i=0;i<recoEle_pt.size();i++){
      tree_.recoEle_pt.push_back(recoEle_pt[i]);
      tree_.recoEle_mass.push_back(recoEle_mass[i]);
      tree_.recoEle_eta.push_back(recoEle_eta[i]);
      tree_.recoEle_phi.push_back(recoEle_phi[i]);   
      tree_.recoEle_charge.push_back(recoEle_charge[i]);   

      tree_.recoEle_vx.push_back(recoEle_vx[i]);   
      tree_.recoEle_vy.push_back(recoEle_vy[i]);   
      tree_.recoEle_vz.push_back(recoEle_vz[i]);   

      tree_.recoEle_dR1.push_back(recoEle_dR1[i]);   
      tree_.recoEle_matchid1.push_back(recoEle_matchid1[i]);   
      tree_.recoEle_dR2.push_back(recoEle_dR2[i]);   
      tree_.recoEle_matchid2.push_back(recoEle_matchid2[i]);   
      tree_.recoEle_E_ecal_preReg.push_back(recoEle_E_ecal_preReg[i]);   
      tree_.recoEle_E_ecal_postReg.push_back(recoEle_E_ecal_postReg[i]);   
      tree_.recoEle_E_ecaltrk_preReg.push_back(recoEle_E_ecaltrk_preReg[i]);   
      tree_.recoEle_E_ecaltrk_postReg.push_back(recoEle_E_ecaltrk_postReg[i]);   

      tree_.recoEle_rawE.push_back(recoEle_rawE[i]);
      tree_.recoEle_corrEcalE.push_back(recoEle_corrEcalE[i]);
      tree_.recoEle_gsfTrkChi2.push_back(recoEle_gsfTrkChi2[i]);
      tree_.recoEle_passConvVeto.push_back(recoEle_passConvVeto[i]);
      tree_.recoEle_isPF.push_back(recoEle_isPF[i]);

      tree_.recoEle_isLowPt.push_back(recoEle_isLowPt[i]);
      tree_.recoEle_mvaID.push_back(recoEle_mvaID[i]);
      tree_.recoEle_isPFoverlap.push_back(recoEle_isPFoverlap[i]);
      

      tree_.recoEle_mvaEGammaID.push_back(recoEle_mvaEGammaID[i]);
      tree_.recoEle_mvaEGammaValue.push_back(recoEle_mvaEGammaValue[i]);

      tree_.recoEle_mvaPFID.push_back(recoEle_mvaPFID[i]);
      tree_.recoEle_mvaPFValue.push_back(recoEle_mvaPFValue[i]);

      tree_.recoEle_mvaID.push_back(recoEle_mvaID[i]);
      tree_.recoEle_mvaValue.push_back(recoEle_mvaValue[i]);

    }


    //filling tree
    tree->Fill();
  
}


// ------------ method called once each job just before starting event loop  ------------
void
TQGenAnalyzer::beginJob()
{
  std::cout << "Starting job" << std::endl;

  // --- set up output tree
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("run", &tree_.run, "run/I");
  tree->Branch("event", &tree_.event, "event/I");
  tree->Branch("lumi", &tree_.lumi, "lumi/I");
  tree->Branch("nEle", &tree_.nEle, "nEle/I");
  tree->Branch("nMu", &tree_.nMu, "nMu/I");
  tree->Branch("vx", &tree_.vx, "vx/F");
  tree->Branch("vy", &tree_.vy, "vy/F");
  tree->Branch("vz", &tree_.vz, "vz/F");

  tree->Branch("TQ_genMass", &tree_.TQ_genMass, "TQ_genMass/F");
  tree->Branch("TQ_recoMass", &tree_.TQ_recoMass, "TQ_recoMass/F");
  tree->Branch("genLep_pt",     &tree_.genLep_pt);
  tree->Branch("genLep_mass",      &tree_.genLep_mass);
  tree->Branch("genLep_eta",&tree_.genLep_eta);
  tree->Branch("genLep_phi",&tree_.genLep_phi);
  tree->Branch("genLep_pdgId",&tree_.genLep_pdgId);
  tree->Branch("genMom_pt",     &tree_.genMom_pt);
  tree->Branch("genMom_mass",      &tree_.genMom_mass);
  tree->Branch("genMom_eta",&tree_.genMom_eta);
  tree->Branch("genMom_phi",&tree_.genMom_phi);
  tree->Branch("genMom_pdgId",&tree_.genMom_pdgId);

  tree->Branch("nMuReco", &tree_.nMuReco, "nMuReco/I");
  tree->Branch("recoMu_pt",     &tree_.recoMu_pt);
  tree->Branch("recoMu_mass",      &tree_.recoMu_mass);
  tree->Branch("recoMu_eta",&tree_.recoMu_eta);
  tree->Branch("recoMu_phi",&tree_.recoMu_phi);
  tree->Branch("recoMu_charge",&tree_.recoMu_charge);
  tree->Branch("recoMu_dR1",&tree_.recoMu_dR1);
  tree->Branch("recoMu_matchid1",&tree_.recoMu_matchid1);
  tree->Branch("recoMu_dR2",&tree_.recoMu_dR2);
  tree->Branch("recoMu_matchid2",&tree_.recoMu_matchid2);
  tree->Branch("recoMu_looseID",&tree_.recoMu_looseID);
  tree->Branch("recoMu_softID",&tree_.recoMu_softID);
  tree->Branch("recoMu_nTrkLayers",&tree_.recoMu_nTrkLayers);
  tree->Branch("recoMu_nPixLayers",&tree_.recoMu_nPixLayers);
  tree->Branch("recoMu_isHQ",&tree_.recoMu_isHQ);
  tree->Branch("recoMu_dxy",&tree_.recoMu_dxy);
  tree->Branch("recoMu_dz",&tree_.recoMu_dz);

  tree->Branch("nEleReco", &tree_.nEleReco, "nEleReco/I");
  tree->Branch("recoEle_pt",     &tree_.recoEle_pt);
  tree->Branch("recoEle_mass",      &tree_.recoEle_mass);
  tree->Branch("recoEle_eta",&tree_.recoEle_eta);
  tree->Branch("recoEle_phi",&tree_.recoEle_phi);
  tree->Branch("recoEle_charge",&tree_.recoEle_charge);

  tree->Branch("recoEle_vx",&tree_.recoEle_vx);
  tree->Branch("recoEle_vy",&tree_.recoEle_vy);
  tree->Branch("recoEle_vz",&tree_.recoEle_vz);

  tree->Branch("recoEle_dR1",&tree_.recoEle_dR1);
  tree->Branch("recoEle_matchid1",&tree_.recoEle_matchid1);
  tree->Branch("recoEle_dR2",&tree_.recoEle_dR2);
  tree->Branch("recoEle_matchid2",&tree_.recoEle_matchid2);
  tree->Branch("recoEle_E_ecal_preReg",&tree_.recoEle_E_ecal_preReg);
  tree->Branch("recoEle_E_ecal_postReg",&tree_.recoEle_E_ecal_postReg);
  tree->Branch("recoEle_E_ecaltrk_preReg",&tree_.recoEle_E_ecaltrk_preReg);
  tree->Branch("recoEle_E_ecaltrk_postReg",&tree_.recoEle_E_ecaltrk_postReg);

  tree->Branch("recoEle_rawE",&tree_.recoEle_rawE);
  tree->Branch("recoEle_corrEcalE", &tree_.recoEle_corrEcalE);
  tree->Branch("recoEle_gsfTrkChi2",&tree_.recoEle_gsfTrkChi2);
  tree->Branch("recoEle_passConvVeto",&tree_.recoEle_passConvVeto);
  tree->Branch("recoEle_isPF",&tree_.recoEle_isPF);

  tree->Branch("recoEle_isLowPt",&tree_.recoEle_isLowPt);

  tree->Branch("recoEle_isPFoverlap",&tree_.recoEle_isPFoverlap);

  tree->Branch("recoEle_mvaEGammaID",&tree_.recoEle_mvaEGammaID);
  tree->Branch("recoEle_mvaEGammaValue",&tree_.recoEle_mvaEGammaValue);

  tree->Branch("recoEle_mvaPFID",&tree_.recoEle_mvaPFID);
  tree->Branch("recoEle_mvaPFValue",&tree_.recoEle_mvaPFValue);

  tree->Branch("recoEle_mvaID",&tree_.recoEle_mvaID);
  tree->Branch("recoEle_mvaValue",&tree_.recoEle_mvaValue);
}

// ------------ method called once each job just after ending the event loop  ------------
void
TQGenAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TQGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


void TQGenAnalyzer::initTreeStructure()
{
  
  tree_.run= -500.;
  tree_.lumi= -500.;
  tree_.event= 0.;


}

void TQGenAnalyzer::clearVectors()
{

  tree_.genLep_pt.clear();
  tree_.genLep_mass.clear();
  tree_.genLep_eta.clear();
  tree_.genLep_phi.clear();
  tree_.genLep_pdgId.clear();
  tree_.genMom_pt.clear();
  tree_.genMom_mass.clear();
  tree_.genMom_eta.clear();
  tree_.genMom_phi.clear();
  tree_.genMom_pdgId.clear();

  tree_.recoMu_pt.clear();
  tree_.recoMu_mass.clear();
  tree_.recoMu_eta.clear();
  tree_.recoMu_phi.clear();
  tree_.recoMu_charge.clear();
  tree_.recoMu_dR1.clear();
  tree_.recoMu_dR2.clear();
  tree_.recoMu_matchid1.clear();
  tree_.recoMu_matchid2.clear();
  tree_.recoMu_looseID.clear();
  tree_.recoMu_softID.clear();
  tree_.recoMu_nTrkLayers.clear();
  tree_.recoMu_nPixLayers.clear();
  tree_.recoMu_isHQ.clear();
  tree_.recoMu_dxy.clear();
  tree_.recoMu_dz.clear();


  tree_.recoEle_pt.clear();
  tree_.recoEle_mass.clear();
  tree_.recoEle_eta.clear();
  tree_.recoEle_phi.clear();
  tree_.recoEle_charge.clear();

  tree_.recoEle_vx.clear();
  tree_.recoEle_vy.clear();
  tree_.recoEle_vz.clear();

  tree_.recoEle_dR1.clear();
  tree_.recoEle_dR2.clear();
  tree_.recoEle_matchid1.clear();
  tree_.recoEle_matchid2.clear();
  tree_.recoEle_E_ecal_preReg.clear();
  tree_.recoEle_E_ecal_postReg.clear();
  tree_.recoEle_E_ecaltrk_preReg.clear();
  tree_.recoEle_E_ecaltrk_postReg.clear();

  tree_.recoEle_rawE.clear();
  tree_.recoEle_corrEcalE.clear();
  tree_.recoEle_gsfTrkChi2.clear();
  tree_.recoEle_passConvVeto.clear();
  tree_.recoEle_isPF.clear();

  tree_.recoEle_isLowPt.clear();
  tree_.recoEle_mvaID.clear();
  tree_.recoEle_isPFoverlap.clear();

  tree_.recoEle_mvaEGammaID.clear();
  tree_.recoEle_mvaEGammaValue.clear();

  tree_.recoEle_mvaPFID.clear();
  tree_.recoEle_mvaPFValue.clear();
  
  tree_.recoEle_mvaID.clear();
  tree_.recoEle_mvaValue.clear();

}
//define this as a plug-in
DEFINE_FWK_MODULE(TQGenAnalyzer);
