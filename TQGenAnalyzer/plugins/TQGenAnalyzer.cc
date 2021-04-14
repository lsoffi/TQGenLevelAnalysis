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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


// root include files
#include "TTree.h"
#include "TLorentzVector.h"

#define MAX_PU_REWEIGHT 100

namespace reco { typedef edm::Ptr<GsfElectron> GsfElectronPtr; }


struct tree_struc_{

  int sampleID;
  float xsec;
  int nvtx;
  float rho;
  int npu;
  float puw_2016;
  float puw_2017;
  float puw_2018;
  float puw_ALL;
  int run;
  int lumi;
  long unsigned int event;
  int nEle;
  int nMu;
  float TQ_genMass;
  int triggerBit;
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

  void SetPuWeights(int year, double* puw_array);
  float GetPUWeight(int npu,int year);

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

  const edm::EDGetTokenT< double > rho_; // MINIAOD
  edm::Handle< double> rhoH_;

  const edm::EDGetTokenT<edm::View<PileupSummaryInfo> > PileUp_;
  edm:: Handle<edm::View< PileupSummaryInfo> > PileUpH_;

  const edm::EDGetTokenT< edm::TriggerResults> triggerBits_;
  edm::Handle<edm::TriggerResults> triggerBitsH_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerFlags_;
  edm::Handle<edm::TriggerResults> triggerFlagsH_;



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

  int sampleIndex_;
  float xsec_;    // pb

  double puweights2016_[100];
  double puweights2017_[100];
  double puweights2018_[100];
  double puweightsALL_[100];

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
  rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  rhoH_(),
  PileUp_(consumes<edm::View<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag> ("PileUp"))),
  PileUpH_(),
  triggerBits_( consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>( "bits" ))),
  triggerBitsH_(),
  triggerFlags_( consumes<edm::TriggerResults>( iConfig.getParameter<edm::InputTag>( "flags" ))),
  triggerFlagsH_(),
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
  

  // Event stuff
  sampleIndex_   = iConfig.getUntrackedParameter<int>("sampleIndex",0);
  xsec_          = iConfig.getUntrackedParameter<double>("sampleXsec",1.);



  if( iConfig.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
    const edm::ParameterSet& iconf = iConfig.getParameterSet("gsfRegressionConfig");
    const std::string& mname = iconf.getParameter<std::string>("modifierName");
    auto cc = consumesCollector();                  
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf, cc);
      regression_.reset(plugin);
      // edm::ConsumesCollector sumes = consumesCollector();
      // regression_->setConsumes(sumes);
  } else {
    regression_.reset(nullptr);
  }


  // Regression stuff - lowPtElectrons
  /*
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
  */

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
    iEvent.getByToken(patElectronsLowPt_,gsfElectronsLowPtH_);
    iEvent.getByToken(vtx_,vtxH_);
    iEvent.getByToken(rho_,rhoH_);
    iEvent.getByToken(PileUp_,PileUpH_);
    iEvent.getByToken( triggerBits_, triggerBitsH_ );
    iEvent.getByToken( triggerFlags_, triggerFlagsH_ );
    
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

    
    // --- sample info (0:signal, <0 background, >0 data)
    int sampleID = sampleIndex_;
    float xsec=1.;
    
    if(sampleID<=0) xsec=xsec_;

    // --- trigger info
    const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBitsH_ );
    //  vector<std::string> const &names = triggerNames.triggerNames();  
    int triggerBit=0;
    for( unsigned index = 0; index < triggerNames.size(); ++index ) {
      //here i print all trigger paths
      //  if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_") && triggerBitsH_->accept( index ) == 1) std::cout << index << " " << triggerNames.triggerName( index ) << " " << triggerBitsH_->accept( index ) << std::endl;
      //here i store trigger bit
      if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Physics") ) triggerBit=triggerBitsH_->accept( index );

    }

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
    float puw_2016 = 1.;
    float puw_2017 = 1.;
    float puw_2018 = 1.;
    float puw_ALL = 1.;
    int npu      = -1;
    
    if(sampleID <=0){
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


      npu = 0;
      for( unsigned int PVI = 0; PVI < PileUpH_->size(); ++PVI ) {
	Int_t pu_bunchcrossing = PileUpH_->ptrAt( PVI )->getBunchCrossing();
	if( pu_bunchcrossing == 0 ) {
	  npu = PileUpH_->ptrAt( PVI )->getTrueNumInteractions();
	}
      }

      puw_2016 = GetPUWeight(npu,2016);
      puw_2017 = GetPUWeight(npu,2017);
      puw_2018 = GetPUWeight(npu,2018);
      puw_ALL = GetPUWeight(npu,0);

    }

    //saving info of primary vertex
    const reco::Vertex & pv = vtxH_->front(); 
    vx=pv.x();
    vy=pv.y();
    vz=pv.z();
    int nvtx=vtxH_->size();

    //saving rho
    float rho = *(rhoH_.product());


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

      std::cout<<" mva: "<<mvaValueEGammaH_->size()<<"  gsf: "<<gsfElectronsH_->size()<<std::endl;
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
    tree_.nvtx=nvtx;
    tree_.rho=rho;
    tree_.npu=npu;
    tree_.puw_2016=puw_2016;
    tree_.puw_2017=puw_2017;
    tree_.puw_2018=puw_2018;
    tree_.puw_ALL=puw_ALL;
    tree_.triggerBit=triggerBit;
    tree_.sampleID=sampleID;
    tree_.xsec=xsec;

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
  SetPuWeights(2016,puweights2016_);
  SetPuWeights(2017,puweights2017_);
  SetPuWeights(2018,puweights2018_);
  SetPuWeights(0,puweightsALL_);
  // --- set up output tree
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("sampleID", &tree_.sampleID, "sampleID/I");
  tree->Branch("xsec", &tree_.xsec, "xsec/F");
  tree->Branch("run", &tree_.run, "run/I");
  tree->Branch("event", &tree_.event, "event/I");
  tree->Branch("lumi", &tree_.lumi, "lumi/I");
  tree->Branch("nEle", &tree_.nEle, "nEle/I");
  tree->Branch("nMu", &tree_.nMu, "nMu/I");
  tree->Branch("vx", &tree_.vx, "vx/F");
  tree->Branch("vy", &tree_.vy, "vy/F");
  tree->Branch("vz", &tree_.vz, "vz/F");
  tree->Branch("nvtx", &tree_.nvtx, "nvtx/I");
  tree->Branch("rho", &tree_.rho, "rho/F");
  tree->Branch("npu", &tree_.npu, "npu/I");
  tree->Branch("puw_2016", &tree_.puw_2016, "puw_2016/F");
  tree->Branch("puw_2017", &tree_.puw_2017, "puw_2017/F");
  tree->Branch("puw_2018", &tree_.puw_2018, "puw_2018/F");
  tree->Branch("puw_ALL", &tree_.puw_ALL, "puw_ALL/F");
  tree->Branch("triggerBit", &tree_.triggerBit, "triggerBit/I");

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

float TQGenAnalyzer::GetPUWeight(int npu,int year) {
  
  float weight=1;
  if (sampleIndex_<=0 && npu<100){
    if(year==2016)weight = puweights2016_[npu];
    if(year==2017)weight = puweights2017_[npu];
    if(year==2018)weight = puweights2018_[npu];
    if(year==0)weight = puweightsALL_[npu];

  }

  return weight;
}


void TQGenAnalyzer::SetPuWeights(int year,double* puw_array) {

  std::string puWeightFile;
  if(year==2016) puWeightFile="PU/pileup_2016.root";
  else   if(year==2017) puWeightFile="PU/pileup_2017.root";
  else   if(year==2018) puWeightFile="PU/pileup_2018.root";
  else   if(year==0) puWeightFile="PU/pileup_ALL.root";

  if (puWeightFile == "") {
    std::cout << "you need a weights file to use this function" << std::endl;
    return;
  }

  std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;

  TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");
  f_pu->cd();

  TH1D* pu_data = (TH1D*) f_pu->Get("pileup");
  pu_data->Scale(1./pu_data->Integral());

  double MC_weights_2016UL[100]={1.00402360149e-05, 5.76498797172e-05, 7.37891400294e-05, 0.000110932895295, 0.000158857714773,
				 0.000368637432599, 0.000893114107873, 0.00189700774575, 0.00358880167437, 0.00636052573486,
				 0.0104173961179, 0.0158122597405, 0.0223785660712, 0.0299186888073, 0.0380275944896,
				 0.0454313901624, 0.0511181088317, 0.0547434577348, 0.0567906239028, 0.0577145461461,
				 0.0578176902735, 0.0571251566494, 0.0555456541498, 0.053134383488, 0.0501519041462,
				 0.0466815838899, 0.0429244592524, 0.0389566776898, 0.0348507152776, 0.0307356862528,
				 0.0267712092206, 0.0229720184534, 0.0193388653099, 0.0159602510813, 0.0129310510552,
				 0.0102888654183, 0.00798782770975, 0.00606651703058, 0.00447820948367, 0.00321589786478,
				 0.0022450422045, 0.00151447388514, 0.000981183695515, 0.000609670479759, 0.000362193408119,
				 0.000211572646801, 0.000119152364744, 6.49133515399e-05, 3.57795801581e-05, 1.99043569043e-05,
				 1.13639319832e-05, 6.49624103579e-06, 3.96626216416e-06, 2.37910222874e-06, 1.50997403362e-06,
				 1.09816650247e-06, 7.31298519122e-07, 6.10398791529e-07, 3.74845774388e-07, 2.65177281359e-07,
				 2.01923536742e-07, 1.39347583555e-07, 8.32600052913e-08, 6.04932421298e-08, 6.52536630583e-08,
				 5.90574603808e-08, 2.29162474068e-08, 1.97294602668e-08, 1.7731096903e-08, 3.57547932012e-09,
				 1.35039815662e-09, 8.50071242076e-09, 5.0279187473e-09, 4.93736669066e-10, 8.13919708923e-10,
				 5.62778926097e-09, 5.15140589469e-10, 8.21676746568e-10, 0.0, 1.49166873577e-09,
				 8.43517992503e-09, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0,0.0};
  

  double MC_weights_2017UL[100]={1.1840841518e-05, 3.46661037703e-05, 8.98772521472e-05, 7.47400487733e-05, 0.000123005176624,
				 0.000156501700614, 0.000154660478659, 0.000177496185603, 0.000324149805611, 0.000737524009713,
				 0.00140432980253, 0.00244424508696, 0.00380027898037, 0.00541093042612, 0.00768803501793,
				 0.010828224552, 0.0146608623707, 0.01887739113, 0.0228418813823, 0.0264817796874,
				 0.0294637401336, 0.0317960986171, 0.0336645950831, 0.0352638818387, 0.036869429333,
				 0.0382797316998, 0.039386705577, 0.0398389681346, 0.039646211131, 0.0388392805703,
				 0.0374195678161, 0.0355377892706, 0.0333383902828, 0.0308286549265, 0.0282914440969,
				 0.0257860718304, 0.02341635055, 0.0213126338243, 0.0195035612803, 0.0181079838989,
				 0.0171991315458, 0.0166377598339, 0.0166445341361, 0.0171943735369, 0.0181980997278,
				 0.0191339792146, 0.0198518804356, 0.0199714909193, 0.0194616474094, 0.0178626975229,
				 0.0153296785464, 0.0126789254325, 0.0100766041988, 0.00773867100481, 0.00592386091874,
				 0.00434706240169, 0.00310217013427, 0.00213213401899, 0.0013996000761, 0.000879148859271,
				 0.000540866009427, 0.000326115560156, 0.000193965828516, 0.000114607606623, 6.74262828734e-05,
				 3.97805301078e-05, 2.19948704638e-05, 9.72007976207e-06, 4.26179259146e-06, 2.80015581327e-06,
				 1.14675436465e-06, 2.52452411995e-07, 9.08394910044e-08, 1.14291987912e-08, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0,0.0};

  double MC_weights_2018UL[100]={8.89374611122e-07, 1.1777062868e-05, 3.99725585118e-05, 0.000129888015252, 0.000265224848687,
				 0.000313088635109, 0.000353781668514, 0.000508787237162, 0.000873670065767, 0.00147166880932,
				 0.00228230649018, 0.00330375581273, 0.00466047608406, 0.00624959203029, 0.00810375867901,
				 0.010306521821, 0.0129512453978, 0.0160303925502, 0.0192913204592, 0.0223108613632,
				 0.0249798930986, 0.0273973789867, 0.0294402350483, 0.031029854302, 0.0324583524255,
				 0.0338264469857, 0.0351267479019, 0.0360320204259, 0.0367489568401, 0.0374133183052,
				 0.0380352633799, 0.0386200967002, 0.039124376968, 0.0394201612616, 0.0394673457109,
				 0.0391705388069, 0.0384758587461, 0.0372984548399, 0.0356497876549, 0.0334655175178,
				 0.030823567063, 0.0278340752408, 0.0246009685048, 0.0212676009273, 0.0180250593982,
				 0.0149129830776, 0.0120582333486, 0.00953400069415, 0.00738546929512, 0.00563442079939,
				 0.00422052915668, 0.00312446316347, 0.00228717533955, 0.00164064894334, 0.00118425084792,
				 0.000847785826565, 0.000603466454784, 0.000419347268964, 0.000291768785963, 0.000199761337863,
				 0.000136624574661, 9.46855200945e-05, 6.80243180179e-05, 4.94806013765e-05, 3.53122628249e-05,
				 2.556765786e-05, 1.75845711623e-05, 1.23828210848e-05, 9.31669724108e-06, 6.0713272037e-06,
				 3.95387384933e-06, 2.02760874107e-06, 1.22535149516e-06, 9.79612472109e-07, 7.61730246474e-07,
				 4.2748847738e-07, 2.41170461205e-07, 1.38701083552e-07, 3.37678010922e-08, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0, 0.0,
				 0.0, 0.0, 0.0, 0.0,0.0};

  double pum=-1.;
  double pud=-1.;
  double puw=1;
  for(int i=0;i<100;i++){
    pud=pu_data->GetBinContent(i+1);
    pum=MC_weights_2016UL[i];
    puw=1;

    if(pum!=0)puw=pud/pum;

    puw_array[i]=puw;
  }

  
  
  
}
//define this as a plug-in
DEFINE_FWK_MODULE(TQGenAnalyzer);
