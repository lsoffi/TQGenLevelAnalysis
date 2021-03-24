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

  std::vector<float>            recoEle_pt;
  std::vector<float>            recoEle_eta;
  std::vector<float>            recoEle_mass;
  std::vector<float>            recoEle_phi;
  std::vector<int>              recoEle_charge;
  std::vector<float>            recoEle_dR1;
  std::vector<float>            recoEle_dR2;
  std::vector<int>              recoEle_matchid1;
  std::vector<int>              recoEle_matchid2;
  std::vector<float>recoEle_E_ecal_preReg;
  std::vector<float>recoEle_E_ecal_postReg;
  std::vector<float>recoEle_E_ecaltrk_preReg;
  std::vector<float>recoEle_E_ecaltrk_postReg;



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

  
  // setup tree;
  TTree* tree;
  tree_struc_ tree_;


  std::unique_ptr<ModifyObjectValueBase> regression_;     // Low pt 
  
};


TQGenAnalyzer::TQGenAnalyzer(const edm::ParameterSet& iConfig)
  :
  prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  genParticlesH_(),
  patMuons_(consumes< edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("patMuons"))),
  patMuonsH_(),
  patElectrons_(consumes< edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("patElectrons"))),
  gsfElectronsH_()

{
  
  // Regression stuff - lowPtElectrons
  if( iConfig.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
    const edm::ParameterSet& iconf = iConfig.getParameterSet("gsfRegressionConfig");
    const std::string& mname = iconf.getParameter<std::string>("modifierName");
    std::cout<<mname<<std::endl;
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

    std::vector<float>recoEle_pt;
    std::vector<float>recoEle_mass;
    std::vector<float>recoEle_eta;
    std::vector<float>recoEle_phi;
    std::vector<float>recoEle_dR1;
    std::vector<float>recoEle_dR2;
    std::vector<int>              recoEle_charge;
    std::vector<int>              recoEle_matchid1;
    std::vector<int>              recoEle_matchid2;
    std::vector<float>recoEle_E_ecal_preReg;
    std::vector<float>recoEle_E_ecal_postReg;
    std::vector<float>recoEle_E_ecaltrk_preReg;
    std::vector<float>recoEle_E_ecaltrk_postReg;


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
    
    
    std::cout<<"-----------------------------------------------------"<<std::endl;

    // save gen particles (only leptons)
    
    for (const auto & genpar_iter : *genParticlesH_){ // loop over genparticles
      if(abs(genpar_iter.pdgId())!=13 && abs(genpar_iter.pdgId())!=11)continue;

      if(genpar_iter.mother(0)->pdgId()<20 || genpar_iter.mother(0)->pdgId()>200) continue;
      
      //      std::cout<<"status: "<<genpar_iter.status()<<" pdgid: "<<genpar_iter.pdgId()<<" pt: "<<genpar_iter.pt()<<" mass: "<<genpar_iter.mass()<<" eta: "<<genpar_iter.eta()<<" phi: "<<genpar_iter.phi()<<std::endl;

      
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

    //looping over muons
    for (const auto & mu_iter : *patMuonsH_){
      nMuReco++;

      float dR=999.;
      float dRMin1=999.;
      float dRMin2=999.;
      int recomatchid1=999;
      int recomatchid2=999;

      float recopt;
      float recoeta;
      float recophi;
      float recomass;
      int recocharge=999;


      float eta=mu_iter.eta();
      float phi=mu_iter.phi();
      

      recopt=mu_iter.pt();
      recoeta=mu_iter.eta();
      recophi=mu_iter.phi();
      recomass=mu_iter.mass();
      recocharge=mu_iter.charge();
      
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

    }





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
      

      
      recoEle_pt.push_back(recopt);
      recoEle_mass.push_back(recomass);
      recoEle_eta.push_back(recoeta);
      recoEle_phi.push_back(recophi);
      recoEle_charge.push_back(recocharge);

      recoEle_dR1.push_back(dRMin1);
      recoEle_dR2.push_back(dRMin2);
      recoEle_matchid1.push_back(recomatchid1);
      recoEle_matchid2.push_back(recomatchid2);

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
    
    std::cout<<" pre ecal: "<<recopre_ecal<<" pre ecaltrk: "<<recopre_ecaltrk<<std::endl;

    //    regression_->modifyObject(newElectron);
    recopost_ecal = newElectron.correctedEcalEnergy();
    recopost_ecaltrk = newElectron.energy();

    std::cout<<" post ecal: "<<recopost_ecal<<" post ecaltrk: "<<recopost_ecaltrk<<std::endl;

    recoEle_E_ecal_preReg.push_back(recopre_ecal);
    recoEle_E_ecal_postReg.push_back(recopost_ecal);
    recoEle_E_ecaltrk_preReg.push_back(recopre_ecaltrk);
    recoEle_E_ecaltrk_postReg.push_back(recopost_ecaltrk);

}
 



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
      
      //      if(i==0)recolep0->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);
      // if(i==1)recolep1->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);

    }



  tree_.nEleReco=nEleReco;
    for(unsigned int i=0;i<recoEle_pt.size();i++){
      tree_.recoEle_pt.push_back(recoEle_pt[i]);
      tree_.recoEle_mass.push_back(recoEle_mass[i]);
      tree_.recoEle_eta.push_back(recoEle_eta[i]);
      tree_.recoEle_phi.push_back(recoEle_phi[i]);   
      tree_.recoEle_charge.push_back(recoEle_charge[i]);   

      tree_.recoEle_dR1.push_back(recoEle_dR1[i]);   
      tree_.recoEle_matchid1.push_back(recoEle_matchid1[i]);   
      tree_.recoEle_dR2.push_back(recoEle_dR2[i]);   
      tree_.recoEle_matchid2.push_back(recoEle_matchid2[i]);   
      tree_.recoEle_E_ecal_preReg.push_back(recoEle_E_ecal_preReg[i]);   
      tree_.recoEle_E_ecal_postReg.push_back(recoEle_E_ecal_postReg[i]);   
      tree_.recoEle_E_ecaltrk_preReg.push_back(recoEle_E_ecaltrk_preReg[i]);   
      tree_.recoEle_E_ecaltrk_postReg.push_back(recoEle_E_ecaltrk_postReg[i]);   
      
      //      if(i==0)recolep0->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);
      // if(i==1)recolep1->SetPtEtaPhiM(recoMu_pt[i],recoMu_eta[i],recoMu_phi[i],recoMu_mass[i]);

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

  tree->Branch("nEleReco", &tree_.nEleReco, "nEleReco/I");
  tree->Branch("recoEle_pt",     &tree_.recoEle_pt);
  tree->Branch("recoEle_mass",      &tree_.recoEle_mass);
  tree->Branch("recoEle_eta",&tree_.recoEle_eta);
  tree->Branch("recoEle_phi",&tree_.recoEle_phi);
  tree->Branch("recoEle_charge",&tree_.recoEle_charge);
  tree->Branch("recoEle_dR1",&tree_.recoEle_dR1);
  tree->Branch("recoEle_matchid1",&tree_.recoEle_matchid1);
  tree->Branch("recoEle_dR2",&tree_.recoEle_dR2);
  tree->Branch("recoEle_matchid2",&tree_.recoEle_matchid2);
  tree->Branch("recoEle_E_ecal_preReg",&tree_.recoEle_E_ecal_preReg);
  tree->Branch("recoEle_E_ecal_postReg",&tree_.recoEle_E_ecal_postReg);
  tree->Branch("recoEle_E_ecaltrk_preReg",&tree_.recoEle_E_ecaltrk_preReg);
  tree->Branch("recoEle_E_ecaltrk_postReg",&tree_.recoEle_E_ecaltrk_postReg);
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


  tree_.recoEle_pt.clear();
  tree_.recoEle_mass.clear();
  tree_.recoEle_eta.clear();
  tree_.recoEle_phi.clear();
  tree_.recoEle_charge.clear();
  tree_.recoEle_dR1.clear();
  tree_.recoEle_dR2.clear();
  tree_.recoEle_matchid1.clear();
  tree_.recoEle_matchid2.clear();
  tree_.recoEle_E_ecal_preReg.clear();
  tree_.recoEle_E_ecal_postReg.clear();
  tree_.recoEle_E_ecaltrk_preReg.clear();
  tree_.recoEle_E_ecaltrk_postReg.clear();



}
//define this as a plug-in
DEFINE_FWK_MODULE(TQGenAnalyzer);
