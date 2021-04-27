// -*-C++ -*-
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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
// #include "RecoVertex/KinematicFitPrimitives/interface/"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
// root include files
#include "TTree.h"
#include "TLorentzVector.h"

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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/Associations/interface/MuonToTrackingParticleAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "RecoEgamma/ElectronIdentification/interface/ElectronIDAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "KinVtxFitter.h"

// root include files
#include "TTree.h"
#include "TLorentzVector.h"



#define MAX_PU_REWEIGHT 100
#define LEP_SIGMA 0.0000001

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


  int nMuReco;
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


  //dimuon pairs
  int nDimuReco;
  std::vector<float>            recoDimu_vx;
  std::vector<float>            recoDimu_vy;
  std::vector<float>            recoDimu_vz;
  std::vector<float>            recoDimu_vtxchi2;
  std::vector<float>            recoDimu_vtxndof;
  std::vector<float>            recoDimu_vtxprob;
  std::vector<float>            recoDimu_pt1;  
  std::vector<float>            recoDimu_eta1;
  std::vector<float>            recoDimu_phi1;
  std::vector<float>            recoDimu_charge1;
  std::vector<float>            recoDimu_mass1;        
  std::vector<float>            recoDimu_pt2;  
  std::vector<float>            recoDimu_eta2;
  std::vector<float>            recoDimu_phi2;
  std::vector<float>            recoDimu_charge2;
  std::vector<float>            recoDimu_mass2;        
  std::vector<float>            recoDimu_pt;        
  std::vector<float>            recoDimu_eta;        
  std::vector<float>            recoDimu_phi;        
  std::vector<float>            recoDimu_mass;        
  std::vector<float>            recoDimu_massErr;        

  //electrons
  int nEleReco;
  int nPFEleReco;
  int nLowPtEleReco;
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
  std::vector<float>recoEle_mvaPFValue;
  std::vector<float>recoEle_isPFoverlap;
  std::vector<float>recoEle_isLowPt;
  std::vector<float>recoEle_mvaValue;

  //dielectron pairs
  int nDieleReco;
  std::vector<float>            recoDiele_vx;
  std::vector<float>            recoDiele_vy;
  std::vector<float>            recoDiele_vz;
  std::vector<float>            recoDiele_vtxchi2;
  std::vector<float>            recoDiele_vtxndof;
  std::vector<float>            recoDiele_vtxprob;
  std::vector<float>            recoDiele_pt1;  
  std::vector<float>            recoDiele_eta1;
  std::vector<float>            recoDiele_phi1;
  std::vector<float>            recoDiele_charge1;
  std::vector<float>            recoDiele_mass1;        
  std::vector<int>            recoDiele_isPF1;        
  std::vector<int>            recoDiele_isLowPt1;        
  std::vector<int>            recoDiele_isPFoverlap1;        
  std::vector<float>            recoDiele_pt2;  
  std::vector<float>            recoDiele_eta2;
  std::vector<float>            recoDiele_phi2;
  std::vector<float>            recoDiele_charge2;
  std::vector<float>            recoDiele_mass2;        
  std::vector<int>            recoDiele_isPF2;        
  std::vector<int>            recoDiele_isLowPt2;        
  std::vector<int>            recoDiele_isPFoverlap2;        
  std::vector<float>            recoDiele_pt;        
  std::vector<float>            recoDiele_eta;        
  std::vector<float>            recoDiele_phi;        
  std::vector<float>            recoDiele_mass;        
  std::vector<float>            recoDiele_massErr;        

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

  // ---------- member data -------------------- //
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<GenEventInfoProduct>genInfoToken_;
  
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > prunedGenParticles_; // MINIAOD
  edm::Handle< edm::View<reco::GenParticle> > genParticlesH_;

  const edm::EDGetTokenT< std::vector<pat::Muon> > patMuons_; // MINIAOD
  edm::Handle< std::vector<pat::Muon> > patMuonsH_;

  const edm::EDGetTokenT<pat::ElectronCollection> lowpt_src_;
   const edm::EDGetTokenT<pat::ElectronCollection> pf_src_;

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

   //setup mva ID for PF electrons

   const edm::EDGetTokenT< edm::ValueMap<float> > mvaValuePF_; 
   edm::Handle< edm::ValueMap<float> > mvaValuePFH_;


   int sampleIndex_;
   float xsec_;    // pb

   double puweights2016_[100];
   double puweights2017_[100];
   double puweights2018_[100];
   double puweightsALL_[100];

   // setup tree;
   TTree* tree;
   tree_struc_ tree_;


};


 TQGenAnalyzer::TQGenAnalyzer(const edm::ParameterSet& iConfig)
   :
   prunedGenParticles_(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
   genParticlesH_(),
   patMuons_(consumes< std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("patMuons"))),
   patMuonsH_(),
   lowpt_src_{ consumes<pat::ElectronCollection>( iConfig.getParameter<edm::InputTag>("lowptSrc") )},
   pf_src_{ consumes<pat::ElectronCollection>( iConfig.getParameter<edm::InputTag>("pfSrc") )},
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
   mvaValuePF_(consumes< edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuePF"))),
   mvaValuePFH_()


 {


   // Event stuff
   sampleIndex_   = iConfig.getUntrackedParameter<int>("sampleIndex",0);
   xsec_          = iConfig.getUntrackedParameter<double>("sampleXsec",1.);

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
     edm::Handle<pat::ElectronCollection> lowpt;
     iEvent.getByToken(lowpt_src_, lowpt);
     edm::Handle<pat::ElectronCollection> pf;
     iEvent.getByToken(pf_src_, pf);
     iEvent.getByToken(vtx_,vtxH_);
     iEvent.getByToken(rho_,rhoH_);
     iEvent.getByToken(PileUp_,PileUpH_);
     iEvent.getByToken( triggerBits_, triggerBitsH_ );
     iEvent.getByToken( triggerFlags_, triggerFlagsH_ );
     iEvent.getByToken(mvaValuePF_, mvaValuePFH_); 

     edm::ESHandle<MagneticField> bFieldHandle;
     iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
     edm::ESHandle<TransientTrackBuilder> theB ;
     iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


     using namespace edm;



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

     std::vector<float>            recoDimu_vx;
     std::vector<float>            recoDimu_vy;
     std::vector<float>            recoDimu_vz;
     std::vector<float>            recoDimu_vtxchi2;
     std::vector<float>            recoDimu_vtxndof;
     std::vector<float>            recoDimu_vtxprob;
     std::vector<float>            recoDimu_pt1;  
     std::vector<float>            recoDimu_eta1;
     std::vector<float>            recoDimu_phi1;
     std::vector<float>            recoDimu_charge1;
     std::vector<float>            recoDimu_mass1;        
     std::vector<float>            recoDimu_pt2;  
     std::vector<float>            recoDimu_eta2;
     std::vector<float>            recoDimu_phi2;
     std::vector<float>            recoDimu_charge2;
     std::vector<float>            recoDimu_mass2;        
     std::vector<float>            recoDimu_pt;        
     std::vector<float>            recoDimu_eta;        
     std::vector<float>            recoDimu_phi;        
     std::vector<float>            recoDimu_mass;        
     std::vector<float>            recoDimu_massErr;        


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
     std::vector<float>recoEle_mvaPFValue;
     std::vector<float>recoEle_mvaValue;
     std::vector<float>recoEle_isPFoverlap;
     std::vector<float>recoEle_isLowPt;
     


     std::vector<float>            recoDiele_vx;
     std::vector<float>            recoDiele_vy;
     std::vector<float>            recoDiele_vz;
     std::vector<float>            recoDiele_vtxchi2;
     std::vector<float>            recoDiele_vtxndof;
     std::vector<float>            recoDiele_vtxprob;
     std::vector<float>            recoDiele_pt1;  
     std::vector<float>            recoDiele_eta1;
     std::vector<float>            recoDiele_phi1;
     std::vector<float>            recoDiele_charge1;
     std::vector<float>            recoDiele_mass1;        
     std::vector<float>            recoDiele_pt2;  
     std::vector<float>            recoDiele_eta2;
     std::vector<float>            recoDiele_phi2;
     std::vector<float>            recoDiele_charge2;
     std::vector<float>            recoDiele_mass2;        
     std::vector<int>            recoDiele_isPF1;        
     std::vector<int>            recoDiele_isLowPt1;        
     std::vector<int>            recoDiele_isPFoverlap1;        
     std::vector<int>            recoDiele_isPF2;        
     std::vector<int>            recoDiele_isLowPt2;        
     std::vector<int>            recoDiele_isPFoverlap2;        
     std::vector<float>            recoDiele_pt;        
     std::vector<float>            recoDiele_eta;        
     std::vector<float>            recoDiele_phi;        
     std::vector<float>            recoDiele_mass;        
     std::vector<float>            recoDiele_massErr;        

     // --- sample info (0:signal, <0 background, >0 data)
     int sampleID = sampleIndex_;
     float xsec=1.;
  
     if(sampleID<=0) xsec=xsec_;

     // ---------- TRIGGER -------------------- //
     const edm::TriggerNames &triggerNames = iEvent.triggerNames( *triggerBitsH_ );

     //  vector<std::string> const &names = triggerNames.triggerNames();  
     int triggerBit=0;
     for( unsigned index = 0; index < triggerNames.size(); ++index ) {
      
       //  if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_") && triggerBitsH_->accept( index ) == 1) std::cout << index << " " << triggerNames.triggerName( index ) << " " << triggerBitsH_->accept( index ) << std::endl;
      
       if( (TString::Format((triggerNames.triggerName( index )).c_str())).Contains("HLT_Physics") ) triggerBit=triggerBitsH_->accept( index );

     }

     // ---------- GENERAL INFOS -------------------- //
     unsigned long int event = iEvent.id().event();   
     int run                 = iEvent.id().run();
     int lumi                = iEvent.luminosityBlock();
     int nEle=0;
     int nMu=0;
     int nMuReco=0;
     int nEleReco=0;
     int nPFEleReco=0;
     int nLowPtEleReco=0;
     int nDimuReco=0;
     int nDieleReco=0;

     TLorentzVector* lep1=new TLorentzVector();
     TLorentzVector* lep2=new TLorentzVector();
     TLorentzVector* lep3=new TLorentzVector();
     TLorentzVector* lep0=new TLorentzVector();
  
     float vx;
     float vy;
     float vz;


     std::cout<<"-----------------------------------------------------"<<std::endl;
	
     
     float puw_2016 = 1.;
     float puw_2017 = 1.;
     float puw_2018 = 1.;
     float puw_ALL = 1.;
     int npu      = -1;
  
     // ---------- GEN LEVEL INFO -------------------- //
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


     // ---------- VERTICES -------------------- //
     const reco::Vertex & pv = vtxH_->front(); 
     vx=pv.x();
     vy=pv.y();
     vz=pv.z();
     int nvtx=vtxH_->size();

     //saving rho
     float rho = *(rhoH_.product());



     // ---------- MUONS -------------------- //
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
       //recontrklayers=mu_iter.innerTrack()->hitPattern().trackerLayersWithMeasurement();
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


     // ---------- DIMUONS -------------------- //			
     size_t imu1=-1;
     size_t imu2=-1;

     std::cout<<" nMuReco:   "<<nMuReco<<std::endl;
     //     for (pat::MuonCollection::const_iterator mu1 = patMuonsH_.begin(); mu1 != patMuonsH_.end(); ++mu1){
     for (const pat::Muon & mu1 : *patMuonsH_){
       imu1++;
       if(mu1.pt()<1.)continue;
       if(abs(mu1.eta())>2.5)continue;

       imu2=-1;
       // for (pat::MuonCollection::const_iterator mu2 = patMuonsH_.begin(); mu2 != patMuonsH_.end(); ++mu2){
       for (const pat::Muon & mu2 : *patMuonsH_){
	 imu2++;
	 if(imu2<=imu1)continue;
	 if(mu2.pt()<1.)continue;
	 if(abs(mu2.eta())>2.5)continue;
	 if(mu1.charge()*mu2.charge()>0)continue;
    
	 std::cout<<"-----------------------------------------"<<std::endl;
	 std::cout<<" pt1: "<<mu1.pt()<<" pt2: "<<mu2.pt()<<" eta1: "<<mu1.eta()<<" eta2: "<<mu2.eta()<<" phi1: "<<mu1.phi()<<" phi2: "<<mu2.phi()<<" ch1: "<<mu1.charge()<<" ch2: "<<mu2.charge()<<std::endl;
	 

       recoDimu_pt1.push_back(mu1.pt());
       recoDimu_eta1.push_back(mu1.eta());
       recoDimu_phi1.push_back(mu1.phi());
       recoDimu_charge1.push_back(mu1.charge());
       recoDimu_mass1.push_back(mu1.mass());


       recoDimu_pt2.push_back(mu2.pt());
       recoDimu_eta2.push_back(mu2.eta());
       recoDimu_phi2.push_back(mu2.phi());
       recoDimu_charge2.push_back(mu2.charge());
       recoDimu_mass2.push_back(mu2.mass());
       
	
       //run kinematic fit
	
       reco::TransientTrack mu1TT = theB->build(mu1.bestTrack());
       reco::TransientTrack mu2TT = theB->build(mu2.bestTrack());
	/*	const reco::TransientTrack mu1TT((*(mu1.bestTrack())),&(*bFieldHandle)); 
	const reco::TransientTrack mu2TT((*(mu2.bestTrack())),&(*bFieldHandle)); 
	*/
       //       std::cout<<mu1TT.isValid()<<" "<<mu2TT.isValid()<<std::endl;
       //       std::cout<<mu1TT.numberOfValidHits()<<" "<<mu2TT.numberOfValidHits()<<std::endl;

	float chi = 0.;
	float ndf = 0.;
	double muon_mass = 0.1056583;
	float muon_sigma = 0.0000001;

	
	KinematicParticleFactoryFromTransientTrack pFactory;
	std::vector<RefCountedKinematicParticle> allParticles;
	allParticles.push_back(pFactory.particle (mu1TT,muon_mass,chi,ndf,muon_sigma));
	allParticles.push_back(pFactory.particle (mu2TT,muon_mass,chi,ndf,muon_sigma));
	KinematicParticleVertexFitter fitter;
	
	/*RefCountedKinematicTree vtx_tree = fitter.fit(allParticles);
	if(vtx_tree->isEmpty() == 1)
	  std::cout << "Kinematic Fit unsuccesfull" << std::endl;
	else{
	vtx_tree->movePointerToTheTop(); 
	RefCountedKinematicVertex fitted_vtx_ = vtx_tree->currentDecayVertex();

	float fitted_chi2=fitted_vtx_->chiSquared();
	std::cout<<" chi2: "<<fitted_chi2<<std::endl;
	}*/


	KinVtxFitter fitter2(
			   {mu1TT,mu2TT},
			   {muon_mass,muon_mass},
			   {muon_sigma,muon_sigma} //some small sigma for the lepton mass
		      );
	
       if(!fitter2.success()) continue; 
       
       std::cout<<" chi2: "<<fitter2.chi2()<<std::endl;
       std::cout<<" dof: "<<fitter2.dof()<<std::endl;
       std::cout<<" prob: "<<fitter2.prob()<<std::endl;
       recoDimu_vx.push_back(fitter2.fitted_vtx().x());
       recoDimu_vy.push_back(fitter2.fitted_vtx().y());
       recoDimu_vz.push_back(fitter2.fitted_vtx().z());
       recoDimu_vtxchi2.push_back(fitter2.chi2());
       recoDimu_vtxndof.push_back(fitter2.dof());
       recoDimu_vtxprob.push_back(fitter2.prob());
       /*       auto fit_p4 = fitter2.fitted_p4();
       recoDimu_pt.push_back(fit_p4.pt());        
       recoDimu_eta.push_back(fit_p4.eta());        
       recoDimu_phi.push_back(fit_p4.phi());        
       //       recoDimu_mass.push_back(fitter2.fitted_candidate().mass());        
       //recoDimu_massErr.push_back(sqrt(fitter2.fitted_candidate().kinematicParametersError().matrix()(6,6)));        
       */

      nDimuReco++;



    }
 }


     // ---------- PF ELECTRONS -------------------- //
    size_t ipfele=-1;
    for(auto ele : *pf) {
      ipfele++;
      nEleReco++;
      nPFEleReco++;
      float dR=999.;
      float dRMin1=999.;
      float dRMin2=999.;
      int recomatchid1=999;
      int recomatchid2=999;

      float recopt=ele.pt();
      float recoeta=ele.eta();
      float recophi=ele.phi();
      float recomass=ele.mass();
      int recocharge=ele.charge();

      float recovx=ele.vx();
      float recovy=ele.vy();
      float recovz=ele.vz();

      
      float recorawe=ele.superCluster()->rawEnergy();
      float recocorrecale=ele.correctedEcalEnergy();
      float recotrkchi2=ele.gsfTrack()->normalizedChi2();
      int recoconvveto=ele.passConversionVeto();
      

      
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


      
      float mvaPF_value = -999.;
      edm::Ref<pat::ElectronCollection> ref(pf,ipfele);
      mvaPF_value = float((*mvaValuePFH_)[ref]);
      
      float mva_value = -999.;

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

      recoEle_mvaPFValue.push_back(mvaPF_value);
      recoEle_mvaValue.push_back(mva_value);

      
      
      float recopre_ecal     = -999.;
      float recopre_ecaltrk  = -999.;
      float recopost_ecal    = -999.;
      float recopost_ecaltrk = -999.;
      
      recopre_ecal = ele.correctedEcalEnergy();
      recopre_ecaltrk = ele.energy();
      recopost_ecal = ele.correctedEcalEnergy();
      recopost_ecaltrk = ele.energy();
      
      recoEle_E_ecal_preReg.push_back(recopre_ecal);
      recoEle_E_ecal_postReg.push_back(recopost_ecal);
      recoEle_E_ecaltrk_preReg.push_back(recopre_ecaltrk);
      recoEle_E_ecaltrk_postReg.push_back(recopost_ecaltrk);

    }
    
    
    unsigned int nPFEle=recoEle_pt.size();

     // ---------- LOW PT ELECTRONS -------------------- //
     if(lowpt.isValid()){
       
       for(auto ele : *lowpt) {

	float dR=999.;
	float dRMin1=999.;
	float dRMin2=999.;
	int recomatchid1=999;
	int recomatchid2=999;

	float recopt=ele.pt();
	float recoeta=ele.eta();
	float recophi=ele.phi();
	float recomass=ele.mass();
	int recocharge=ele.charge();

	float recovx=ele.vx();
	float recovy=ele.vy();
	float recovz=ele.vz();
	
	int recoispfoverlap=0;

	// ---------- REMOVE OVERLAP WITH PF -------------------- //
	   
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
 
      if(recoispfoverlap==0) nEleReco++;
      if(recoispfoverlap==0) nLowPtEleReco++;
      
      float recorawe=ele.superCluster()->rawEnergy();
      float recocorrecale=ele.correctedEcalEnergy();
      float recotrkchi2=ele.gsfTrack()->normalizedChi2();
      int recoconvveto=ele.passConversionVeto();
      

      
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



      
      float mvaPF_value = -999.;      
      float mva_value = -999.;
      mva_value=ele.isElectronIDAvailable("ID") ? ele.electronID("ID") : -100. ;
      
      
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

      recoEle_mvaPFValue.push_back(mvaPF_value);

      recoEle_mvaValue.push_back(mva_value);


      // ---------- REGRESSION: TO BE APPLIED -------------------- //
    float recopre_ecal     = -999.;
    float recopre_ecaltrk  = -999.;
    float recopost_ecal    = -999.;
    float recopost_ecaltrk = -999.;


    recopre_ecal = ele.correctedEcalEnergy();
    recopre_ecaltrk = ele.energy();
    recopost_ecal = ele.correctedEcalEnergy();
    recopost_ecaltrk = ele.energy();

    recoEle_E_ecal_preReg.push_back(recopre_ecal);
    recoEle_E_ecal_postReg.push_back(recopost_ecal);
    recoEle_E_ecaltrk_preReg.push_back(recopre_ecaltrk);
    recoEle_E_ecaltrk_postReg.push_back(recopost_ecaltrk);

       }
 

    }


     // ----------------- DIELECTRON PF-PF -------------------- //
     size_t iel1=-1;
     size_t iel2=-1;
     for (const pat::Electron & el1 : *pf){
       iel1++;
       if(el1.pt()<1.)continue;
       if(abs(el1.eta())>2.5)continue;

       iel2=-1;
       // for (pat::MuonCollection::const_iterator mu2 = patMuonsH_.begin(); mu2 != patMuonsH_.end(); ++mu2){
       for (const pat::Electron & el2 : *pf){
	 iel2++;
	 if(iel2<=iel1)continue;
	 if(el2.pt()<1.)continue;
	 if(abs(el2.eta())>2.5)continue;
	 if(el1.charge()*el2.charge()>0)continue;
    

	 //	std::cout<<" pt1: "<<el1.pt()<<" pt2: "<<el2.pt()<<" eta1: "<<el1.eta()<<" eta2: "<<el2.eta()<<" phi1: "<<el1.phi()<<" phi2: "<<el2.phi()<<" ch1: "<<el1.charge()<<" ch2: "<<el2.charge()<<std::endl;
	

       recoDiele_pt1.push_back(el1.pt());
       recoDiele_eta1.push_back(el1.eta());
       recoDiele_phi1.push_back(el1.phi());
       recoDiele_charge1.push_back(el1.charge());
       recoDiele_mass1.push_back(el1.mass());
       recoDiele_isPF1.push_back(1);
       recoDiele_isLowPt1.push_back(0);
       recoDiele_isPFoverlap1.push_back(0);


       recoDiele_pt2.push_back(el2.pt());
       recoDiele_eta2.push_back(el2.eta());
       recoDiele_phi2.push_back(el2.phi());
       recoDiele_charge2.push_back(el2.charge());
       recoDiele_mass2.push_back(el2.mass());
       recoDiele_isPF2.push_back(1);
       recoDiele_isLowPt2.push_back(0);
       recoDiele_isPFoverlap2.push_back(0);
	
	
       //run kinematic fit
	
       reco::TransientTrack el1TT = (*theB).build( el1.gsfTrack() );
       reco::TransientTrack el2TT = (*theB).build( el2.gsfTrack() );

	float chi = 0.;
	float ndf = 0.;
	double elon_mass = 0.005056583;
	float elon_sigma = 0.0000001;

	KinVtxFitter fitter(
			   {el1TT,el2TT},
			   {elon_mass,elon_mass},
			   {elon_sigma,elon_sigma} //some small sigma for the lepton mass
		      );
	
       if(!fitter.success()) continue; 
       //       std::cout<<" prob: "<<fitter.prob()<<std::endl;
       recoDiele_vx.push_back(fitter.fitted_vtx().x());
       recoDiele_vy.push_back(fitter.fitted_vtx().y());
       recoDiele_vz.push_back(fitter.fitted_vtx().z());
       recoDiele_vtxchi2.push_back(fitter.chi2());
       recoDiele_vtxndof.push_back(fitter.dof());
       recoDiele_vtxprob.push_back(fitter.prob());
       
       nDieleReco++;



    }
 }



     // ----------------- DLP_IELECTRON LP-LP -------------------- //
     size_t lp_iel1=-1;
     size_t lp_iel2=-1;
     for (const pat::Electron & el1 : *lowpt){
       lp_iel1++;
       if(el1.pt()<1.)continue;
       if(abs(el1.eta())>2.5)continue;

       lp_iel2=-1;
       // for (pat::MuonCollection::const_iterator mu2 = patMuonsH_.begin(); mu2 != patMuonsH_.end(); ++mu2){
       for (const pat::Electron & el2 : *pf){
	 lp_iel2++;
	 if(lp_iel2<=lp_iel1)continue;
	 if(el2.pt()<1.)continue;
	 if(abs(el2.eta())>2.5)continue;
	 if(el1.charge()*el2.charge()>0)continue;
    

	 //	std::cout<<" pt1: "<<el1.pt()<<" pt2: "<<el2.pt()<<" eta1: "<<el1.eta()<<" eta2: "<<el2.eta()<<" phi1: "<<el1.phi()<<" phi2: "<<el2.phi()<<" ch1: "<<el1.charge()<<" ch2: "<<el2.charge()<<std::endl;
	bool clean_out1 = false;
	bool clean_out2 = false;
	for(unsigned int iEle=0; iEle<nPFEle; ++iEle) {
	  
	  clean_out1 |= (
			 fabs(recoEle_vz[iEle] - el1.vz()) < dz_cleaning_ &&
			 reco::deltaR(el1.eta(), el1.phi(), recoEle_eta[iEle], recoEle_phi[iEle]) < dr_cleaning_   );
	  clean_out2 |= (
			 fabs(recoEle_vz[iEle] - el2.vz()) < dz_cleaning_ &&
			 reco::deltaR(el2.eta(), el2.phi(), recoEle_eta[iEle], recoEle_phi[iEle]) < dr_cleaning_   );
	  
	}
	
            
      int recoisoverlap1=0;
      if(clean_out1) recoisoverlap1=1;
      else recoisoverlap1=0;
      int recoisoverlap2=0;
      if(clean_out2) recoisoverlap2=1;
      else recoisoverlap2=0;
	

       recoDiele_pt1.push_back(el1.pt());
       recoDiele_eta1.push_back(el1.eta());
       recoDiele_phi1.push_back(el1.phi());
       recoDiele_charge1.push_back(el1.charge());
       recoDiele_mass1.push_back(el1.mass());
       recoDiele_isPF1.push_back(0);
       recoDiele_isLowPt1.push_back(1);
       recoDiele_isPFoverlap1.push_back(recoisoverlap1);


       recoDiele_pt2.push_back(el2.pt());
       recoDiele_eta2.push_back(el2.eta());
       recoDiele_phi2.push_back(el2.phi());
       recoDiele_charge2.push_back(el2.charge());
       recoDiele_mass2.push_back(el2.mass());
       recoDiele_isPF2.push_back(0);
       recoDiele_isLowPt2.push_back(1);
       recoDiele_isPFoverlap2.push_back(recoisoverlap2);
	
	
       //run kinematic fit
	
       reco::TransientTrack el1TT = (*theB).build( el1.gsfTrack() );
       reco::TransientTrack el2TT = (*theB).build( el2.gsfTrack() );

	float chi = 0.;
	float ndf = 0.;
	double elon_mass = 0.005056583;
	float elon_sigma = 0.0000001;

	KinVtxFitter fitter(
			   {el1TT,el2TT},
			   {elon_mass,elon_mass},
			   {elon_sigma,elon_sigma} //some small sigma for the lepton mass
		      );
	
       if(!fitter.success()) continue; 
       //       std::cout<<" prob: "<<fitter.prob()<<std::endl;
       recoDiele_vx.push_back(fitter.fitted_vtx().x());
       recoDiele_vy.push_back(fitter.fitted_vtx().y());
       recoDiele_vz.push_back(fitter.fitted_vtx().z());
       recoDiele_vtxchi2.push_back(fitter.chi2());
       recoDiele_vtxndof.push_back(fitter.dof());
       recoDiele_vtxprob.push_back(fitter.prob());
       
       nDieleReco++;



    }
 }





    
    initTreeStructure();
    clearVectors();


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
      
    }


    
    tree_.nDimuReco=nDimuReco;
    for(unsigned int i=0;i<recoDimu_vx.size();i++){
       
      tree_.recoDimu_vx.push_back(           recoDimu_vx[i]);	      
      tree_.recoDimu_vy.push_back(	     recoDimu_vy[i]);	      
      tree_.recoDimu_vz.push_back(	     recoDimu_vz[i]);	      
      tree_.recoDimu_vtxchi2.push_back(     recoDimu_vtxchi2[i]);  
      tree_.recoDimu_vtxndof.push_back(     recoDimu_vtxndof[i]);  
      tree_.recoDimu_vtxprob.push_back(     recoDimu_vtxprob[i]);  
      tree_.recoDimu_pt1.push_back(  	     recoDimu_pt1[i]);      
      tree_.recoDimu_eta1.push_back(	     recoDimu_eta1[i]);     
      tree_.recoDimu_phi1.push_back(	     recoDimu_phi1[i]);     
      tree_.recoDimu_charge1.push_back(     recoDimu_charge1[i]);  
      tree_.recoDimu_mass1.push_back(       recoDimu_mass1[i]);     
      tree_.recoDimu_pt2.push_back(  	     recoDimu_pt2[i]);      
      tree_.recoDimu_eta2.push_back(	     recoDimu_eta2[i]);     
      tree_.recoDimu_phi2.push_back(	     recoDimu_phi2[i]);     
      tree_.recoDimu_charge2.push_back(     recoDimu_charge2[i]);  
      tree_.recoDimu_mass2.push_back(       recoDimu_mass2[i]);           
      /*      tree_.recoDimu_pt.push_back(       recoDimu_pt[i]);           
      tree_.recoDimu_eta.push_back(       recoDimu_eta[i]);           
      tree_.recoDimu_phi.push_back(       recoDimu_phi[i]);           
      tree_.recoDimu_mass.push_back(       recoDimu_mass[i]);           
      tree_.recoDimu_massErr.push_back(       recoDimu_massErr[i]);           
      */       
      }


    tree_.nEleReco=nEleReco;
    tree_.nPFEleReco=nPFEleReco;
    tree_.nLowPtEleReco=nLowPtEleReco;
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
      tree_.recoEle_isPFoverlap.push_back(recoEle_isPFoverlap[i]);
      tree_.recoEle_mvaPFValue.push_back(recoEle_mvaPFValue[i]);
      tree_.recoEle_mvaValue.push_back(recoEle_mvaValue[i]);

    }

    tree_.nDieleReco=nDieleReco;
    for(unsigned int i=0;i<recoDiele_vx.size();i++){
       
      tree_.recoDiele_vx.push_back(           recoDiele_vx[i]);	      
      tree_.recoDiele_vy.push_back(	     recoDiele_vy[i]);	      
      tree_.recoDiele_vz.push_back(	     recoDiele_vz[i]);	      
      tree_.recoDiele_vtxchi2.push_back(     recoDiele_vtxchi2[i]);  
      tree_.recoDiele_vtxndof.push_back(     recoDiele_vtxndof[i]);  
      tree_.recoDiele_vtxprob.push_back(     recoDiele_vtxprob[i]);  
      tree_.recoDiele_pt1.push_back(  	     recoDiele_pt1[i]);      
      tree_.recoDiele_eta1.push_back(	     recoDiele_eta1[i]);     
      tree_.recoDiele_phi1.push_back(	     recoDiele_phi1[i]);     
      tree_.recoDiele_charge1.push_back(     recoDiele_charge1[i]);  
      tree_.recoDiele_mass1.push_back(       recoDiele_mass1[i]);     
      tree_.recoDiele_isPF1.push_back(       recoDiele_isPF1[i]);     
      tree_.recoDiele_isLowPt1.push_back(       recoDiele_isLowPt1[i]);     
      tree_.recoDiele_isPFoverlap1.push_back(       recoDiele_isPFoverlap1[i]);     
      tree_.recoDiele_pt2.push_back(  	     recoDiele_pt2[i]);      
      tree_.recoDiele_eta2.push_back(	     recoDiele_eta2[i]);     
      tree_.recoDiele_phi2.push_back(	     recoDiele_phi2[i]);     
      tree_.recoDiele_charge2.push_back(     recoDiele_charge2[i]);  
      tree_.recoDiele_mass2.push_back(       recoDiele_mass2[i]);           
      tree_.recoDiele_isPF2.push_back(       recoDiele_isPF2[i]);     
      tree_.recoDiele_isLowPt2.push_back(       recoDiele_isLowPt2[i]);     
      tree_.recoDiele_isPFoverlap2.push_back(       recoDiele_isPFoverlap2[i]);     
      /*      tree_.recoDiele_pt.push_back(       recoDiele_pt[i]);           
      tree_.recoDiele_eta.push_back(       recoDiele_eta[i]);           
      tree_.recoDiele_phi.push_back(       recoDiele_phi[i]);           
      tree_.recoDiele_mass.push_back(       recoDiele_mass[i]);           
      tree_.recoDiele_massErr.push_back(       recoDiele_massErr[i]);           
      */       
    }


    tree->Fill();

 }


void TQGenAnalyzer::beginJob()
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

  tree->Branch("nDimuReco", &tree_.nDimuReco, "nDimuReco/I");
  tree->Branch("recoDimu_vx",     &tree_.recoDimu_vx);
  tree->Branch("recoDimu_vy",     &tree_.recoDimu_vy);
  tree->Branch("recoDimu_vz",     &tree_.recoDimu_vz);
  tree->Branch("recoDimu_vtxchi2",     &tree_.recoDimu_vtxchi2);
  tree->Branch("recoDimu_vtxndof",     &tree_.recoDimu_vtxndof);
  tree->Branch("recoDimu_vtxprob",     &tree_.recoDimu_vtxprob);
  tree->Branch("recoDimu_pt1",     &tree_.recoDimu_pt1);
  tree->Branch("recoDimu_eta1",     &tree_.recoDimu_eta1);
  tree->Branch("recoDimu_phi1",     &tree_.recoDimu_phi1);
  tree->Branch("recoDimu_charge1",     &tree_.recoDimu_charge1);
  tree->Branch("recoDimu_mass1",     &tree_.recoDimu_mass1);
  tree->Branch("recoDimu_pt2",     &tree_.recoDimu_pt2);
  tree->Branch("recoDimu_eta2",     &tree_.recoDimu_eta2);
  tree->Branch("recoDimu_phi2",     &tree_.recoDimu_phi2);
  tree->Branch("recoDimu_charge2",     &tree_.recoDimu_charge2);
  tree->Branch("recoDimu_mass2",     &tree_.recoDimu_mass2);
  tree->Branch("recoDimu_pt",     &tree_.recoDimu_pt);
  tree->Branch("recoDimu_eta",     &tree_.recoDimu_eta);
  tree->Branch("recoDimu_phi",     &tree_.recoDimu_phi);
  tree->Branch("recoDimu_mass",     &tree_.recoDimu_mass);
  tree->Branch("recoDimu_massErr",     &tree_.recoDimu_massErr);


  tree->Branch("nEleReco", &tree_.nEleReco, "nEleReco/I");
  tree->Branch("nPFEleReco", &tree_.nPFEleReco, "nPFEleReco/I");
  tree->Branch("nLowPtEleReco", &tree_.nLowPtEleReco, "nLowPtEleReco/I");
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

  tree->Branch("recoEle_mvaPFValue",&tree_.recoEle_mvaPFValue);

  tree->Branch("recoEle_mvaValue",&tree_.recoEle_mvaValue);


  tree->Branch("nDieleReco", &tree_.nDieleReco, "nDieleReco/I");
  tree->Branch("recoDiele_vx",     &tree_.recoDiele_vx);
  tree->Branch("recoDiele_vy",     &tree_.recoDiele_vy);
  tree->Branch("recoDiele_vz",     &tree_.recoDiele_vz);
  tree->Branch("recoDiele_vtxchi2",     &tree_.recoDiele_vtxchi2);
  tree->Branch("recoDiele_vtxndof",     &tree_.recoDiele_vtxndof);
  tree->Branch("recoDiele_vtxprob",     &tree_.recoDiele_vtxprob);
  tree->Branch("recoDiele_pt1",     &tree_.recoDiele_pt1);
  tree->Branch("recoDiele_eta1",     &tree_.recoDiele_eta1);
  tree->Branch("recoDiele_phi1",     &tree_.recoDiele_phi1);
  tree->Branch("recoDiele_charge1",     &tree_.recoDiele_charge1);
  tree->Branch("recoDiele_mass1",     &tree_.recoDiele_mass1);
  tree->Branch("recoDiele_isPF1",     &tree_.recoDiele_isPF1);
  tree->Branch("recoDiele_isLowPt1",     &tree_.recoDiele_isLowPt1);
  tree->Branch("recoDiele_isPFoverlap1",     &tree_.recoDiele_isPFoverlap1);
  tree->Branch("recoDiele_pt2",     &tree_.recoDiele_pt2);
  tree->Branch("recoDiele_eta2",     &tree_.recoDiele_eta2);
  tree->Branch("recoDiele_phi2",     &tree_.recoDiele_phi2);
  tree->Branch("recoDiele_charge2",     &tree_.recoDiele_charge2);
  tree->Branch("recoDiele_mass2",     &tree_.recoDiele_mass2);
  tree->Branch("recoDiele_isPF2",     &tree_.recoDiele_isPF2);
  tree->Branch("recoDiele_isLowPt2",     &tree_.recoDiele_isLowPt2);
  tree->Branch("recoDiele_isPFoverlap2",     &tree_.recoDiele_isPFoverlap2);
  tree->Branch("recoDiele_pt",     &tree_.recoDiele_pt);
  tree->Branch("recoDiele_eta",     &tree_.recoDiele_eta);
  tree->Branch("recoDiele_phi",     &tree_.recoDiele_phi);
  tree->Branch("recoDiele_mass",     &tree_.recoDiele_mass);
  tree->Branch("recoDiele_massErr",     &tree_.recoDiele_massErr);


}

void
TQGenAnalyzer::endJob()
{
}

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

  tree_.recoDimu_vx.clear();
  tree_.recoDimu_vy.clear() ;
  tree_.recoDimu_vz.clear();
  tree_.recoDimu_vtxchi2.clear();
  tree_.recoDimu_vtxndof.clear();
  tree_.recoDimu_vtxprob.clear();
  tree_.recoDimu_pt1.clear();
  tree_.recoDimu_eta1.clear();
  tree_.recoDimu_phi1.clear();
  tree_.recoDimu_charge1.clear();
  tree_.recoDimu_mass1.clear();
  tree_.recoDimu_pt2.clear();
  tree_.recoDimu_eta2.clear();
  tree_.recoDimu_phi2.clear();
  tree_.recoDimu_charge2.clear();
  tree_.recoDimu_mass2.clear();
  tree_.recoDimu_pt.clear();
  tree_.recoDimu_eta.clear();
  tree_.recoDimu_phi.clear();
  tree_.recoDimu_mass.clear();
  tree_.recoDimu_massErr.clear();


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
  tree_.recoEle_isPFoverlap.clear();

  tree_.recoEle_mvaPFValue.clear();
  
  tree_.recoEle_mvaValue.clear();


  tree_.recoDiele_vx.clear();
  tree_.recoDiele_vy.clear() ;
  tree_.recoDiele_vz.clear();
  tree_.recoDiele_vtxchi2.clear();
  tree_.recoDiele_vtxndof.clear();
  tree_.recoDiele_vtxprob.clear();
  tree_.recoDiele_pt1.clear();
  tree_.recoDiele_eta1.clear();
  tree_.recoDiele_phi1.clear();
  tree_.recoDiele_charge1.clear();
  tree_.recoDiele_mass1.clear();
  tree_.recoDiele_isPF1.clear();
  tree_.recoDiele_isLowPt1.clear();
  tree_.recoDiele_isPFoverlap1.clear();
  tree_.recoDiele_pt2.clear();
  tree_.recoDiele_eta2.clear();
  tree_.recoDiele_phi2.clear();
  tree_.recoDiele_charge2.clear();
  tree_.recoDiele_mass2.clear();
  tree_.recoDiele_isPF2.clear();
  tree_.recoDiele_isLowPt2.clear();
  tree_.recoDiele_isPFoverlap2.clear();
  tree_.recoDiele_pt.clear();
  tree_.recoDiele_eta.clear();
  tree_.recoDiele_phi.clear();
  tree_.recoDiele_mass.clear();
  tree_.recoDiele_massErr.clear();

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
  if(year==2016) puWeightFile="pileup_2016.root";
  else   if(year==2017) puWeightFile="pileup_2017.root";
  else   if(year==2018) puWeightFile="pileup_2018.root";
  else   if(year==0) puWeightFile="pileup_ALL.root";

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

DEFINE_FWK_MODULE(TQGenAnalyzer);
