// -*- C++ -*-
//
// Package:    RecoEgamma/ElectronIdentification
// Class:      ElectronMVANtuplizer
//
/**\class ElectronMVANtuplizer ElectronMVANtuplizer.cc RecoEgamma/ElectronIdentification/plugins/ElectronMVANtuplizer.cc

 Description: Ntuplizer for training and testing electron MVA IDs.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jonas REMBSER
//         Created:  Thu, 22 Mar 2018 14:54:24 GMT
//
//


// user include files
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "RecoEgamma/EgammaTools/interface/MVAVariableManager.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <TTree.h>
#include <TFile.h>
#include <Math/VectorUtil.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//

class ElectronMVANtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronMVANtuplizer(const edm::ParameterSet&);
      ~ElectronMVANtuplizer() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      void beginJob() override;
      void analyze(const edm::Event&, const edm::EventSetup&) override;
      void endJob() override;

      void findFirstNonElectronMother2(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);

      template<class T, class V>
      int matchToTruth(const T &el, const V &genParticles, int &genIdx);

      // ----------member data ---------------------------

      // for AOD case
      const edm::EDGetToken src_;
      const edm::EDGetToken vertices_;
      const edm::EDGetToken pileup_;
      const edm::EDGetToken genParticles_;

      // for miniAOD case
      const edm::EDGetToken  srcMiniAOD_;
      const edm::EDGetToken verticesMiniAOD_;
      const edm::EDGetToken pileupMiniAOD_;
      const edm::EDGetToken genParticlesMiniAOD_;

      // other
      TTree* tree_;

      std::unique_ptr<ModifyObjectValueBase> regressionGsf_;  // Gsf
      MVAVariableManager<reco::GsfElectron> mvaVarMngr_;
      std::vector<float> vars_;
      int nVars_;

      //global variables
      int nEvent_, nRun_, nLumi_;
      int genNpu_;
      int vtxN_;

      // electron variables
      float eleQ_;
      int ele3Q_;
      int matchedToGenEle_;
      int matchedGenIdx_;
      double trackPAtVtx_;
      double eSuperClusterOverP_;
      double IoEmIop_;
      double CorreSuperClusterOverP_;
      double CorrIoEmIop_;
      double elePt_;

      // gap variables
      bool eleIsEB_;
      bool eleIsEE_;
      bool eleIsEBEtaGap_;
      bool eleIsEBPhiGap_;
      bool eleIsEBEEGap_;
      bool eleIsEEDeeGap_;
      bool eleIsEERingGap_;

      // to hold ID decisions and categories
      std::vector<int> mvaPasses_;
      std::vector<float> mvaValues_;
      std::vector<int> mvaCats_;

      // config
      const bool isMC_;
      const double deltaR_;
      const double ptThreshold_;

      // ID decisions objects
      const std::vector< std::string > eleMapTags_;
      std::vector< edm::EDGetTokenT< edm::ValueMap<bool> > > eleMapTokens_;
      const std::vector< std::string > eleMapBranchNames_;
      const size_t nEleMaps_;

      // MVA values and categories (optional)
      const std::vector< std::string > valMapTags_;
      std::vector< edm::EDGetTokenT<edm::ValueMap<float> > > valMapTokens_;
      const std::vector< std::string > valMapBranchNames_;
      const size_t nValMaps_;

      const std::vector< std::string > mvaCatTags_;
      std::vector< edm::EDGetTokenT<edm::ValueMap<int> > > mvaCatTokens_;
      const std::vector< std::string > mvaCatBranchNames_;
      const size_t nCats_;
};

//
// constants, enums and typedefs
//

enum ElectronMatchType {
                        UNMATCHED,
                        TRUE_PROMPT_ELECTRON,
                        TRUE_ELECTRON_FROM_TAU,
                        TRUE_NON_PROMPT_ELECTRON,
                       }; // The last does not include tau parents

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronMVANtuplizer::ElectronMVANtuplizer(const edm::ParameterSet& iConfig)
 :
  src_                   (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("src") )),
  vertices_              (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
  pileup_                (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileup"))),
  genParticles_          (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  srcMiniAOD_            (consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("srcMiniAOD"))),
//  srcMiniAOD_		 (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("srcMiniAOD") )),
  verticesMiniAOD_       (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesMiniAOD"))),
  pileupMiniAOD_         (consumes<std::vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("pileupMiniAOD"))),
  genParticlesMiniAOD_   (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticlesMiniAOD"))),
  mvaVarMngr_            (iConfig.getParameter<std::string>("variableDefinition")),
  isMC_                  (iConfig.getParameter<bool>("isMC")),
  deltaR_                (iConfig.getParameter<double>("deltaR")),
  ptThreshold_           (iConfig.getParameter<double>("ptThreshold")),
  eleMapTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAs")),
  eleMapBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVALabels")),
  nEleMaps_              (eleMapBranchNames_.size()),
  valMapTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAValMaps")),
  valMapBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAValMapLabels")),
  nValMaps_              (valMapBranchNames_.size()),
  mvaCatTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVACats")),
  mvaCatBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVACatLabels")),
  nCats_                 (mvaCatBranchNames_.size())
{
    // Regression stuff - GSF electrons
    if( iConfig.existsAs<edm::ParameterSet>("gsfRegressionConfig") ) {
      const edm::ParameterSet& iconf = iConfig.getParameterSet("gsfRegressionConfig");
      const std::string& mname = iconf.getParameter<std::string>("modifierName");
      ModifyObjectValueBase* plugin =
        ModifyObjectValueFactory::get()->create(mname,iconf);
      regressionGsf_.reset(plugin);
      edm::ConsumesCollector sumes = consumesCollector();
      regressionGsf_->setConsumes(sumes);
    } else {
      regressionGsf_.reset(nullptr);
    }
    // eleMaps
    for (size_t k = 0; k < nEleMaps_; ++k) {

        eleMapTokens_.push_back(consumes<edm::ValueMap<bool> >(edm::InputTag(eleMapTags_[k])));

        // Initialize vectors for holding ID decisions
        mvaPasses_.push_back(0);
    }

    // valMaps
    for (size_t k = 0; k < nValMaps_; ++k) {
        valMapTokens_.push_back(consumes<edm::ValueMap<float> >(edm::InputTag(valMapTags_[k])));

        // Initialize vectors for holding MVA values
        mvaValues_.push_back(0.0);
    }

    // categories
    for (size_t k = 0; k < nCats_; ++k) {
        mvaCatTokens_.push_back(consumes<edm::ValueMap<int> >(edm::InputTag(mvaCatTags_[k])));

        // Initialize vectors for holding MVA values
        mvaCats_.push_back(0);
    }

   // Book tree
   usesResource(TFileService::kSharedResource);
   edm::Service<TFileService> fs ;
   tree_  = fs->make<TTree>("tree","tree");

    // Regression stuff - lowPtElectrons
 //  if( cfg.existsAs<edm::ParameterSet>("lowPtRegressionConfig") ) {
 //    const edm::ParameterSet& iconf = cfg.getParameterSet("lowPtRegressionConfig");
 //    const std::string& mname = iconf.getParameter<std::string>("modifierName");
 //    ModifyObjectValueBase* plugin =
 //      ModifyObjectValueFactory::get()->create(mname,iconf);
 //    regression_.reset(plugin);
 //    edm::ConsumesCollector sumes = consumesCollector();
 //    regression_->setConsumes(sumes);
 //  } else {
 //    regression_.reset(nullptr);
 //  }

   nVars_ = mvaVarMngr_.getNVars();

   tree_->Branch("nEvent",  &nEvent_);
   tree_->Branch("nRun",    &nRun_);
   tree_->Branch("nLumi",   &nLumi_);
   if (isMC_) tree_->Branch("genNpu", &genNpu_);
   tree_->Branch("vtxN",   &vtxN_);

   tree_->Branch("ele_q",&eleQ_);
   tree_->Branch("ele_3q",&ele3Q_);

   if (isMC_) {
       tree_->Branch("matchedToGenEle",   &matchedToGenEle_);
   }

   // Has to be in two different loops
   for (int i = 0; i < nVars_; ++i) {
       vars_.push_back(0.0);
   }
   for (int i = 0; i < nVars_; ++i) {
       tree_->Branch(mvaVarMngr_.getName(i).c_str(), &vars_[i]);
   }

   tree_->Branch("ele_isEB",&eleIsEB_);
   tree_->Branch("ele_isEE",&eleIsEE_);
   tree_->Branch("ele_isEBEtaGap",&eleIsEBEtaGap_);
   tree_->Branch("ele_isEBPhiGap",&eleIsEBPhiGap_);
   tree_->Branch("ele_isEBEEGap", &eleIsEBEEGap_);
   tree_->Branch("ele_isEEDeeGap",&eleIsEEDeeGap_);
   tree_->Branch("ele_isEERingGap",&eleIsEERingGap_);
   tree_->Branch("ele_trackPAtVtx", &trackPAtVtx_);
   tree_->Branch("ele_RegSuperClusterOverP",&CorreSuperClusterOverP_);
   tree_->Branch("ele_RegIoEmIop",&CorrIoEmIop_);
   tree_->Branch("ele_eSuperClusterOverP",&eSuperClusterOverP_);
   tree_->Branch("ele_IoEmIopR",&IoEmIop_);
   tree_->Branch("ele_Regpt",&elePt_);

   // IDs
   for (size_t k = 0; k < nValMaps_; ++k) {
       tree_->Branch(valMapBranchNames_[k].c_str() ,  &mvaValues_[k]);
   }

   for (size_t k = 0; k < nEleMaps_; ++k) {
       tree_->Branch(eleMapBranchNames_[k].c_str() ,  &mvaPasses_[k]);
   }

   for (size_t k = 0; k < nCats_; ++k) {
       tree_->Branch(mvaCatBranchNames_[k].c_str() ,  &mvaCats_[k]);
   }

   // All tokens for event content needed by this MVA
   // Tags from the variable helper
   for (auto &tag : mvaVarMngr_.getHelperInputTags()) {
       consumes<edm::ValueMap<float>>(tag);
   }
   for (auto &tag : mvaVarMngr_.getGlobalInputTags()) {
       consumes<double>(tag);
   }
}


ElectronMVANtuplizer::~ElectronMVANtuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronMVANtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool useEnergyRegression=1;
  bool debug =1;
    // Fill global event info
    nEvent_ = iEvent.id().event();
    nRun_   = iEvent.id().run();
    nLumi_  = iEvent.luminosityBlock();

  if (useEnergyRegression) {
  //  regression_->setEvent(event);
   // regression_->setEventContent(setup);
    regressionGsf_->setEvent(iEvent);
    regressionGsf_->setEventContent(iSetup);
  }

    // Retrieve Vertecies
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vertices_, vertices);
    if( !vertices.isValid() ){
      iEvent.getByToken(verticesMiniAOD_,vertices);
      if( !vertices.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD vertex collection " << std::endl;
    }

    vtxN_ = vertices->size();

    // Retrieve Pileup Info
    if(isMC_) {
        edm::Handle<std::vector< PileupSummaryInfo > >  pileup;
        iEvent.getByToken(pileup_, pileup);
        if( !pileup.isValid() ){
          iEvent.getByToken(pileupMiniAOD_,pileup);
          if( !pileup.isValid() )
            throw cms::Exception(" Collection not found: ")
              << " failed to find a standard AOD or miniAOD pileup collection " << std::endl;
        }

        // Fill with true number of pileup
       for(const auto& pu : *pileup)
       {
           int bx = pu.getBunchCrossing();
           if(bx == 0)
           {
               genNpu_ = pu.getPU_NumInteractions();
               break;
           }
       }
    }

    // Retrieve genParticles
    edm::Handle<edm::View<reco::GenParticle> >  genParticles;
    if(isMC_) {
        iEvent.getByToken(genParticles_, genParticles);
        if( !genParticles.isValid() ){
          iEvent.getByToken(genParticlesMiniAOD_, genParticles);
          if( !genParticles.isValid() )
            throw cms::Exception(" Collection not found: ")
              << " failed to find a standard AOD or miniAOD genParticle collection " << std::endl;
        }
    }


    edm::Handle<edm::View<reco::GsfElectron> > src;

    // Retrieve the collection of particles from the event.
    // If we fail to retrieve the collection with the standard AOD
    // name, we next look for the one with the stndard miniAOD name.
    iEvent.getByToken(src_, src);
    if( !src.isValid() ){
      iEvent.getByToken(srcMiniAOD_,src);
      if( !src.isValid() )
        throw cms::Exception(" Collection not found: ")
          << " failed to find a standard AOD or miniAOD particle collection " << std::endl;
    }

        std::cout <<" Look up and save the ID decisions_0" << std::endl;
    // Get MVA decisions
    edm::Handle<edm::ValueMap<bool> > decisions[nEleMaps_];
    for (size_t k = 0; k < nEleMaps_; ++k) {
        std::cout <<" elemaps_" <<k  << std::endl;
        iEvent.getByToken(eleMapTokens_[k],decisions[k]);
    }

    // Get MVA values
    edm::Handle<edm::ValueMap<float> > values[nValMaps_];
    for (size_t k = 0; k < nValMaps_; ++k) {
        std::cout <<" mvaVals_" <<k  << std::endl;
        iEvent.getByToken(valMapTokens_[k],values[k]);
    }

    // Get MVA categories
    edm::Handle<edm::ValueMap<int> > mvaCats[nCats_];
    for (size_t k = 0; k < nCats_; ++k) {
        std::cout <<" mvaCats_" <<k  << std::endl;
        iEvent.getByToken(mvaCatTokens_[k],mvaCats[k]);
    }
    int iEle=0;
    int nEle = src->size();
//	for(auto ele : *src) {
   for(int iEle = 0; iEle < nEle; ++iEle) {
      // const auto electron  = ele.Ref();
      	auto ele = src->at(iEle);
        eleQ_ = ele.charge();
        ele3Q_ = ele.chargeInfo().isGsfCtfScPixConsistent;

       // std::cout <<" after charges "  << iEle<< std::endl;
        if (ele.pt() < ptThreshold_) {
            continue;
        }


	trackPAtVtx_ = ele.superCluster()->energy()/ele.eSuperClusterOverP();
	eSuperClusterOverP_ = ele.eSuperClusterOverP();
	IoEmIop_= 1.0/ele.superCluster()->energy()-1.0/trackPAtVtx_;
       if(debug) std::cout <<" before mva manager "  << ele.superCluster()->rawEnergy() << " "<<ele.ecalEnergy() << " " << ele.trackMomentumAtVtx().R() <<" "<<  ele.eSuperClusterOverP()<< " " <<  IoEmIop_ << std::endl;
    	if (useEnergyRegression) {
//  	  reco::GsfElectron newElectron(ele);
   //	  pre_ecal = newElectron.correctedEcalEnergy();
    	 // pre_ecaltrk = newElectron.energy();
    	  regressionGsf_->modifyObject(ele);
    //	  post_ecal = newElectron.correctedEcalEnergy();
  //	  post_ecaltrk = newElectron.energy();
           }
	CorreSuperClusterOverP_= eSuperClusterOverP_*ele.ecalEnergy()/ele.superCluster()->energy();
	CorrIoEmIop_= 1.0/ele.ecalEnergy()-1.0/trackPAtVtx_;
	elePt_ = ele.pt();
   	//edm::Ref<pat::ElectronCollection> ref(src,iEle);
//	auto ePtr = ele.sourceCandidatePtr(iEle)->gsfElectronRef();
//	edm::Ptr<pat::Electron> ePtr = ele;
         auto electron =  src->ptrAt(iEle);
	std::cout <<"   electron pointer " << electron<< std::endl;
       // std::cout <<" before charges "  << iEle<< std::endl;

  if(debug)   { 
  std::cout <<" before mva manager "  << ele.superCluster()->rawEnergy() << " "<<ele.ecalEnergy() << " " << ele.trackMomentumAtVtx().R() <<" "<<  ele.eSuperClusterOverP()<< " " <<  IoEmIop_ << std::endl;
  std::cout <<" before mva manager check "  << ele.superCluster()->rawEnergy() << " "<<ele.ecalEnergy() << " " << trackPAtVtx_ <<" "<<  CorreSuperClusterOverP_<< " " <<  CorrIoEmIop_ << std::endl;
  std::cout <<" before mva manager check "  << ele.superCluster()->rawEnergy() << " "<<ele.pt() << std::endl;
}
        for (int iVar = 0; iVar < nVars_; ++iVar) {
            vars_[iVar] = mvaVarMngr_.getValue(iVar,electron, iEvent);
        }

        std::cout <<" after mva manager"  << iEle<< std::endl;
        if (isMC_) {
            matchedToGenEle_ = matchToTruth( &ele, genParticles, matchedGenIdx_);
        }

        // gap variables
        eleIsEB_ = ele.isEB();
        eleIsEE_ = ele.isEE();
        eleIsEBEEGap_ = ele.isEBEEGap();
        eleIsEBEtaGap_ = ele.isEBEtaGap();
        eleIsEBPhiGap_ = ele.isEBPhiGap();
        eleIsEEDeeGap_ = ele.isEEDeeGap();
        eleIsEERingGap_ = ele.isEERingGap();

        //
        std::cout <<" Look up and save the ID decisions" << std::endl;
        //
        for (size_t k = 0; k < nEleMaps_; ++k) {
          mvaPasses_[k] = (int)(*decisions[k])[electron];
        }

        for (size_t k = 0; k < nValMaps_; ++k) {
          mvaValues_[k] = (*values[k])[electron];
        }

        for (size_t k = 0; k < nCats_; ++k) {
          mvaCats_[k] = (*mvaCats[k])[electron];
        }


        tree_->Fill();
   	iEle++;
    }

}

void ElectronMVANtuplizer::findFirstNonElectronMother2(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == nullptr ){
    edm::LogError  ("ElectronNtuplizer") << "ElectronNtuplizer: ERROR! null candidate pointer, this should never happen";
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother2(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

template<class T, class V>
int ElectronMVANtuplizer::matchToTruth(const T &el, const V &prunedGenParticles, int &genIdx){

  //
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  double deltaR = deltaR_;
  const reco::Candidate *closestElectron = nullptr;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
      genIdx = i;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
 // std::cout << " " << dR << " "<< deltaR  <<" " << genIdx <<std::endl;
  if( closestElectron == nullptr || dR > deltaR  ) {
    return UNMATCHED;
  }

  
  int ancestorPID = -999;
  int ancestorStatus = -999;
 // findFirstNonElectronMother2(closestElectron, ancestorPID, ancestorStatus);

 // std::cout <<closestElectron->pdgId()<< " " << closestElectron->mother()->pdgId()<< " " << dR << " "<< deltaR_ <<std::endl;

/*
  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    edm::LogError  ("ElectronNtuplizer") << "ElectronNtuplizer: ERROR! null candidate pointer, this should never happen";
    return UNMATCHED;
  }

  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons*/
  //

//  std::cout <<closestElectron->pdgId()<< " " << closestElectron->mother()->pdgId() << "" << fabs(ancestorPID) << std::endl;
 // if ( fabs(ancestorPID) == 9900015){ 
  if ( fabs(closestElectron->pdgId()) == 11){ 
 // std::cout <<closestElectron->pdgId()<< " " << closestElectron->mother()->pdgId()<< " " << dR << " "<< genIdx  <<std::endl;
  return TRUE_PROMPT_ELECTRON;
  } 
 else return UNMATCHED;
}

// ------------ method called once each job just before starting event loop  ------------
void
ElectronMVANtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronMVANtuplizer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
ElectronMVANtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src");
    desc.add<edm::InputTag>("vertices");
    desc.add<edm::InputTag>("pileup");
    desc.add<edm::InputTag>("genParticles");
    desc.add<edm::InputTag>("srcMiniAOD");
    desc.add<edm::InputTag>("verticesMiniAOD");
    desc.add<edm::InputTag>("pileupMiniAOD");
    desc.add<edm::InputTag>("genParticlesMiniAOD");
    desc.add<std::string>("variableDefinition");
    desc.add<bool>("isMC");
    desc.add<double>("deltaR", 0.1);
    desc.add<double>("ptThreshold", 2.0);
    desc.addUntracked<std::vector<std::string>>("eleMVAs");
    desc.addUntracked<std::vector<std::string>>("eleMVALabels");
    desc.addUntracked<std::vector<std::string>>("eleMVAValMaps");
    desc.addUntracked<std::vector<std::string>>("eleMVAValMapLabels");
    desc.addUntracked<std::vector<std::string>>("eleMVACats");
    desc.addUntracked<std::vector<std::string>>("eleMVACatLabels");
  //  desc.add<edm::ParameterSet>("gsfRegressionConfig");
    desc.setAllowAnything(); 
    descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronMVANtuplizer);
