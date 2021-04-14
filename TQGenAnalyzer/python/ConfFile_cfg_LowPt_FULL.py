import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.VarParsing as VarParsing

from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
#from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP

process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            'file:/afs/cern.ch/user/s/soffi/public/10C23D4F-94BD-E811-9588-E0071B7B2320.root'
                                '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/854B1DC0-2F71-694D-A3F5-8DC1CDE1EF18.root'


                )
                            )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_tetraquarks_wLowPt.root")
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] :
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplizer_seq = cms.Sequence()

#loading low pt mva
process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
#process.ntuplizer_seq *= process.lowPtGsfElectronID

#loading pf retrained mva
mvaConfigsForEleProducer = cms.VPSet( )
from python.mvaElectronID_BParkRetrain_cff import mvaEleID_BParkRetrain_producer_config
mvaConfigsForEleProducer.append( mvaEleID_BParkRetrain_producer_config )

# The producer to compute the MVA input variables which are not accessible with the cut parser

from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma8XObjectUpdateModifier,egamma9X105XUpdateModifier,prependEgamma8XObjectUpdateModifier
ele9X105XUpdateModifier=egamma9X105XUpdateModifier.clone(
    phoPhotonIso = "",
    phoNeutralHadIso = "",
    phoChargedHadIso = "",
    phoChargedHadWorstVtxIso = "",
    phoChargedHadWorstVtxConeVetoIso = "",
    phoChargedHadPFPVIso = ""
)
#we have dataformat changes to 106X so to read older releases we use egamma updators
slimmedElectronsTo106X = cms.EDProducer("ModifiedElectronProducer",
    src = cms.InputTag("slimmedElectrons"),
    modifierConfig = cms.PSet( modifications = cms.VPSet(ele9X105XUpdateModifier) )
)
electronMVAVariableHelper = cms.EDProducer('GsfElectronMVAVariableHelper',
  # The module automatically detects AOD vs miniAOD, so we configure both
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),
  vertexCollection = cms.InputTag("offlinePrimaryVertices"),
  beamSpot         = cms.InputTag("offlineBeamSpot"),
conversions = cms.InputTag("allConversions"),
  # miniAOD case
#  srcMiniAOD              = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
  srcMiniAOD              = cms.InputTag('slimmedElectronsTo106X'),
  vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
  conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

electronMVAValueMapProducer = cms.EDProducer(
  'ElectronMVAValueMapProducer',
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),  
  # miniAOD case
#  srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
  srcMiniAOD = cms.InputTag('slimmedElectrons'),  
  #srcMiniAOD = cms.InputTag('regressionForEle:regressedElectrons'),
    
  # MVA configurations
  mvaConfigurations = mvaConfigsForEleProducer
)

egmGsfElectronIDs = cms.EDProducer(
    "VersionedGsfElectronIdProducer",
    physicsObjectSrc = cms.InputTag('gedGsfElectrons'),
    physicsObjectIDs = cms.VPSet( )
)

egmGsfElectronIDTask = cms.Task(
    electronMVAVariableHelper,
    electronMVAValueMapProducer,
    egmGsfElectronIDs,
)

egmGsfElectronIDSequence = cms.Sequence(egmGsfElectronIDTask)





# this is for the LowPt energy regression
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To50_mean_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To50_mean_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To50_sigma_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To50_sigma_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To50_mean_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To50_mean_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To50_sigma_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To50_sigma_2017UL"),
         connect = cms.string("sqlite_file:lowPtEleReg_2017UL_25112020.db")))

#setup input variables
options = VarParsing.VarParsing('analysis')
options.register ('isMC',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,           # string, int, float, bool
                  "Bool isMC")
options.register ('isSignal',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,           # string, int, float, bool
                  "Bool isSignal")
options.parseArguments()

if options.isMC and options.isSignal :
    print'Sample is MC Signal'
    xsec=1.

if options.isMC and not options.isSignal : print'Sample is MC Background'
if options.isMC == 0 : print'Sample is Data'

if options.isMC and options.isSignal : index = 0
if options.isMC and not options.isSignal : index = -100
if options.isMC ==0 : index==100


#setup json file
if (options.isMC==False):
    print "applying json"
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    #change json below once you run on new data
    JSONfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
    process.source.lumisToProcess.extend(myLumis)
#    print myLumis



process.GenAnalysis = cms.EDAnalyzer('TQGenAnalyzer',
                                     generatorInfo= cms.InputTag("generator"),
                                     prunedGenParticles    = cms.InputTag("prunedGenParticles"),
                                     patMuons = cms.InputTag("slimmedMuons"),
                                     patElectrons = cms.InputTag("slimmedElectrons"),                 # MINIAOD
                                     patElectronsLowPt = cms.InputTag("slimmedLowPtElectrons"), #change when running on our old signals
                                     vtx=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     rho= cms.InputTag('fixedGridRhoAll'),
                                     PileUp = cms.InputTag('slimmedAddPileupInfo'),
                                     bits         = cms.InputTag("TriggerResults::HLT"),
                                     flags        = cms.InputTag("TriggerResults::SIM"),
                                     sampleIndex  = cms.untracked.int32(index),
                                     sampleXsec  = cms.untracked.double(xsec),
                                     drForCleaning = cms.double(0.03),
                                     dzForCleaning = cms.double(0.5), ##keep tighter dZ to check overlap of pfEle with lowPt (?)
                                     mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
                                     mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'), 
                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
#                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),
                                     mvaIdPF = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'), 
#                                     mvaValue = cms.InputTag('lowPtGsfElectronID'),
                                     mvaValue = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
                                     mvaId = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90'), 

                                     gsfRegressionConfig = cms.PSet(
                                         modifierName = cms.string('EGRegressionModifierV3'),
                                         rhoTag = cms.string('fixedGridRhoFastjetAll'),
                                         useClosestToCentreSeedCrysDef = cms.bool(False),
                                         maxRawEnergyForLowPtEBSigma = cms.double(-1),
                                         maxRawEnergyForLowPtEESigma = cms.double(1200.),
                                         eleRegs = cms.PSet(
                                             ecalOnlyMean = cms.PSet(
                                                 rangeMinLowEt = cms.double(0.2),
                                                 rangeMaxLowEt = cms.double(2.0),
                                                 rangeMinHighEt = cms.double(-1.),
                                                 rangeMaxHighEt = cms.double(3.0),
                                                 forceHighEnergyTrainingIfSaturated = cms.bool(True),
                                                 lowEtHighEtBoundary = cms.double(999999.),
                                                 ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_mean"),
                                                 ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_mean"),
                                                 eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_mean"),
                                                 eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_mean"),
                                             ),
                                                 ecalOnlySigma = cms.PSet(
                                                     rangeMinLowEt = cms.double(0.0002),
                                                     rangeMaxLowEt = cms.double(0.5),
                                                     rangeMinHighEt = cms.double(0.0002),
                                                     rangeMaxHighEt = cms.double(0.5),
                                                     forceHighEnergyTrainingIfSaturated = cms.bool(True),
                                                     lowEtHighEtBoundary = cms.double(999999.),
                                                     ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_sigma"),
                                                     ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To50_sigma"),
                                                     eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_sigma"),
                                                     eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To50_sigma"),
                                                 ),
                                             epComb = cms.PSet(
                                                 ecalTrkRegressionConfig = cms.PSet(
                                                     rangeMinLowEt = cms.double(0.2),
                                                     rangeMaxLowEt = cms.double(2.0),
                                                     rangeMinHighEt = cms.double(0.2),
                                                     rangeMaxHighEt = cms.double(2.0),
                                                     lowEtHighEtBoundary = cms.double(999999.),
                                                     forceHighEnergyTrainingIfSaturated = cms.bool(False),
                                                     ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_mean'),
                                                     ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_mean'),
                                                     eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_mean'),
                                                     eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_mean'),
                                                ),
                                                 ecalTrkRegressionUncertConfig = cms.PSet(
                                                     rangeMinLowEt = cms.double(0.0002),
                                                     rangeMaxLowEt = cms.double(0.5),
                                                     rangeMinHighEt = cms.double(0.0002),
                                                     rangeMaxHighEt = cms.double(0.5),
                                                     lowEtHighEtBoundary = cms.double(999999.),
                                                     forceHighEnergyTrainingIfSaturated = cms.bool(False),
                                                     ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_sigma'),
                                                     ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To50_sigma'),
                                                     eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_sigma'),
                                                     eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To50_sigma'),
                                                 ),
                                                 maxEcalEnergyForComb=cms.double(200.),
                                                 minEOverPForComb=cms.double(0.025),
                                                 maxEPDiffInSigmaForComb=cms.double(15.),
                                                 maxRelTrkMomErrForComb=cms.double(10.),
                                             )
                                         )

                                     )
                                )

#process.ntuplizer_seq *= process.lowPtGsfElectronID
process.ntuplizer_seq *= process.GenAnalysis

process.ntuplizer_path = cms.Path(#process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.schedule = cms.Schedule(process.ntuplizer_path)
