import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

#1. setting GT
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            'file:/afs/cern.ch/user/s/soffi/public/10C23D4F-94BD-E811-9588-E0071B7B2320.root'
#                                '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/854B1DC0-2F71-694D-A3F5-8DC1CDE1EF18.root'
#                                    'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/GenProduction/BPH-RunIISummer20UL18MiniAODv2-00008.root'
#                                    '/store/cmst3/group/bpark/BToKmumu_1000Events_MINIAOD.root'
                                    'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/GenProduction/root_files_23042021/BPH-RunIISummer20UL18MiniAOD_2.root'
                )
                            )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_tetraquarks_wLowPt.root")
)

#2. setting input variables
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



#3. setting up post reco tools from EGM to have all IDs and SS
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2018-Prompt')  





process.ntuplizer_seq = cms.Sequence()

#4. setting json file
if (options.isMC==False):
    print "applying json"
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    #change json below once you run on new data
    JSONfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
    process.source.lumisToProcess.extend(myLumis)
#    print myLumis



#5. setting Regression configs
from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL
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





#6. setting Analyzer
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
 
                                     mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values'),
#                                     mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),
                                     mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90'), 
                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values'),
#                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),
                                     mvaIdPF = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90'), 
                                     mvaValue = cms.InputTag('lowPtGsfElectronID'),
#                                     mvaValue = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),
#                                     mvaValue = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values'),
                                     mvaId = cms.InputTag('egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90'), 

                                     gsfRegressionConfig = cms.PSet(
                                         modifierName = cms.string('EGRegressionModifierV3'),
                                         rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
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


#8. setting full path

process.p = cms.Path((process.egammaPostRecoSeq*process.egmGsfElectronIDSequence) * process.GenAnalysis)
#process.p = cms.Path(process.egammaPostRecoSeq+process.lowPtGsfElectronID+ process.GenAnalysis)
#process.p = cms.Path((process.egmGsfElectronIDSequence) * process.GenAnalysis)


#process.ntuplizer_seq *= process.egammaPostRecoSeq

#process.ntuplizer_seq *=process.GenAnalysis
#process.ntuplizer_path = cms.Path(
#                                  process.ntuplizer_seq)
#process.schedule = cms.Schedule(process.ntuplizer_path)
