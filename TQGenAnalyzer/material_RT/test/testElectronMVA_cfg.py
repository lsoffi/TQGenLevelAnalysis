import FWCore.ParameterSet.Config as cms
# from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import electronMVAVariableHelper
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("ElectronMVANtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# File with the ID variables to include in the Ntuplizer
mvaVariablesFile = "RecoEgamma/ElectronIdentification/data/ElectronIDVariables.txt"

outputFile = "electron_validation_ntuple.root"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it///store/user/ratramon/HNLGen_ntuples/private_BtoDMuN_NtoE_step4/privateHNLstep4_2021Jan22/210122_100537/0000/BPH-mini_10.root'
	'/store/mc/RunIIAutumn18MiniAOD/DoubleElectron_FlatPt-HalfTo50/MINIAODSIM/FlatPU0to70IdealECALforBParking_102X_upgrade2018_realistic_v15-v2/250000/028D33B4-276E-6342-B4C2-2A0EC03F6283.root'
	#'/store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/20000/00D13F2E-6F44-E811-923E-001E0BED0560.root'
   #      '/store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/00000/0293A280-B5F3-E711-8303-3417EBE33927.root'
    )
)

useAOD = False

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
#       'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
        'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_BParkRetrain_cff',
                 ]

process.mvaConfigsForEleProducer = cms.VPSet( )
# Import and add all desired MVAs
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff \
    import mvaEleID_Fall17_noIso_V2_producer_config
process.mvaConfigsForEleProducer.append( mvaEleID_Fall17_noIso_V2_producer_config )
  
  
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff \
    import mvaEleID_Fall17_iso_V2_producer_config
process.mvaConfigsForEleProducer.append( mvaEleID_Fall17_iso_V2_producer_config )

from PhysicsTools.BParkingNano.mvaElectronID_BParkRetrain_cff \
    import mvaEleID_BParkRetrain_producer_config
process.mvaConfigsForEleProducer.append( mvaEleID_BParkRetrain_producer_config )
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#add them to the VID producer
#for idmod in my_id_modules:
 #   setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
    
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP



from Configuration.AlCa.GlobalTag import GlobalTag
#@process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')
# this is for the LowPt energy regression
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_05To20_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_05To20_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_mean"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_eb_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_eb_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("lowPtElectron_ee_ecalTrk_20To50_sigma"),
         tag = cms.string("lowPtElectron_ee_ecalTrk_20To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))




#from PhysicsTools.BParkingNano.electronsBPark_cff import regressionForEle
process.electronMVAVariableHelper = cms.EDProducer('GsfElectronMVAVariableHelper',
  # The module automatically detects AOD vs miniAOD, so we configure both
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),
  vertexCollection = cms.InputTag("offlinePrimaryVertices"),
  beamSpot         = cms.InputTag("offlineBeamSpot"),
conversions = cms.InputTag("allConversions"),
  # miniAOD case
  srcMiniAOD              = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
#  srcMiniAOD              = cms.InputTag('regressionForEle:regressedElectrons'),
  vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpotMiniAOD         = cms.InputTag("offlineBeamSpot"),
  conversionsMiniAOD      = cms.InputTag("reducedEgamma:reducedConversions"),
)

process.electronMVAValueMapProducer = cms.EDProducer(
  'ElectronMVAValueMapProducer',
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),  
  # miniAOD case
  srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),  
  #srcMiniAOD = cms.InputTag('regressionForEle:regressedElectrons'),
    
  # MVA configurations
  mvaConfigurations =process.mvaConfigsForEleProducer
)

process.egmGsfElectronIDs = cms.EDProducer(
    "VersionedGsfElectronIdProducer",
    physicsObjectSrc = cms.InputTag('gedGsfElectrons'),
    physicsObjectIDs = cms.VPSet( )
)

process.egmGsfElectronIDTask = cms.Task(
   process.electronMVAVariableHelper,
   process.electronMVAValueMapProducer,
   process.egmGsfElectronIDs,
)
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDTask)


#process.regressionForEle = cms.EDProducer(
#  'ElectronRegresser1',
#  lowptSrc = cms.InputTag('slimmedLowPtElectrons'),
#  pfSrc    = cms.InputTag('slimmedElectrons'),
#    lowPtRegressionConfig = cms.PSet(
#      modifierName = cms.string('EGRegressionModifierLPV1'),
#      rhoTag = cms.string('fixedGridRhoFastjetAll'),
#      useClosestToCentreSeedCrysDef = cms.bool(False),
#      maxRawEnergyForLowPtEBSigma = cms.double(-1),
#      maxRawEnergyForLowPtEESigma = cms.double(1200.),
#      eleRegs = cms.PSet(
#        ecalOnlyMean = cms.PSet(
#            rangeMinLowEt = cms.double(0.2),
#            rangeMaxLowEt = cms.double(2.0),
#            rangeMinHighEt = cms.double(-1.),
#            rangeMaxHighEt = cms.double(3.0),
#            forceHighEnergyTrainingIfSaturated = cms.bool(True),
#            lowEtHighEtBoundary = cms.double(20.),
#            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_mean"),
#            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_mean"),
#            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_mean"),
#            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_mean"),
#            ),
#        ecalOnlySigma = cms.PSet(
#            rangeMinLowEt = cms.double(0.0002),
#            rangeMaxLowEt = cms.double(0.5),
#            rangeMinHighEt = cms.double(0.0002),
#            rangeMaxHighEt = cms.double(0.5),
#            forceHighEnergyTrainingIfSaturated = cms.bool(True),
#            lowEtHighEtBoundary = cms.double(20.),
#            ebLowEtForestName = cms.string("lowPtElectron_eb_ecalOnly_05To20_sigma"),
#            ebHighEtForestName = cms.string("lowPtElectron_eb_ecalOnly_20To50_sigma"),
#            eeLowEtForestName = cms.string("lowPtElectron_ee_ecalOnly_05To20_sigma"),
#            eeHighEtForestName = cms.string("lowPtElectron_ee_ecalOnly_20To50_sigma"),
#            ),
#        epComb = cms.PSet(
#            ecalTrkRegressionConfig = cms.PSet(
#                rangeMinLowEt = cms.double(0.2),
#                rangeMaxLowEt = cms.double(2.0),
#                rangeMinHighEt = cms.double(0.2),
#                rangeMaxHighEt = cms.double(2.0),
#                lowEtHighEtBoundary = cms.double(20.),
#                forceHighEnergyTrainingIfSaturated = cms.bool(False),
#                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_mean'),
#                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_mean'),
#                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_mean'),
#                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_mean'),
#                ),
#            ecalTrkRegressionUncertConfig = cms.PSet(
#                rangeMinLowEt = cms.double(0.0002),
#                rangeMaxLowEt = cms.double(0.5),
#                rangeMinHighEt = cms.double(0.0002),
#                rangeMaxHighEt = cms.double(0.5),
#                lowEtHighEtBoundary = cms.double(20.),
#                forceHighEnergyTrainingIfSaturated = cms.bool(False),
#                ebLowEtForestName = cms.string('lowPtElectron_eb_ecalTrk_05To20_sigma'),
#                ebHighEtForestName = cms.string('lowPtElectron_eb_ecalTrk_20To50_sigma'),
#                eeLowEtForestName = cms.string('lowPtElectron_ee_ecalTrk_05To20_sigma'),
#                eeHighEtForestName = cms.string('lowPtElectron_ee_ecalTrk_20To50_sigma'),
#                ),
#            maxEcalEnergyForComb=cms.double(200.),
#            minEOverPForComb=cms.double(0.025),
#            maxEPDiffInSigmaForComb=cms.double(15.),
#            maxRelTrkMomErrForComb=cms.double(10.),
#            )
#        ),
#      phoRegs = regressionModifier106XUL.phoRegs.clone()
#    ),
#    gsfRegressionConfig = cms.PSet(
#      modifierName = cms.string('EGRegressionModifierV3'),
#      rhoTag = cms.string('fixedGridRhoFastjetAll'),
#      useClosestToCentreSeedCrysDef = cms.bool(False),
#      maxRawEnergyForLowPtEBSigma = cms.double(-1),
#      maxRawEnergyForLowPtEESigma = cms.double(1200.),
#      eleRegs = cms.PSet(
#        ecalOnlyMean = cms.PSet(
#            rangeMinLowEt = cms.double(0.2),
#            rangeMaxLowEt = cms.double(2.0),
#            rangeMinHighEt = cms.double(-1.),
#            rangeMaxHighEt = cms.double(3.0),
#            forceHighEnergyTrainingIfSaturated = cms.bool(True),
#            lowEtHighEtBoundary = cms.double(999999.),
#            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
#            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
#            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
#            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
#            ),
#        ecalOnlySigma = cms.PSet(
#            rangeMinLowEt = cms.double(0.0002),
#            rangeMaxLowEt = cms.double(0.5),
#            rangeMinHighEt = cms.double(0.0002),
#            rangeMaxHighEt = cms.double(0.5),
#            forceHighEnergyTrainingIfSaturated = cms.bool(True),
#            lowEtHighEtBoundary = cms.double(999999.),
#            ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
#            ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
#            eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
#            eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
#            ),
#        epComb = cms.PSet(
#            ecalTrkRegressionConfig = cms.PSet(
#                rangeMinLowEt = cms.double(0.2),
#                rangeMaxLowEt = cms.double(2.0),
#                rangeMinHighEt = cms.double(0.2),
#                rangeMaxHighEt = cms.double(2.0),
#                lowEtHighEtBoundary = cms.double(999999.),
#                forceHighEnergyTrainingIfSaturated = cms.bool(False),
#                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
#                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
#                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
#                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
#                ),
#            ecalTrkRegressionUncertConfig = cms.PSet(
#                rangeMinLowEt = cms.double(0.0002),
#                rangeMaxLowEt = cms.double(0.5),
#                rangeMinHighEt = cms.double(0.0002),
#                rangeMaxHighEt = cms.double(0.5),
#                lowEtHighEtBoundary = cms.double(999999.),
#                forceHighEnergyTrainingIfSaturated = cms.bool(False),
#                ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
#                ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
#                eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
#                eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
#                ),
#            maxEcalEnergyForComb=cms.double(200.),
#            minEOverPForComb=cms.double(0.025),
#            maxEPDiffInSigmaForComb=cms.double(15.),
#            maxRelTrkMomErrForComb=cms.double(10.),
#            )
#        ),
#      phoRegs = regressionModifier106XUL.phoRegs.clone()
#    )
#
#)


process.ntuplizer = cms.EDAnalyzer('ElectronMVANtuplizer',
        # AOD case
        src                  = cms.InputTag('gedGsfElectrons'),
        vertices             = cms.InputTag('offlinePrimaryVertices'),
        pileup               = cms.InputTag('addPileupInfo'),
        genParticles         = cms.InputTag('genParticles'),
        # miniAOD case
        srcMiniAOD           = cms.InputTag('slimmedElectrons'),
        verticesMiniAOD      = cms.InputTag('offlineSlimmedPrimaryVertices'),
        pileupMiniAOD        = cms.InputTag('slimmedAddPileupInfo'),
        genParticlesMiniAOD  = cms.InputTag('prunedGenParticles'),
        #
        eleMVAs             = cms.untracked.vstring(
                                          ),
        eleMVALabels        = cms.untracked.vstring(
                                          ),
        eleMVAValMaps        = cms.untracked.vstring(
#                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values",
 #                                          "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values",
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values",
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values",
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues",
                                           ),
        eleMVAValMapLabels   = cms.untracked.vstring(
#                                          "Spring16GPVals",
#                                          "Spring16HZZVals",
                                           "Fall17IsoV2Vals",
                                           "Fall17NoIsoV2Vals",
                                           "BParkRetrainVals",
                                           ),
        eleMVACats           = cms.untracked.vstring(
                                           "electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Categories",
                                           ),
        eleMVACatLabels      = cms.untracked.vstring(
                                           "EleMVACats",
                                           ),
        #
        variableDefinition   = cms.string(mvaVariablesFile),
        isMC                 = cms.bool(True),
        deltaR               = cms.double(0.01),
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
	           ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
	           ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_mean"),
	           eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
	           eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_mean"),
	           ),
	       ecalOnlySigma = cms.PSet(
	           rangeMinLowEt = cms.double(0.0002),
	           rangeMaxLowEt = cms.double(0.5),
	           rangeMinHighEt = cms.double(0.0002),
	           rangeMaxHighEt = cms.double(0.5),
	           forceHighEnergyTrainingIfSaturated = cms.bool(True),
	           lowEtHighEtBoundary = cms.double(999999.),
	           ebLowEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
	           ebHighEtForestName = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
	           eeLowEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
	           eeHighEtForestName = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
	           ),
	       epComb = cms.PSet(
	           ecalTrkRegressionConfig = cms.PSet(
	               rangeMinLowEt = cms.double(0.2),
	               rangeMaxLowEt = cms.double(2.0),
	               rangeMinHighEt = cms.double(0.2),
	               rangeMaxHighEt = cms.double(2.0),
	               lowEtHighEtBoundary = cms.double(999999.),
	               forceHighEnergyTrainingIfSaturated = cms.bool(False),
	               ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
	               ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_mean'),
	               eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
	               eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_mean'),
	               ),
	           ecalTrkRegressionUncertConfig = cms.PSet(
	               rangeMinLowEt = cms.double(0.0002),
	               rangeMaxLowEt = cms.double(0.5),
	               rangeMinHighEt = cms.double(0.0002),
	               rangeMaxHighEt = cms.double(0.5),
	               lowEtHighEtBoundary = cms.double(999999.),
	               forceHighEnergyTrainingIfSaturated = cms.bool(False),
	               ebLowEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
	               ebHighEtForestName = cms.string('gsfElectron_eb_ecalTrk_05To50_sigma'),
	               eeLowEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
	               eeHighEtForestName = cms.string('gsfElectron_ee_ecalTrk_05To50_sigma'),
	               ),
	           maxEcalEnergyForComb=cms.double(200.),
	           minEOverPForComb=cms.double(0.025),
	           maxEPDiffInSigmaForComb=cms.double(15.),
	           maxRelTrkMomErrForComb=cms.double(10.),
	           )
	       ),
	     phoRegs = regressionModifier106XUL.phoRegs.clone()
	   )
        )



process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )
#process.regressionForEle= regressionForEle

process.RegEleSequence = cms.Sequence(
# process.regressionForEle
 process.egmGsfElectronIDSequence
)
#EleSequence =cms.Sequence(process.regressionForEle+process.egmGsfElectronIDSequence+process.ntuplizer)
process.RedoEle = cms.Path(process.RegEleSequence)
process.p = cms.Path(process.ntuplizer)

process.endjob_step = cms.EndPath(process.endOfProcess)
#process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process.schedule = cms.Schedule(
				process.RedoEle,
				process.p,
#                               process.nanoAOD_Kee_step, 
#                               process.nanoAOD_KstarMuMu_step,
#                               process.nanoAOD_KstarEE_step,
                                process.endjob_step, 
##                                process.NANOAODoutput_step
                               )


