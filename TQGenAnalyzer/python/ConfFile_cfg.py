import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XUL
#from RecoEgamma.EgammaTools.regressionModifierNew_cfi import regressionModifier106XULLP


process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            'file:/afs/cern.ch/user/s/soffi/public/10C23D4F-94BD-E811-9588-E0071B7B2320.root'
#                            '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/854B1DC0-2F71-694D-A3F5-8DC1CDE1EF18.root'
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/file_Livia/01823A93-C4C6-4049-9A01-72D2BCA68238.root'
#running in local on psudoscalar signals old
                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/14GeV/14262656-E2BD-E811-B72F-20040FE8ECAC.root',
                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/14GeV/1A96C0D8-03BE-E811-8ED2-509A4C83EF52.root',
                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/14GeV/28930861-18C4-E811-A60C-008CFA111190.root',
                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/14GeV/28FC491B-7DB9-E811-8FDD-0CC47AA992AE.root',
                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/14GeV/10C23D4F-94BD-E811-9588-E0071B7B2320.root'

#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/18GeV/160B5697-46ED-E811-AD55-001E675A6D10.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/18GeV/367F64E8-94F1-E811-8F2D-B4969109F684.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/18GeV/AA8A78D1-ABE5-E811-9F72-0242AC130002.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/18GeV/DA652E66-D1EC-E811-A686-0025905C2CE8.root',
#                                    'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/18GeV/2C888440-D2EC-E811-9D6E-20CF307C98DC.root'

#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/22GeV/1E757B89-32B9-E811-9333-0CC47AD98D6C.root',
#                               'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/22GeV/24EAB9CE-C9B9-E811-B244-5065F3810301.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/22GeV/5A398598-28B9-E811-BFDB-00000086FE80.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/22GeV/367F64E8-94F1-E811-8F2D-B4969109F684.root'
#
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/26GeV/24AA2372-BBB9-E811-8CEB-0025904B7C40.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/26GeV/30467328-81B9-E811-B011-0CC47A4C8E2A.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/26GeV/6E864B97-BBB9-E811-BDD7-0CC47AA98D60.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/26GeV/8831B77A-BBB9-E811-A8CA-001E67E6F918.root',
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/26GeV/9052B279-BBB9-E811-9DFD-0242AC1C0502.root'


                )
                            )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_tetraquarks_pseudoscalar_14GeV.root")
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)

#for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] :
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff'] :
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.ntuplizer_seq = cms.Sequence()

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
#process.ntuplizer_seq *= process.lowPtGsfElectronID

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



process.GenAnalysis = cms.EDAnalyzer('TQGenAnalyzer',
                                     generatorInfo= cms.InputTag("generator"),
                                     prunedGenParticles    = cms.InputTag("prunedGenParticles"),
                                     patMuons = cms.InputTag("slimmedMuons"),
                                     patElectrons = cms.InputTag("slimmedElectrons"),                 # MINIAOD
                                     patElectronsLowPt = cms.InputTag("slimmedLowPtElectrons"), #change when running on our old signals
                                     vtx=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     drForCleaning = cms.double(0.03),
                                     dzForCleaning = cms.double(0.5), ##keep tighter dZ to check overlap of pfEle with lowPt (?)
                                     mvaValueEGamma = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values'),
                                     mvaIdEGamma = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90'), 
                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values'),
                                     mvaIdPF = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90'), 
                     #                mvaValue = cms.InputTag('lowPtGsfElectronID'),
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

process.ntuplizer_seq *= process.GenAnalysis

#process.p = cms.Path(process.egmGsfElectronIDSequence* process.ntuplizer_seq)

process.p = cms.Path(process.ntuplizer_seq)

process.schedule = cms.Schedule(process.p)
