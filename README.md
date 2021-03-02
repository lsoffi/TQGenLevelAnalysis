# TQGenLevelAnalysis

#setup release

scram project -n CMSSW-10218-GenAnalysis CMSSW CMSSW_10_2_18

cd CMSSW-10218-GenAnalysis/src

cmsenv

#clone github repo (better first to fork in your own github and then clone yours

git clone https://github.com/lsoffi/TQGenLevelAnalysis.git


#compile

cd TQGenLevelAnalysis/TQGenAnalyzer/

scram b


#run in local

cmsRun python/ConfFile_cfg.py




#usefule links:

https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py

https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule

https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
