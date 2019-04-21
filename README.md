# Run2VBSWGamma_ntuple

cmsrel CMSSW_8_0_32

cd CMSSW_8_0_32/src

cmsenv 

git cms-init

git cms-merge-topic cms-egamma:EGM_gain_v1

cd EgammaAnalysis/ElectronTools/data

git clone https://github.com/ECALELFS/ScalesSmearings.git

git checkout -b Moriond17_gainSwitch_unc/Moriond17_23Jan_v2

cd ../../../

sed -i '1ccorrectionType = "Moriond17_23Jan"' EgammaAnalysis/ElectronTools/python/calibrationTablesRun2.py

git cms-merge-topic lathomas:L1Prefiring_8_0_32

git cms-merge-topic cms-met:METRecipe_8020 -u

mv L1Prefiring/EventWeightProducer/files L1Prefiring/EventWeightProducer/data

git-cms-addpkg CommonTools/UtilAlgos/

git-cms-addpkg CommonTools/UtilAlgos/

git-cms-addpkg CondFormats/JetMETObjects

git-cms-addpkg DataFormats/Candidate/

git-cms-addpkg DataFormats/Common/   

git-cms-addpkg DataFormats/EgammaCandidates/

git-cms-addpkg DataFormats/HepMCCandidate/  

git-cms-addpkg DataFormats/Math/          

git-cms-addpkg DataFormats/PatCandidates/

git-cms-addpkg FWCore/Framework

git-cms-addpkg FWCore/ParameterSet

git-cms-addpkg FWCore/PluginManager

git-cms-addpkg FWCore/Utilities    

git-cms-addpkg GeneratorInterface/LHEInterface

git-cms-addpkg HLTrigger/HLTcore              

git-cms-addpkg JetMETCorrections/Algorithms

git-cms-addpkg JetMETCorrections/Modules   

git-cms-addpkg JetMETCorrections/Objects

git-cms-addpkg MuonAnalysis/MuonAssociators

git-cms-addpkg PhysicsTools/Heppy          

git-cms-addpkg PhysicsTools/PatAlgos

git-cms-addpkg PhysicsTools/PatExamples

git-cms-addpkg PhysicsTools/PatUtils   

git-cms-addpkg RecoEcal/EgammaCoreTools

git-cms-addpkg RecoEgamma/EgammaTools  

git-cms-addpkg SimDataFormats/JetMatching/



git clone https://github.com/qianminghuang/VBSWGamma.git VAJets

scram b (-j 8) 
