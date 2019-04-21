
/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Chayanit Asawatangtrakuldee chayanit@cern.ch 
 *
 * Description:
 *   - Selects "loose" and "tight" electrons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of electrons passing the 
 *     required electron ID.
 * History:
 *   
 *
 *****************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class ElectronIdSelector : public edm::EDProducer
{
public:
  // construction/destruction
  ElectronIdSelector(const edm::ParameterSet& iConfig);
  virtual ~ElectronIdSelector();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:  
  // member data
  //edm::InputTag  src_;
  std::string    moduleLabel_;
  std::string    idLabel_;  
  bool           useDetectorIsolation_;
  bool           applyTightID_;
  bool           applyMediumID_;
  bool           applyLooseID_;
  bool           applyVetoID_;
  unsigned int nTot_;
  unsigned int nPassed_;
  edm::EDGetTokenT<pat::ElectronCollection> ElectronToken_;
  edm::EDGetTokenT<reco::VertexCollection> VertexToken_;
  edm::EDGetTokenT<double> RhoToken_;
  EffectiveAreas effectiveAreas_;
};



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
ElectronIdSelector::ElectronIdSelector(const edm::ParameterSet& iConfig)
//  : src_    (iConfig.getParameter<edm::InputTag>     ("src"))
  : moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  , idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  , useDetectorIsolation_(iConfig.existsAs<bool>("useDetectorIsolation") ? iConfig.getParameter<bool>("useDetectorIsolation") : false)
  , nTot_(0)
  , nPassed_(0)
  , ElectronToken_ (consumes<pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>( "src" ) ) ) 
  , VertexToken_ (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>( "vertex" ) ) )
  , RhoToken_ (consumes<double> (iConfig.getParameter<edm::InputTag>( "rho") ) )
  , effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{
  produces<std::vector<pat::Electron> >();

  /// ------- Decode the ID criteria --------
  applyTightID_ = false;
  applyMediumID_ = false;
  applyLooseID_ = false;
  applyVetoID_ = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) ||
      (idLabel_.compare("WP70")==0) ||
      (idLabel_.compare("wp70")==0) )  
    applyTightID_ = true;
  else if( (idLabel_.compare("medium")==0) ||
      (idLabel_.compare("Medium")==0) ||
      (idLabel_.compare("MEDIUM")==0) ||
      (idLabel_.compare("WP80")==0) ||
      (idLabel_.compare("wp80")==0) )  applyMediumID_ = true;
  else if( (idLabel_.compare("loose")==0) || 
      (idLabel_.compare("Loose")==0) || 
      (idLabel_.compare("LOOSE")==0) ||
      (idLabel_.compare("WP90")==0) ||
      (idLabel_.compare("wp90")==0) )  applyLooseID_ = true;
  else if( (idLabel_.compare("veto")==0) || 
      (idLabel_.compare("Veto")==0) || 
      (idLabel_.compare("VETO")==0) ||
      (idLabel_.compare("VETOid")==0) ||
      (idLabel_.compare("VetoId")==0) )  applyVetoID_ = true;
}

 
//______________________________________________________________________________
ElectronIdSelector::~ElectronIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void ElectronIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

  // edm::Handle<reco::VertexCollection> vtxs;
  //iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtxs);
   edm::Handle<reco::VertexCollection> vtxs;
   iEvent.getByToken(VertexToken_, vtxs);
//  edm::Handle<reco::ConversionCollection> conversions;
//  iEvent.getByLabel("allConversions", conversions);

//  edm::Handle<reco::BeamSpot> beamspot_h;
//  iEvent.getByLabel("offlineBeamSpot", beamspot_h);
//  const reco::BeamSpot &beamspot = *(beamspot_h.product());

  std::auto_ptr<std::vector<pat::Electron> > passingElectrons(new std::vector<pat::Electron >);

//  edm::Handle<edm::View<pat::Electron> > electrons;
//  iEvent.getByLabel(src_,electrons);
  
   edm::Handle<pat::ElectronCollection > electrons;
   iEvent.getByToken(ElectronToken_, electrons);  

  bool* isPassing = new bool[electrons->size()];

  double rhoVal_;
  rhoVal_=-99.;
  edm::Handle<double> rho;
  iEvent.getByToken(RhoToken_,rho);
  rhoVal_ = *rho;

  for(unsigned int iElec=0; iElec<electrons->size(); iElec++) { 

    isPassing[iElec]=false;

    const pat::Electron& ele = electrons->at(iElec);

    // -------- Make sure that the electron is within acceptance ------
    float eta = ele.superCluster()->eta();
    bool isEB = ele.isEB() && fabs(eta) < 1.479;
    bool isEE = ele.isEE() && fabs(eta) > 1.479 && fabs(eta) < 2.5;
    //bool inAcceptance = (isEB || isEE) && (ele.ecalDrivenSeed()==1);
    float pt  = ele.pt();

    // -------- Compute Detector isolation ------
//    const double PI = 4.0*atan(1.0);
//    float detector_isolation = (ele.dr03TkSumPt() + 
//			       std::max(0.,ele.dr03EcalRecHitSumEt()-1.0) + 
//			       ele.dr03HcalTowerSumEt() - 
//			       PI*0.3*0.3*rhoVal_) / pt;

//    float EffArea = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , eta , ElectronEffectiveArea::kEleEAData2012);
//    float pfIso03EA = (ele.chargedHadronIso() + std::max(0.,ele.neutralHadronIso() + ele.photonIso() - EffArea*rhoVal_)) / pt;

    float eA = effectiveAreas_.getEffectiveArea(fabs(eta));
    float pf_isolation = ( ( ele.pfIsolationVariables().sumChargedHadronPt
                                   + std::max( 0., ele.pfIsolationVariables().sumNeutralHadronEt + ele.pfIsolationVariables().sumPhotonEt - eA*rhoVal_) )
                                 / pt );

    float isolation = 100.;
    isolation = pf_isolation;

//    if(useDetectorIsolation_) isolation = detector_isolation;
//    else isolation = pfIso03EA;

    // -------- Compute ID ------
    double sigmaIEtaIEta   = ele.full5x5_sigmaIetaIeta();
    double dPhiIn    = fabs(ele.deltaPhiSuperClusterTrackAtVtx());
//    double dEtaIn    = fabs(ele.deltaEtaSuperClusterTrackAtVtx());
    double dEtaIn    = fabs(ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ? ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max());
    double hoe     = ele.hadronicOverEm();
    double ooemoop = fabs((1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy()));

    // impact parameter variables
    float d0vtx         = 0.0;
    float dzvtx         = 0.0;
    if (vtxs->size() > 0) {
        reco::VertexRef vtx(vtxs, 0);    
        d0vtx = ele.gsfTrack()->dxy(vtx->position());
        dzvtx = ele.gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = ele.gsfTrack()->dxy();
        dzvtx = ele.gsfTrack()->dz();
    }

    // conversion rejection variables
//    bool vtxFitConversion = ConversionTools::hasMatchedConversion(ele, conversions, beamspot.position());
      bool vtxFitConversion = !(ele.passConversionVeto()==1);
//    float mHits = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    float mHits=ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);  

    bool isTight  = false;  /////// <--- equivalent to WP70
    bool isMedium = false;  /////// <--- equivalent to WP80
    bool isLoose  = false;  /////// <--- equivalent to WP90
    bool isVeto    = false;  /////// <--- the loosest cut for veto

    // ---------- cut-based ID -----------------
      ///Spring15
/*
    isTight = (pt>20.)  &&
      (!vtxFitConversion) &&
      ((isEB && mHits<=2 && isolation<0.0354 && sigmaIEtaIEta<0.0101 && dPhiIn<0.0336 && dEtaIn<0.00926 && hoe<0.0597 && ooemoop<0.012 && fabs(d0vtx)<0.0111 && fabs(dzvtx)<0.0466 )  ||
      (isEE && mHits<=1 && isolation<0.0646 && sigmaIEtaIEta<0.0279 && dPhiIn<0.0918 && dEtaIn<0.00724 && hoe<0.0615 && ooemoop<0.00999 && fabs(d0vtx)<0.0351 && fabs(dzvtx)<0.417));

    isMedium = (pt>20.)  &&
        (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.0766 && sigmaIEtaIEta<0.0101 && dPhiIn<0.0336 && dEtaIn<0.0103 && hoe<0.0876 && ooemoop<0.0174 && fabs(d0vtx)<0.0118 && fabs(dzvtx)<0.373) ||
         (isEE && mHits<=1 && isolation<0.0678 && sigmaIEtaIEta<0.0283 && dPhiIn<0.114 && dEtaIn<0.00733 && hoe<0.0678 && ooemoop<0.0898 && fabs(d0vtx)<0.0739 && fabs(dzvtx)<0.602));

    isLoose = (pt>20.)  &&
         (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.0893 && sigmaIEtaIEta<0.0103 && dPhiIn<0.115  && dEtaIn<0.0105  && hoe<0.104 && ooemoop<0.102 && fabs(d0vtx)<0.0261 && fabs(dzvtx)<0.41) ||
         (isEE && mHits<=1 && isolation<0.121 && sigmaIEtaIEta<0.0301  && dPhiIn<0.182 && dEtaIn<0.00814 && hoe<0.0897 && ooemoop<0.126 && fabs(d0vtx)<0.118 && fabs(dzvtx)<0.822));

    isVeto = (pt>20.) && 
         (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.126 && sigmaIEtaIEta<0.0114 && dPhiIn<0.216 && dEtaIn<0.0152 && hoe<0.181 && ooemoop<0.207  && fabs(d0vtx)<0.0564 && fabs(dzvtx)<0.472 ) ||
         (isEE && mHits<=3 && isolation<0.144 && sigmaIEtaIEta<0.0352 && dPhiIn<0.237 && dEtaIn<0.0113 && hoe<0.116 && ooemoop<0.174  && fabs(d0vtx)<0.222 && fabs(dzvtx)<0.921));
*/

	//2016 data
//    isTight = (pt>20.)  &&
//      (!vtxFitConversion) &&
//      ((isEB && mHits<=1 && isolation<0.0588 && sigmaIEtaIEta<0.00998 && dPhiIn<0.0816 && dEtaIn<0.00308 && hoe<0.0414 && ooemoop<0.0129 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 )  || (isEE && mHits<=1 && isolation<0.0571 && sigmaIEtaIEta<0.0292 && dPhiIn<0.0394 && dEtaIn<0.00605 && hoe<0.0641 && ooemoop<0.0129 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

//fail tight ID, pass veto ID
    isTight = (pt>20.)  &&
      (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.175 && sigmaIEtaIEta<0.0115 && dPhiIn<0.228 && dEtaIn<0.00749 && hoe<0.356 && ooemoop<0.299  && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 && !(mHits<=1 && isolation<0.0588 && sigmaIEtaIEta<0.00998 && dPhiIn<0.0816 && dEtaIn<0.00308 && hoe<0.0414 && ooemoop<0.0129)) ||
         (isEE && mHits<=3 && isolation<0.159 && sigmaIEtaIEta<0.037 && dPhiIn<0.213 && dEtaIn<0.00895 && hoe<0.211 && ooemoop<0.15  && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20 && !(mHits<=1 && isolation<0.0571 && sigmaIEtaIEta<0.0292 && dPhiIn<0.0394 && dEtaIn<0.00605 && hoe<0.0641 && ooemoop<0.0129)));

    isMedium = (pt>20.)  &&
        (!vtxFitConversion) &&
        ((isEB && mHits<=1 && isolation<0.0695 && sigmaIEtaIEta<0.00998 && dPhiIn<0.103 && dEtaIn<0.00311 && hoe<0.253 && ooemoop<0.134 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         (isEE && mHits<=1 && isolation<0.0821 && sigmaIEtaIEta<0.0298 && dPhiIn<0.045 && dEtaIn<0.00609 && hoe<0.0878 && ooemoop<0.13 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    isLoose = (pt>20.)  &&
         (!vtxFitConversion) &&
        ((isEB && mHits<=1 && isolation<0.0994 && sigmaIEtaIEta<0.011 && dPhiIn<0.222  && dEtaIn<0.00477  && hoe<0.298 && ooemoop<0.241 && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10) ||
         (isEE && mHits<=1 && isolation<0.107 && sigmaIEtaIEta<0.0314  && dPhiIn<0.213 && dEtaIn<0.00868 && hoe<0.101 && ooemoop<0.14 && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));

    isVeto = (pt>20.) && 
         (!vtxFitConversion) &&
        ((isEB && mHits<=2 && isolation<0.175 && sigmaIEtaIEta<0.0115 && dPhiIn<0.228 && dEtaIn<0.00749 && hoe<0.356 && ooemoop<0.299  && fabs(d0vtx)<0.05 && fabs(dzvtx)<0.10 ) ||
         (isEE && mHits<=3 && isolation<0.159 && sigmaIEtaIEta<0.037 && dPhiIn<0.213 && dEtaIn<0.00895 && hoe<0.211 && ooemoop<0.15  && fabs(d0vtx)<0.10 && fabs(dzvtx)<0.20));


    /// ------- Finally apply selection --------
    if(applyTightID_ && isTight)   isPassing[iElec]= true;
    if(applyMediumID_ && isMedium) isPassing[iElec]= true;
    if(applyLooseID_ && isLoose)   isPassing[iElec]= true;
    if(applyVetoID_ && isVeto) isPassing[iElec]= true;
    
 }
  

 for (unsigned int iElectron = 0; iElectron < electrons -> size(); iElectron ++)
   {     if(isPassing[iElectron]) passingElectrons->push_back( electrons -> at(iElectron) );       
  }
 
/*  unsigned int counter=0;
  edm::View<pat::Electron>::const_iterator tIt, endcands = electrons->end();
  for (tIt = electrons->begin(); tIt != endcands; ++tIt, ++counter) {
    if(isPassing[counter]) passingElectrons->push_back( *tIt );  
  }
*/
  nTot_  +=electrons->size();
  nPassed_+=passingElectrons->size();

  delete [] isPassing;  
  iEvent.put(passingElectrons);
}

 
//______________________________________________________________________________
void ElectronIdSelector::endJob()
{
  std::stringstream ss;
  ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   <<"\n"<<moduleLabel_<<"(ElectronIdSelector) SUMMARY:\n"<<ss.str()
	   <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef ElectronIdSelector   			    PATElectronIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATElectronIdSelector);
