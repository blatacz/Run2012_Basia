
#include <sstream>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "HTTAnalyzer"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(&aDir, selectionFlavours_);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT& myEvent){

  float genWeight = myEvent.sampleWeight;
  
  std::string sampleName = "MC";
  if(genWeight==1) sampleName = "Data";

  if(fabs(fabs(myEvent.genDecay/24.0)-13)<1E-5 ||
     fabs(fabs(myEvent.genDecay/24.0)-15)<1E-5){
     sampleName = "WJets";
  }
  if(fabs(fabs(myEvent.genDecay/23.0)-13)<1E-5 ||
     fabs(fabs(myEvent.genDecay/23.0)-15)<1E-5){
    sampleName = "DY";
  }
  if(fabs(genWeight-0.0326672)<1E-5) sampleName = "TT";
  if(fabs(genWeight-0.0031596)<1E-5) sampleName = "Other";

  return sampleName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  
  float puWeight = myEventProxy.puWeight;
  float genWeight = myEventProxy.sampleWeight;
  float eventWeight = puWeight*genWeight;

  std::string sampleName = getSampleName(myEventProxy);

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  bool Selection = myEventProxy.ptL1>20 && myEventProxy.isPFMuon && myEventProxy.isTightMuon &&
    myEventProxy.ptL2>20 && myEventProxy.muFlag!=1 && myEventProxy.vetoEvent==0 &&
    myEventProxy.tightestHPSMVAWP>=0 && myEventProxy.pairIndex<1 && 
    myEventProxy.HLTx==1 && (myEventProxy.run>=163269 || myEventProxy.run==1) &&
    myEventProxy.HLTmatch==1;

  bool wSelection = Selection && myEventProxy.diTauCharge==0 && myEventProxy.MtLeg1MVA>60 && myEventProxy.combRelIsoLeg1DBetav2<0.1;

  bool qcdSelectionSS = Selection && (myEventProxy.diTauCharge==2 || myEventProxy.diTauCharge==-2);

  bool qcdSelectionOS = Selection && myEventProxy.diTauCharge==0; // && myEventProxy.MtLeg1MVA<40, wywalone aby mozna bylo ogladac obrazki ;

  bool baselineSelection = Selection && myEventProxy.diTauCharge==0 && myEventProxy.combRelIsoLeg1DBetav2<0.1;

  if(wSelection) {
	  ///Fill transverse mass
  	myHistos_->fill1DHistogram("h1DMt"+sampleName+"wsel",myEventProxy.MtLeg1MVA,eventWeight);
  }

  if(qcdSelectionSS){
	myHistos_->fill1DHistogram("h1DIso"+sampleName+"qcdselSS",myEventProxy.combRelIsoLeg1DBetav2,eventWeight);
  }

  if(qcdSelectionOS){
	myHistos_->fill1DHistogram("h1DIso"+sampleName+"qcdselOS",myEventProxy.combRelIsoLeg1DBetav2,eventWeight);
  }

// to receive QCD background
  if(qcdSelectionSS && myEventProxy.combRelIsoLeg1DBetav2<0.1){
	myHistos_->fill1DHistogram("h1DSVfit"+sampleName+"qcdselSS",myEventProxy.diTauNSVfitMass,eventWeight);
  }

  if(!baselineSelection) return true;
    
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.numPV,eventWeight);

  ///Fill SVfit mass
  myHistos_->fill1DHistogram("h1DSVfit"+sampleName,myEventProxy.diTauNSVfitMass,eventWeight);

  ///Fill transverse mass
  myHistos_->fill1DHistogram("h1DMt"+sampleName,myEventProxy.MtLeg1MVA,eventWeight);

  ///Fill muon pt histogram
  myHistos_->fill1DHistogram("h1DPtMu"+sampleName,myEventProxy.ptL1,eventWeight);

  ///Fill iso pt histogram
  myHistos_->fill1DHistogram("h1DIso"+sampleName,myEventProxy.ptL1,eventWeight);
  
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

