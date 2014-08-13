// -*- C++ -*-
//
// Package:    VHbb/HbbProducer
// Class:      HbbProducer
// 
/**\class HbbProducer HbbProducer.cc VHbb/HbbProducer/plugins/HbbProducer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  john stupak
//         Created:  Fri, 27 Jun 2014 22:16:05 GMT
//
//


// system include files
#include <memory>
#include <vector>

// cmssw include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// user include files
#include "VHbb/HbbProducer/interface/HbbTuple.h"
#include "VHbb/HbbProducer/src/telescope.cc"
#include "VHbb/HbbProducer/src/groom.cc"

using namespace std;

typedef edm::Handle< edm::View< pat::Jet > > h_patJets;
typedef vector<Hbb::Jet> v_HbbJets;

//
// class declaration
//

class HbbProducer : public edm::EDProducer {
public:
  explicit HbbProducer(const edm::ParameterSet&);
  ~HbbProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  edm::InputTag _rhoSource;
  edm::InputTag _packedCandidateSource;
  edm::InputTag _AK4Source, _AK8Source, _AK10Source, _AK12Source, _AK15Source;
  edm::InputTag _AK8PackedSource, _AK10PackedSource, _AK12PackedSource, _AK15PackedSource;
  edm::InputTag _muonSource;
  
  Hbb::Tuple _output;
  
  
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HbbProducer::HbbProducer(const edm::ParameterSet& iConfig)
{
  
  _rhoSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("rhoSource"));

  _packedCandidateSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("packedCandidateSource"));

  _AK4Source=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK4Source"));
  _AK8Source=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK8Source"));
  _AK10Source=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK10Source"));
  _AK12Source=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK12Source"));
  _AK15Source=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK15Source"));

  _AK8PackedSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK8PackedSource"));
  _AK10PackedSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK10PackedSource"));
  _AK12PackedSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK12PackedSource"));
  _AK15PackedSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("AK15PackedSource"));

  _muonSource=edm::InputTag(iConfig.getParameter<edm::InputTag>("muonSource"));

  //register your products
  produces<Hbb::Tuple>();
}


HbbProducer::~HbbProducer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HbbProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<double> rho;
  iEvent.getByLabel(_rhoSource,rho);

  edm::Handle< edm::PtrVector<reco::Candidate> > packedCandidates;
  iEvent.getByLabel(_packedCandidateSource, packedCandidates);
  
  map< string, h_patJets > fatJetInputs=map< string, h_patJets >();
  map< string, v_HbbJets* > fatJetOutputs=map< string, v_HbbJets* >();

  h_patJets AK4jets;
  iEvent.getByLabel(_AK4Source,AK4jets);
  
  h_patJets AK8jets;
  iEvent.getByLabel(_AK8PackedSource,AK8jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK8"),AK8jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK8"),&_output.AK8PFCHS));

  h_patJets AK10jets;
  iEvent.getByLabel(_AK10PackedSource,AK10jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK10"),AK10jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK10"),&_output.AK10PFCHS));
  
  h_patJets AK12jets;
  iEvent.getByLabel(_AK12PackedSource,AK12jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK12"),AK12jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK12"),&_output.AK12PFCHS));
  
  h_patJets AK15jets;
  iEvent.getByLabel(_AK15PackedSource,AK15jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK15"),AK15jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK15"),&_output.AK15PFCHS));
  /*
  h_patJets AK8jets;
  iEvent.getByLabel(_AK8Source,AK8jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK8"),AK8jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK8"),&_output.AK8PFCHS));

  h_patJets AK10jets;
  iEvent.getByLabel(_AK10Source,AK10jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK10"),AK10jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK10"),&_output.AK10PFCHS));

  h_patJets AK12jets;
  iEvent.getByLabel(_AK12Source,AK12jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK12"),AK12jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK12"),&_output.AK12PFCHS));

  h_patJets AK15jets;
  iEvent.getByLabel(_AK15Source,AK15jets);
  fatJetInputs.insert(pair<string,h_patJets >(string("AK15"),AK15jets));
  fatJetOutputs.insert(pair<string, v_HbbJets* >(string("AK15"),&_output.AK15PFCHS));

  h_patJets AK8PackedJets;
  iEvent.getByLabel(_AK8PackedSource,AK8PackedJets);

  h_patJets AK10PackedJets;
  iEvent.getByLabel(_AK10PackedSource,AK10PackedJets);

  h_patJets AK12PackedJets;
  iEvent.getByLabel(_AK12PackedSource,AK12PackedJets);

  h_patJets AK15PackedJets;
  iEvent.getByLabel(_AK15PackedSource,AK15PackedJets);
  */

  edm::Handle< edm::View< pat::Muon > > inputMuons;
  iEvent.getByLabel(_muonSource,inputMuons);

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //Event quantities

  _output.rho=*rho;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //AK4 Jets

  pat::Jet d1, d2;
  double maxPT2=0;
  for(auto jet1=AK4jets->begin(); jet1!=AK4jets->end()-1; ++jet1){
    for(auto jet2=jet1+1; jet2!=AK4jets->end(); ++jet2){
      
      double pT2=pow(jet1->px()+jet2->px(),2)+pow(jet1->py()+jet2->py(),2);
      if (pT2>maxPT2){
	d1=*jet1;
	d2=*jet2;
      }
    }
  }

  vector<Hbb::Higgs> theHiggses=telescope(d1, d2, packedCandidates, iSetup);
  _output.Higgses=theHiggses;
  /*
  for(int i=0; i<nInterpretations; i++){
    float R=Rmin+(i*(Rmax-Rmin)/nInterpretations);

    Hbb::Higgs H;
    Hbb::Jet j1, j2;
    j1.lv.SetPxPyPzE(0,0,0,0);
    j2.lv.SetPxPyPzE(0,0,0,0);

    for(auto candidateItr=packedCandidates->begin(); candidateItr!=packedCandidates->end(); ++candidateItr){
      edm::Ptr<reco::Candidate> candidate=*candidateItr;

      TLorentzVector p4=TLorentzVector();
      p4.SetPxPyPzE(candidate->px(), candidate->py(), candidate->pz(), candidate->energy());

      double dR1Squared=deltaR2(d1, *candidate);
      double dR2Squared=deltaR2(d2, *candidate);

      if(dR1Squared<R*R || dR2Squared<R*R){
        if(dR1Squared<dR2Squared) j1.lv+=p4;
        else j2.lv+=p4;
      }
    }

    double delta=deltaR(d1,d2);
    double A=(TMath::Pi()*R*R) - (R*R*TMath::ACos(delta/(2*R))) + (delta*R*sqrt(1-pow(delta/(2*R),2)));
    j1.area=A;
    j2.area=A;
    

    double corrR=1;
    if (R<corrR) corrR=R;
    char c[99];
    sprintf(c, "ak%iPFchsL1L2L3", int(10*corrR));
    const JetCorrector* corrector = JetCorrector::getJetCorrector(string(c),iSetup);

    //JetCorrector objects only work with real jets - So create them

    XYZTLorentzVectorD bs_j1;
    bs_j1.SetPx(j1.lv.Px());
    bs_j1.SetPy(j1.lv.Py());
    bs_j1.SetPz(j1.lv.Pz());
    bs_j1.SetE(j1.lv.E());

    XYZTLorentzVectorD bs_j2;
    bs_j2.SetPx(j2.lv.Px());
    bs_j2.SetPy(j2.lv.Py());
    bs_j2.SetPz(j2.lv.Pz());
    bs_j2.SetE(j2.lv.E());

    pat::Jet p_j1;
    p_j1.setJetArea(A);
    p_j1.setP4(bs_j1);

    pat::Jet p_j2;
    p_j2.setJetArea(A);
    p_j2.setP4(bs_j2);

    double f1=corrector->correction(p_j1);
    j1.lv*=f1;

    double f2=corrector->correction(p_j2);
    j2.lv*=f2;
    
    H.lv=j1.lv+j2.lv;
    _output.Higgses.push_back(H);
    }
  */

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //fat jets

  for(auto fatJetCollection=fatJetInputs.begin(); fatJetCollection!=fatJetInputs.end(); ++fatJetCollection){
    string name=fatJetCollection->first;
    h_patJets inputJets=fatJetCollection->second;

    double radius=.1*std::atof(name.substr(2).c_str());

    v_HbbJets outputJets;
    for(auto inputJet=inputJets->begin(); inputJet!=inputJets->end(); ++inputJet){
      
      Hbb::Jet jet=Hbb::Jet(inputJet->pt(), inputJet->eta(), inputJet->phi(), inputJet->mass());
      
      groom(*inputJet, &jet, radius);

      outputJets.push_back(jet);
    }
    *fatJetOutputs[name]=outputJets;
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //Muons

  for(auto inputMuon=inputMuons->begin(); inputMuon!=inputMuons->end(); ++inputMuon){
    Hbb::Muon muon=Hbb::Muon();
    muon.lv.SetPtEtaPhiM(inputMuon->pt(), inputMuon->eta(), inputMuon->phi(), inputMuon->mass());
    muon.charge=inputMuon->charge();
    _output.Muons.push_back(muon);
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto_ptr<Hbb::Tuple> pOut(new Hbb::Tuple(_output));
  iEvent.put(pOut);

}
  
// ------------ method called once each job just before starting event loop  ------------
void 
HbbProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HbbProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  HbbProducer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
  void
  HbbProducer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  HbbProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  HbbProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HbbProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HbbProducer);
  
