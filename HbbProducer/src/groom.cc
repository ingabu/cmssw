#include "VHbb/HbbProducer/interface/HbbTuple.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Njettiness.hh"

#include <vector>

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//---------------------------------------------------------------------------------------------------------------------------------
//grooming parameters

//trimming
double Rtrim = 0.2;
double ptfrac = 0.05;

//filtering
double Rfilt = 0.2;
double Nfilt = 3;

//pruning
double zcut=0.1;
double rcut_factor=0.5;

//n-subjettiness
double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means                                                    
double Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff) 
//---------------------------------------------------------------------------------------------------------------------------------

void groom(pat::Jet iJet, Hbb::Jet* oJet, double R){
  
  reco::Jet::Constituents constituents=iJet.getJetConstituents();

  vector<PseudoJet> fjConstituents;
  for(auto constituentItr=constituents.begin(); constituentItr!=constituents.end(); ++constituentItr){
    edm::Ptr<reco::Candidate> constituent=*constituentItr;

    fjConstituents.push_back(PseudoJet(constituent->px(), constituent->py(), constituent->pz(), constituent->energy()));
  }
  
  JetDefinition jet_def(cambridge_algorithm, R);
  // the use of a ClusterSequenceArea (instead of a plain ClusterSequence)
  // is only needed because we will later combine filtering with area-based
  // subtraction
  ClusterSequenceArea clust_seq(fjConstituents, jet_def, AreaDefinition(active_area_explicit_ghosts));
  vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(0));
  PseudoJet theJet=inclusive_jets[0];
  
  //------------------------------------
  // Trimming
  //------------------------------------
  
  Filter trimmer(JetDefinition(cambridge_algorithm, Rtrim), SelectorPtFractionMin(ptfrac) );
  PseudoJet trimmedJet = trimmer(theJet);
  
  oJet->trimmedMass=trimmedJet.m();

  vector<PseudoJet> fjSubjets = trimmedJet.pieces();
  for (size_t i = 0; i < fjSubjets.size(); i++){
    oJet->trimmedSubjets.push_back(new Hbb::Jet(fjSubjets[i].pt(), fjSubjets[i].eta(), fjSubjets[i].phi(), fjSubjets[i].m()));
  }

  //------------------------------------
  // Filtering
  //------------------------------------

  Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(Nfilt));
  PseudoJet filteredJet = filter(theJet);

  oJet->filteredMass=filteredJet.m();

  fjSubjets = filteredJet.pieces();
  for (size_t i = 0; i < fjSubjets.size(); i++){
    oJet->filteredSubjets.push_back(new Hbb::Jet(fjSubjets[i].pt(), fjSubjets[i].eta(), fjSubjets[i].phi(), fjSubjets[i].m()));
  }

  //------------------------------------
  //Pruning
  //------------------------------------

  Pruner pruner(cambridge_algorithm, zcut, rcut_factor);
  PseudoJet prunedJet = pruner(theJet);

  oJet->prunedMass=prunedJet.m();
  
  fjSubjets = prunedJet.pieces();
  for (size_t i = 0; i < fjSubjets.size(); i++){
    oJet->prunedSubjets.push_back(new Hbb::Jet(fjSubjets[i].pt(), fjSubjets[i].eta(), fjSubjets[i].phi(), fjSubjets[i].m()));
  }

  //------------------------------------
  // NSubJettiness
  //------------------------------------
  NsubParameters paraNsub(beta, R, Rcut);
  
  Njettiness nSubOnePass(Njettiness::onepass_kt_axes,paraNsub);
  
  oJet->tau1 = nSubOnePass.getTau(1,fjConstituents);
  oJet->tau2 = nSubOnePass.getTau(2,fjConstituents);
  oJet->tau3 = nSubOnePass.getTau(3,fjConstituents);
}
