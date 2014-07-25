#ifndef VHbb_HbbProducer_HbbTuple_h
#define VHbb_HbbProducer_HbbTuple_h

#include "TLorentzVector.h"

#include <string>

namespace Hbb
{

  //---------------------------------------------------------------------------------

  struct Object
  {
    TLorentzVector lv;
    //Object *daughter1, *daughter2;

  Object() :
    lv(TLorentzVector())
    {
    }

  Object(TLorentzVector theLV) : lv(TLorentzVector())
    {
      lv.SetPtEtaPhiM(theLV.Pt(), theLV.Eta(), theLV.Phi(), theLV.M());
    }

  Object(double pT, double eta, double phi, double m) : lv(TLorentzVector())
    {
      lv.SetPtEtaPhiM(pT, eta, phi, m);
    }

  };
  
  //---------------------------------------------------------------------------------
  
  struct Jet:Object
  {
    float tau1, tau2, tau3;
    float prunedMass, trimmedMass, filteredMass;
    float qJetsVolatility;
  
  Jet() : Object(),
      tau1(-9999), tau2(-9999), tau3(-9999),
      prunedMass(-9999), trimmedMass(-9999), filteredMass(-9999),
      qJetsVolatility(-9999)
    {
    }

  };
  
  //---------------------------------------------------------------------------------

  struct Lepton:Object
  {
    int charge;
    
  Lepton() : Object(),
      charge(-9999)
      {
      }

  };

  //---------------------------------------------------------------------------------
  
  struct Electron:Lepton
  {
  Electron() : Lepton()
      {
      }

  };
  
  //---------------------------------------------------------------------------------
  
  struct Muon:Lepton
  {
  Muon() : Lepton()
      {
      }

  };

  //---------------------------------------------------------------------------------

  struct Tau:Lepton
  {
  Tau() : Lepton()
      {
      }

  };

  //---------------------------------------------------------------------------------

  struct V:Object
  {
  };

  //---------------------------------------------------------------------------------

  struct Higgs:Object
  {
  };

  //---------------------------------------------------------------------------------
  struct MET:Object
  {
  };

  //---------------------------------------------------------------------------------

  struct Tuple
  {
    int eventClassification;
    float rho;
    
    std::vector<Jet> AK4PFCHS;
    std::vector<Jet> AK8PFCHS;
    std::vector<Jet> AK10PFCHS;
    std::vector<Jet> AK12PFCHS;
    std::vector<Jet> AK15PFCHS;

    std::vector<Electron> Electrons;
    std::vector<Muon> Muons;
    std::vector<Tau> Taus;

    std::vector<Higgs> Higgses;
    
  Tuple() : 
    rho(-9999), 
      AK4PFCHS(std::vector<Jet>()), AK8PFCHS(std::vector<Jet>()), AK10PFCHS(std::vector<Jet>()), AK12PFCHS(std::vector<Jet>()), AK15PFCHS(std::vector<Jet>()),
      Electrons(std::vector<Electron>()), Muons(std::vector<Muon>()),Taus(std::vector<Tau>()),
      Higgses(std::vector<Higgs>())
    {
    }
  };
  
  typedef std::vector<Object> ObjectCollection;
  typedef std::vector<Jet> JetCollection;
  typedef std::vector<Electron> ElectronCollection;
  typedef std::vector<Muon> MuonCollection;
  typedef std::vector<Tau> TauCollection;
  typedef std::vector<Tuple> TupleCollection;
  typedef std::vector<Higgs> HiggsCollection;
}

#endif
