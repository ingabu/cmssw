#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TList.h"
#include "TKey.h"
#include "TObject.h"

//
// NOTE: Contains dirty hacks for things like 0P and 0PH so thy don't interfere.  not completely general!!
//

void makeNew(TFile * myFile, TString oldname, TString newname) {
  myFile->cd();
  TH1F* h1 = (TH1F*)myFile->Get(oldname);
  TH1F* h2 = (TH1F*)h1->Clone(newname);
  //cout << oldname << " " << h1->Integral() << endl;
  h2->Write();
}

void renormalize_to_other(TFile * myFile, TString integralHistName, TString changeHistName) {
  cout << "Renormalizing " << changeHistName << endl;
  myFile->cd();
  TH1F* hIntegral = (TH1F*)myFile->Get(integralHistName);
  TH1F* hChange = (TH1F*)myFile->Get(changeHistName);
  hChange->Scale(hIntegral->Integral()/hChange->Integral());
  hChange->Write();
}

void renormalize_by_ratio_of_hists(TFile * myFile, TString integralHistName_Numerator, TString integralHistName_Denominator, TString changeHistName) {
  cout << "Renormalizing " << changeHistName << endl;
  myFile->cd();
  TH1F* hIntegral_N = (TH1F*)myFile->Get(integralHistName_Numerator);
  //cout << hIntegral_N << endl;
  TH1F* hIntegral_D = (TH1F*)myFile->Get(integralHistName_Denominator);
  //cout << hIntegral_D << endl;
  TH1F* hChange = (TH1F*)myFile->Get(changeHistName);
  //cout << hChange << endl;
  hChange->Scale(hIntegral_N->Integral()/hIntegral_D->Integral());
  hChange->Write();
}

void renormalize_by_factor(TFile * myFile, double factor, TString changeHistName) {
  cout << "Renormalizing " << changeHistName;
  myFile->cd();
  TH1F* hChange = (TH1F*)myFile->Get(changeHistName);
  cout << setprecision(15) << " from " << hChange->Integral();
  hChange->Scale(factor);
  cout << setprecision(15) << " to " << hChange->Integral() << endl;
  hChange->Write();
}

void convertRootFile(TString sigName = "Wh_125p6_0P", TString sigAltName = "Wh_125p6_0M", double altMu=-1.0) {

  TFile * myFile = new TFile("plots.root", "update");
  cout << "Changing file " << myFile << endl;
  cout << sigName << " is sig, and " << sigAltName << " is sig_ALT" << endl;
  cout << endl;

  TList * list = myFile->GetListOfKeys();
  for (int i = 0; i < list->GetSize(); i++) {
    TKey *key = dynamic_cast<TKey*>(list->At(i));
    TObject* obj = key->ReadObj();

    TString oldName = obj->GetName(); 
    TString newName = oldName;
    
    if( oldName.Contains(sigName) ) {
      if( !(sigName=="Wh_125p6_0P" && oldName.Contains("0PH")) ){ //dirty hack
	newName.ReplaceAll(sigName,"sig");
	cout << oldName << " -> " << newName << endl;
	makeNew(myFile, oldName, newName);
      }//dirty hack
    }
    if( oldName.Contains(sigAltName) ) {
      if( !(sigAltName=="Wh_125p6_0M" && oldName.Contains("0Mf05ph0")) ){ //dirty hack
	newName.ReplaceAll(sigAltName,"sig_ALT");
	cout << oldName << " -> " << newName << endl;
	makeNew(myFile, oldName, newName);
	
	//Renormalize by factor
	if(altMu >= 0) {
	  renormalize_by_factor(myFile, altMu, newName);
	}
      }//dirty hack
    }
    
  }//end loop over keys
  
  myFile->Close();

}

