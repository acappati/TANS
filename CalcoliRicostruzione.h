#ifndef CALCOLIRICOSTRUZIONE_H
#define CALCOLIRICOSTRUZIONE_H
#include "TObject.h"
#include <vector>
#include <utility>
#include "TH1D.h"
#include "TH1F.h"

using namespace std;

class CalcoliRicostruzione : public TObject {

 public:
  CalcoliRicostruzione(); 
  virtual ~CalcoliRicostruzione(); 

  //FUNZIONI PER RICOSTRUZIONE VERTICE DI COLLISIONE
  //________________________________________________

  //fornisce FWHM per criterio di matching
  Double_t CalcoloFWHM();
  //esegue intersezione tra tracklet e asse del fascio
  Double_t CalcoloCandidatiVertice(Double_t z1,Double_t z2); 
  //riempie istogramma con i possibili valori candidati della z del vertice
  void GeneraIstoCandidatiVertice(Double_t FWHM,TH1F *hist_vertice, vector <pair<Double_t,Double_t> > vL1,vector <pair<Double_t,Double_t> > vL2, vector <pair<Double_t, Double_t> >::iterator it1, vector <pair<Double_t, Double_t> >::iterator it2, vector<Double_t> &vCandidatiVertice);
  //a partire dall'istogramma ricostruisce z per ogni evento (ricostruzione avvenuta --> return kTRUE)
  Bool_t GetZRicostruita_Hist(Double_t &z, TH1F *hist, Int_t &Count_NonRec, vector<Double_t> vCandidatiVertice);

  

  //FUNZIONI PER GENERARE GRAFICI
  //_____________________________

  void GeneraIstoRisoluzione(TH1F *hist, vector<Double_t> vZVertice_Ricostruita, vector<Double_t> vZVertice_MC);
  void GraficoRisoluzioneVSMolteplicita(vector<Double_t> vZVertice_Ricostruita, vector<Double_t>  vZVertice_MC, vector<UInt_t> vMult);
  void GraficoRisoluzioneVSZGenerata(vector<Double_t> vZVertice_Ricostruita, vector<Double_t>  vZVertice_MC);
  void GraficoEfficienzaVSMolteplicita(vector <UInt_t> vMult_Tot, vector <UInt_t> vMult_Rec);

  private:
  
  ClassDef(CalcoliRicostruzione,1) 
};



#endif
