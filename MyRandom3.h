#ifndef MYRANDOM3_H
#define MYRANDOM3_H
#include "TObject.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TAxis.h"

class Punto;
class Retta;
class Cilindro;
class Hit;

class MyRandom3 : public TObject {

 public:

  MyRandom3();
  MyRandom3(UInt_t seed);
  virtual ~MyRandom3();
    
  Punto EstraiVertice();
  Double_t EstraiPhi();
  Double_t EstraiTheta();
  Double_t EstraiThetaScattering(Double_t X0,Double_t x);
  UInt_t EstraiMolteplicita();
 
  void contohit(Punto P0, Hit &h1, Hit &h2, Bool_t &controlloL1, Bool_t &controlloL2, Bool_t ScatteringMultiploAttivo);
  void rotate(Double_t Theta,Double_t Phi,Double_t Theta_primo,Double_t Phi_primo,Double_t *CosDir);
  


 private:
  UInt_t fSeed;
  TH1F* fPseudorapidita;
  TH1F* fMolteplicita;
  
  ClassDef(MyRandom3,1) 
};



#endif
