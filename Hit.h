#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "Punto.h"
#include "Retta.h"
#include "Cilindro.h"
#include "TRandom3.h"


class Hit: public TObject
{

 public:
  Hit();
  Hit(Double_t Z,Double_t Phi, Int_t Etichetta);
  virtual ~Hit();

  //funzioni Get
  Double_t GetZ() const {return fZ;}
  Double_t GetPhi() const {return fPhi;}
  Int_t GetEtichetta() const {return fEtichetta;}
  
  //funzioni Set
  void SetZ(Double_t Z) {fZ=Z;}
  void SetPhi(Double_t Phi) {fPhi=Phi;}
  void SetEtichetta(Int_t Etichetta) {fEtichetta=Etichetta;}
  
  void Reset()
	{ fZ=fPhi=0.;
	  fEtichetta=0;
	}

 private:
  Double_t fZ;
  Double_t fPhi;
  Int_t fEtichetta;  

  ClassDef(Hit,1)
  
};

#endif
