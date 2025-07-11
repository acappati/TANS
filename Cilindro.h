#ifndef CILINDRO_H
#define CILINDRO_H
#include "TObject.h"

class Retta;
class Punto;

class Cilindro : public TObject
{

 public:
  Cilindro();
  Cilindro(Double_t Raggio,Double_t Altezza);  
  virtual ~Cilindro();

  //funzioni Get per accedere ai data member
  Double_t GetRaggio() const {return fRaggio;}
  Double_t GetAltezza() const {return fAltezza;}

  //funzione per il calcolo dell'intersezione tra retta e cilindro
  Bool_t Intersezione(Retta &retta,Punto &pintersezione);


 private:   
  Double_t fRaggio;
  Double_t fAltezza;


ClassDef(Cilindro,1)

};

#endif
