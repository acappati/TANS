#ifndef RETTA_H
#define RETTA_H

#include "TObject.h"
#include "Punto.h"

class Retta : public TObject
{
 public:
  Retta();
  Retta(Punto p, Double_t Theta, Double_t Phi);
  Retta(Punto p, Double_t CosDir[]);
  virtual ~Retta();

  Punto GetPunto() const {return fP0;}
  Punto GetPosizioneAttuale() const;
  void GetDirezione(Double_t &Theta, Double_t &Phi) const {Theta=fTheta; Phi=fPhi;}
  void GetDirezione(Double_t VetCos[]) const { for(Int_t i=0;i<3;i++) VetCos[i]=fCosDir[i];}
 
  
  void SetParametro(Double_t t) {fParametro=t;}
  void UpdateDirezione(Double_t Theta, Double_t Phi);
  void UpdateDirezione(Double_t VetCos[]);

  // UpdatePunto imposta come nuovo vertice della semiretta il punto di hit contro il layer
  void UpdatePunto(Punto p) {fP0=p;}   

 
 
 private:
  //sono metodi privati perche' non si deve mai cambiare gli angoli in modo indipendente dai coseni direttori
  void DefCoseni();
  void DefThetaPhi();

  Punto fP0;                //punto appartenente alla retta
  Double_t fTheta;          //angolo polare
  Double_t fPhi;            //angolo azimutale
  Double_t fParametro;      //! parametro della retta
  
  Double_t fCosDir[3];

  ClassDef(Retta,1)

};

#endif
