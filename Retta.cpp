#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Retta.h"
#include "TMatrixD.h"

ClassImp(Retta)

//////////////////////////
//COSTRUTTORI/////////////
//////////////////////////

// default constructor
Retta::Retta():TObject(),
  fP0(),
  fTheta(0.),
  fPhi(0.),
  fParametro(0.){ 
  fP0=Punto(0.,0.,0.);
  for(Int_t k=0;k<3;k++)
    {fCosDir[k]=0.;}  
 }

//standard constructor con angoli
Retta::Retta(Punto p, Double_t Theta, Double_t Phi):TObject(),
   fP0(p),fTheta(Theta),fPhi(Phi),fParametro(0.)
{
  DefCoseni();
}

//Distruttore
Retta::~Retta()
{
  
}

//standard constructor con coseni direttori
Retta::Retta(Punto p, Double_t CosDir[]):TObject(),
  fP0(p),fTheta(0.),fPhi(0.),fParametro(0.)
{
  UpdateDirezione(CosDir);
}

//////////////////////////////////////////////////////////////////////


//______________________________________________________________
void Retta::UpdateDirezione(Double_t Theta, Double_t Phi)
{
  fTheta=Theta;
  fPhi=Phi;
  DefCoseni();
}
//______________________________________________________________
void Retta::UpdateDirezione(Double_t VetCos[])
{
  for(Int_t i=0;i<3;i++)
    {fCosDir[i]=VetCos[i];}

  DefThetaPhi();
}
//______________________________________________________________
Punto Retta::GetPosizioneAttuale() const{
  
 Punto P(fP0.GetX()+fCosDir[0]*fParametro,fP0.GetY()+fCosDir[1]*fParametro, fP0.GetZ()+fCosDir[2]*fParametro);
 return P;
  
}


////////////////////////////////////////////////////////////////////

//__________funzioni private_________________

void Retta::DefCoseni()
{
  fCosDir[0]=TMath::Sin(fTheta)*TMath::Cos(fPhi);
  fCosDir[1]=TMath::Sin(fTheta)*TMath::Sin(fPhi);
  fCosDir[2]=TMath::Cos(fTheta);

}

void Retta::DefThetaPhi()
{
  fTheta=TMath::ACos(fCosDir[2]);
  fPhi=TMath::ATan2(fCosDir[1],fCosDir[0]);
  if(fCosDir[1]<0.) fPhi=2.*TMath::ACos(-1.)+fPhi;

}

