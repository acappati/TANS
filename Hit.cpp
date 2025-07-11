#include <Riostream.h>
#include "TMath.h"
#include "TObject.h"
#include "Punto.h"
#include "Retta.h" 
#include "Cilindro.h"
#include "TRandom3.h"
#include "Hit.h"

//________________________________________________________________________
Hit::Hit():TObject(),
 fZ(0.),
 fPhi(0.),
 fEtichetta(0)
 {
   // default constructor
 }

//___________________________________________________________________________
Hit::Hit(Double_t Z, Double_t Phi, Int_t Etichetta):TObject(),
 fZ(Z),
 fPhi(Phi),
 fEtichetta(Etichetta)
 {
  //standard constructor 
 }	     

//___________________________________________________________________________
Hit::~Hit()	 
{
  // destructor
}
