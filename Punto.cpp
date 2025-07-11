#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Punto.h"

ClassImp(Punto)

//________________________________________________________________________
Punto::Punto():TObject(),
 fX(0.),
 fY(0.),
 fZ(0.){
   // default constructor
 }


//___________________________________________________________________________
Punto::Punto(Double_t X, Double_t Y, Double_t Z):TObject(),
 fX(X),
 fY(Y),
 fZ(Z){
	//standard constructor 
}	     

//___________________________________________________________________________
Punto::~Punto()	 {
  // destructor
}

