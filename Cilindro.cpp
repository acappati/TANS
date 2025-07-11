#include "Cilindro.h"
#include "Retta.h"
#include "Punto.h"
#include "TMath.h"
#include <Riostream.h>


ClassImp(Cilindro)

  //implementazione di costruttori (default e standard) e distruttore
  Cilindro::Cilindro():TObject(),
 fRaggio(0.),fAltezza(0.)
    {
    }

  Cilindro::Cilindro(Double_t Raggio, Double_t Altezza):TObject(),
    fRaggio(Raggio),fAltezza(Altezza)
    {
    }

  Cilindro::~Cilindro()
    {
    }


//implementazione funzione Intersezione
//la funzione restituisce un bool TRUE solo se l'intersezione produce due valori t1 e t2 (doppia intersezione)
Bool_t Cilindro::Intersezione(Retta &retta,Punto &pintersezione)
{ 
   pintersezione.Reset();
   Double_t CosDir[3];
   retta.GetDirezione(CosDir);
   Punto p0 = retta.GetPunto();

   //calcolo discriminante
   Double_t a =CosDir[0]*CosDir[0]+CosDir[1]*CosDir[1];
   Double_t b =p0.GetX()*CosDir[0]+p0.GetY()*CosDir[1];
   Double_t c =p0.GetX()*p0.GetX()+p0.GetY()*p0.GetY()-fRaggio*fRaggio;
   Double_t Discriminante = b*b -a*c;
 
 if(Discriminante<0.)
      {
	//printf("Errore: retta esterna al cilindro"); debug
	return kFALSE;
      }
  
    Double_t t1=(-b-TMath::Sqrt(Discriminante))/a;
    Double_t t2=(-b+TMath::Sqrt(Discriminante))/a;
    Double_t prodotto=t1*t2;

	if(prodotto>0.)
	  { //printf("Errore: P0 della retta esterno al cilindro "); debug
	    return kFALSE;
	  }

	if(t2>0.)
	  { retta.SetParametro(t2); }
	else
	  { retta.SetParametro(t1); }

	//si definiscono il punto di intersezione e la sua coordinata z
	pintersezione=retta.GetPosizioneAttuale();
	Double_t z=pintersezione.GetZ();
	
	//si controlla se la z è interna ai rivelatori
       if(z<(-fAltezza/2.)|| z>(fAltezza/2.))
	  {	    
	    //printf("Errore: z del punto fuori dal range \n");
	    return kFALSE;
	  }
	
	return kTRUE;

} //fine funzione intersezione
