#include "TRandom3.h"
#include "MyRandom3.h"
#include "Punto.h"
#include "Retta.h"
#include "Cilindro.h"
#include "Hit.h"
#include "TMath.h"
#include <Riostream.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TCanvas.h"

class Hit;

ClassImp(MyRandom3)

  //DATI DEL PROBLEMA
  //spessori rivelatori e lunghezze di radiazione (cm)
  Double_t X0_berillio=35.28;
  Double_t x_berillio=0.08;  
  Double_t X0_silicio=9.36;  
  Double_t x_silicio=0.02; 
 
  //oggetti della classe Cilindro: beampipe, rivelatore R1 ed R2 (raggio e altezza in cm)
  Cilindro beampipe=Cilindro(3.,100.); 
  Cilindro L1=Cilindro(4.,27.);
  Cilindro L2=Cilindro(7.,27.);

//__________________________________________________________________
MyRandom3::MyRandom3():TObject(),
  fSeed(0),fPseudorapidita(0),fMolteplicita(0)
  {
    // Default constructor
  }


//___________________________________________________________________
MyRandom3::MyRandom3(UInt_t seed):TObject(),
  fSeed(seed),fPseudorapidita(),fMolteplicita()
 {
   // Standard constructor

   //viene impostato il seed
   gRandom->SetSeed(seed); 
 
   //si legge il file e si inizializzano i data member
   TFile F("kinem.root");
   fPseudorapidita=dynamic_cast<TH1F*>(F.Get("heta"));
   fMolteplicita=dynamic_cast<TH1F*>(F.Get("hmul")); 
   fPseudorapidita->SetDirectory(0);
   fMolteplicita->SetDirectory(0);
   F.Close();
 }


//___________________________________________________________________
MyRandom3::~MyRandom3()
 {
  // Default destructor
  delete fPseudorapidita;
  fPseudorapidita=NULL;
  delete fMolteplicita;
  fMolteplicita=NULL;  
 }

//___________________________________________________________________
Punto MyRandom3::EstraiVertice()
{    
     Double_t Sigmax=0.01; //cm
     Double_t Sigmay=0.01;
     Double_t Sigmaz=5.3;
     Double_t Xvertice, Yvertice, Zvertice;     
     Xvertice=gRandom->Gaus(0.,Sigmax);
     Yvertice=gRandom->Gaus(0.,Sigmay);
     //Zvertice=-13.5 + 27*gRandom->Rndm(); 
     Zvertice=gRandom->Gaus(0.,Sigmaz); 
     Punto Vertice(Xvertice, Yvertice, Zvertice);
     return Vertice;
}
 
//_________________________________________________________________
Double_t MyRandom3::EstraiPhi()
{
  Double_t pi=TMath::ACos(-1.);
  Double_t Phi=2.*pi*gRandom->Rndm();
  return Phi;
}


//_________________________________________________________________
Double_t MyRandom3::EstraiTheta()
{  
  Double_t Eta;
  do{Eta=fPseudorapidita->GetRandom();}
  while(Eta<-2. || Eta>2.);

  Double_t Theta=2*TMath::ATan(TMath::Exp(-Eta));
  return Theta;
}


//________________________________________________________________
Double_t MyRandom3::EstraiThetaScattering(Double_t X0,Double_t x)
{
  Double_t Theta0= 0.0136*TMath::Sqrt(x/X0)*(1+0.038*TMath::Log(x/X0));
  Double_t Theta=gRandom->Gaus(0.,Theta0);
  return Theta;  
}


//_________________________________________________________________
UInt_t MyRandom3::EstraiMolteplicita()
{
  UInt_t Molteplicita=fMolteplicita->GetRandom();
  return Molteplicita;
}


//_________________________________________________________________
void MyRandom3::rotate(Double_t Theta,Double_t Phi,Double_t Theta_primo,Double_t Phi_primo,Double_t *CosDir)
{
  
    Double_t mat[3][3]; 
    Double_t pi=TMath::ACos(-1.);

    mat[0][0]=-TMath::Sin(Phi);
    mat[1][0]=TMath::Cos(Phi);
    mat[2][0]=0.;
    mat[0][1]=-TMath::Cos(Phi)*TMath::Cos(Theta);
    mat[1][1]=-TMath::Cos(Theta)*TMath::Sin(Phi);
    mat[2][1]=TMath::Sin(Theta);
    mat[0][2]=TMath::Sin(Theta)*TMath::Cos(Phi);
    mat[1][2]=TMath::Sin(Theta)*TMath::Sin(Phi);
    mat[2][2]=TMath::Cos(Theta);

    Double_t CosDir_primo[3];

    CosDir_primo[0]=TMath::Sin(Theta_primo)*TMath::Cos(Phi_primo);
    CosDir_primo[1]=TMath::Sin(Theta_primo)*TMath::Sin(Phi_primo);
    CosDir_primo[2]=TMath::Cos(Theta_primo);

    for(Int_t i=0;i<3;i++)
        {
	  CosDir[i]=static_cast<Double_t>(TMath::Sin(pi)); 
	  for(Int_t j=0;j<3;j++)
	    { CosDir[i]+=mat[i][j]*CosDir_primo[j];  }
	}    
}//fine funzione rotate



//******************************************************************************************************************
//la funzione "contohit" propaga la particella, fa scattering (se si vuole) e calcola intersezioni con i rivelatori
//******************************************************************************************************************
void MyRandom3::contohit(Punto P0, Hit &h1, Hit &h2, Bool_t &controlloL1, Bool_t &controlloL2, Bool_t ScatteringMultiploAttivo)
{
  h1.Reset();
  h2.Reset();
  controlloL1=kFALSE;
  controlloL2=kFALSE;

  //SI CREANO GLI OGGETTI DEL PROBLEMA  
  //oggetti della classe Punto: punti di hit con gli oggetti Cilindro
  Punto Phitbeam;
  Punto Phit1;
  Punto Phit2;
  //oggetti della classe retta (serve solo per chiamare le funzioni)
  Retta retta;
  //angoli theta e phi iniziali
  Double_t Theta;
  Double_t Phi;
  //vettore coseni direttori
  Double_t CosDir[3];   

  //estrazione degli angoli theta e phi 
  Phi=EstraiPhi();
  Theta=EstraiTheta(); 
  
  //creazione e inizializzazione oggetto della classe Retta
  retta=Retta(P0,Theta,Phi); 

  ////////////////////////////////
  // intersezione con beampipe  //
  ////////////////////////////////
   beampipe.Intersezione(retta,Phitbeam);
 
  ///////////////////////////////////////////////
  // eventuale scattering multiplo su beampipe //
  ///////////////////////////////////////////////
  //si estraggono angoli in SR particella e li si ruota in SR lab
  //si aggiornano i data member
  //OSS: se scattering multiplo è silenziato, si mantengono gli angoli Theta e Phi inalterati e si procede con il calcolo dell'intersezione
   
   if(ScatteringMultiploAttivo==kTRUE)
     {
       //scattering multiplo attivo
       //si estraggono gli angoli nuovi nel SR della particella
       //Phitemp tra 0 e 2Pi, Thetatemp con distribuzione gaussiana di dev.st. Theta0 per il berillio      
       Double_t Phitemp=EstraiPhi();
       Double_t Thetatemp=EstraiThetaScattering(X0_berillio,x_berillio);

       //rotate: scrive angoli rispetto al SR del laboratorio e aggiorna i coseni direttori
       rotate(Theta,Phi,Thetatemp,Phitemp,CosDir);
       
       //update direzione: ha come argomento i coseni direttori e aggiorna gli angoli
       retta.UpdateDirezione(CosDir);
       //calcolo intersezione della retta con la beampipe
       retta.UpdatePunto(Phitbeam);
     }

   else
     {
       //scattering multiplo non attivo
       //non vengono estratti angoli nuovi
       //si calcola solo intersezione della retta con la beampipe
       retta.UpdatePunto(Phitbeam);      
     }

  //si definiscono due nuovi angoli Phi1 e Theta1, angoli nuovi nel SR del LAB
  //la funzione GetDirezione copia i valori dei data member sulle nuove variabili theta1 e phi1
  //OSS se lo scattering multiplo è spento gli angoli rimangono invariati
  Double_t Theta1;
  Double_t Phi1;
  retta.GetDirezione(Theta1,Phi1);
  


  ////////////////////////////////
  //  intersezione con LAYER 1  //
  ////////////////////////////////
  controlloL1=L1.Intersezione(retta,Phit1); 
  Double_t Z1=Phit1.GetZ(); 

  //oggetto hit h1 (si usano le funzioni Set)
  h1.SetZ(Z1);
  h1.SetPhi(Phi1);
  h1.SetEtichetta(1);
 
   

  ///////////////////////////////////////////////
  //  eventuale scattering multiplo su LAYER1  //
  ///////////////////////////////////////////////
  //si estraggono angoli in SR particella e li si ruota in SR lab, infine si aggiornano data member
 
  if(ScatteringMultiploAttivo==kTRUE)
     {  
       //scattering multiplo attivo: estraggo gli angoli nuovi nel SR della particella
       Double_t Phitemp1=EstraiPhi();
       Double_t Thetatemp1=EstraiThetaScattering(X0_silicio,x_silicio);
       rotate(Theta1,Phi1,Thetatemp1,Phitemp1,CosDir);
       retta.UpdateDirezione(CosDir);
       retta.UpdatePunto(Phit1);
     }
  else
    {
      //scattering multiplo non attivo
      retta.UpdatePunto(Phit1);
    }
  
  Double_t Theta2;
  Double_t Phi2;
  retta.GetDirezione(Theta2,Phi2);


  ////////////////////////////////
  //  intersezione con LAYER 2  //
  ////////////////////////////////
  controlloL2=L2.Intersezione(retta,Phit2);
  Double_t Z2=Phit2.GetZ();
  
  h2.SetZ(Z2);
  h2.SetPhi(Phi2);
  h2.SetEtichetta(1);
 
}//fine funzione contohit
