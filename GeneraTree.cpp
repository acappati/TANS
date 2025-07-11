#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TObject.h"
#include "Punto.h"
#include "Retta.h" 
#include "Cilindro.h"
#include "TRandom3.h"
#include "MyRandom3.h"
#include "Hit.h"
#include "TH1F.h"
#include "TStopwatch.h"

//    !     !     !     !    !     !     !     !     !     !     !     !     !
//scegliere se attivare lo scattering multiplo (kTRUE) o disattivarlo (kFALSE)
//    !     !     !     !    !     !     !     !     !     !     !     !     !
Bool_t ScatteringMultiploAttivo=kTRUE;


//costanti globali
Double_t pi=TMath::ACos(-1.);
//raggi dei rivelatori
Double_t Rbeam=3.;
Double_t R1=4.;
Double_t R2=7.;
Double_t LunghezzaRiv=27.;
//si imposta il seed chiamando il costruttore della classe MyRandom
MyRandom3 pippo=MyRandom3(2126);



//Si crea e si riempie un TTree con vertice, molteplicità e hits di ogni evento di collisione

void GeneraTree(UInt_t NumeroPuntiRumore=0)
{
  //si monitora il tempo di esecuzione
  TStopwatch timer;
  timer.Start();

  //Si creano file e istogramma per gli scarti azimutali in phi
  TFile file_scarti("file_scarti.root","RECREATE");
  TH1F *hist_scarti=new TH1F("hist_scarti","istogramma scarti tra le phi",200.,-0.02,0.02);

  
  //TTREE________________________________________________________________________________
  // Apertura del file di output per il TTree
  TFile hfile("htree.root","RECREATE"); 

  // dichiarazione del TTree
  TTree *tree = new TTree("T","TTree con 3 branches");

  // Dichiarazione dei TClonesArray di Hit che ci serviranno per il secondo e terzo branch
  TClonesArray *ptrhits1 = new TClonesArray("Hit",300);
  TClonesArray &hits1 = *ptrhits1;
  TClonesArray *ptrhits2 = new TClonesArray("Hit",300);
  TClonesArray &hits2 = *ptrhits2;

  // Definizione di una struct POINT che ci servira' per il primo branch
  typedef struct {
     Double_t X,Y,Z;
     Int_t mult;} POINT;      
  static POINT point; 

  // Dichiarazione dei 3 branch del TTree
  tree->Branch("VertMult",&point.X,"X/D:Y:Z:mult/I");
  tree->Branch("Hits1",&ptrhits1);
  tree->Branch("Hits2",&ptrhits2);

  //_________________________________________________________________________________________

  //Oggetti Hit: contengono Z, Phi ed etichetta degli hit
  Hit h1;
  Hit h2;
  //Oggetti Bool per validare le intersezioni con i layer (kFALSE se l'intersezione non avviene)
  Bool_t controlloL1;
  Bool_t controlloL2;
  //Etichette particelle
  UInt_t E1=0;
  UInt_t E2=0;  
  //variabile di loop
  UInt_t j=0;
  //Int_t k=0; per debug
   
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // FOR LOOP sul numero di eventi di collisione
  for(Int_t i=0; i<1000000;i++)
  { 
    E1=0;
    E2=0;
    //si generano vertice e molteplicita'
    Punto P0=pippo.EstraiVertice();
    UInt_t molteplicita=pippo.EstraiMolteplicita();
   

    //PRIMO BRANCH
    //assegnazione delle variabili per il primo branch del TTree (per ogni evento ci sarà un P0 e una molteplicità)  
    point.mult=molteplicita;
    point.X=P0.GetX();
    point.Y=P0.GetY();
    point.Z=P0.GetZ();
 
    //SECONDO e TERZO BRANCH
    //variabili per lo smearing
    Double_t Z1_smearing, Phi1_smearing, Z2_smearing, Phi2_smearing;
    Double_t SigmaZ_smearing=0.012; //![cm] in direzione z
    Double_t Delta=0.003; //![cm] in direzione r*Phi
    Double_t SigmaPhi1_smearing=Delta/R1; 
    Double_t SigmaPhi2_smearing=Delta/R2;

    //loop sulla molteplicita'
    for (j=0; j<molteplicita; j++)
      {	
	//funzione per il trasporto delle particelle	
	pippo.contohit(P0,h1,h2,controlloL1,controlloL2, ScatteringMultiploAttivo);

	//smearing dei punti di impatto
	Z1_smearing=h1.GetZ()+gRandom->Gaus(0.,SigmaZ_smearing);
	Phi1_smearing=h1.GetPhi()+gRandom->Gaus(0.,SigmaPhi1_smearing);
	Z2_smearing=h2.GetZ()+gRandom->Gaus(0.,SigmaZ_smearing);
        Phi2_smearing=h2.GetPhi()+gRandom->Gaus(0.,SigmaPhi2_smearing);

	//___________________________________________________________________
	//il TTree viene riempito solo se le particelle hanno generato un hit
	//hit su layer1
	if(controlloL1==kTRUE)
	  {  new(hits1[E1]) Hit(Z1_smearing,Phi1_smearing,E1); 
	     E1++;  }
	//hit su layer2	
	if(controlloL2==kTRUE)
	  {  new(hits2[E2]) Hit(Z2_smearing,Phi2_smearing,E2);
	     E2++;   }
	//____________________________________________________________________
 

	//scarti phi1-phi2
	if(controlloL1==kTRUE && controlloL2==kTRUE)
	  {
	    hist_scarti->Fill(Phi1_smearing-Phi2_smearing);
	  }

       }//fine del for molteplicita'

    //RUMORE
    //loop  sul rumore, separato per L1 e L2 
    for (j=0; j<NumeroPuntiRumore; j++)
      { new(hits1[E1]) Hit(-LunghezzaRiv/2.+LunghezzaRiv*gRandom->Rndm(),2.*pi*gRandom->Rndm(),-1);
	E1++;
      }
    for (j=0; j<NumeroPuntiRumore; j++)
      { new(hits2[E2]) Hit(-LunghezzaRiv/2.+LunghezzaRiv*gRandom->Rndm(),2.*pi*gRandom->Rndm(),-1);
	E2++;
      }

     /*
     // Debug
     /////////////////////////////////////////////////////////////////
     printf("Evento %d - moltepl: %d\n",i,molteplicita);
     printf("x= %f ; y= %f; z= %f \n",point.X,point.Y,point.Z);
     printf("Entries nel 1 TClonesArray: %d\n",ptrhits1->GetEntries());
     printf("Entries nel 2 TClonesArray: %d\n",ptrhits2->GetEntries());
     
     for (k=0; k<ptrhits1->GetEntries(); k++)
       {
    	Hit *tst1=(Hit*)ptrhits1->At(k);
    	cout<<"Hit "<<k<<") Z1, Phi1, etichetta = "<<tst1->GetZ()<<setw(10)<<tst1->GetPhi()<<setw(5)<<tst1->GetEtichetta()<<endl;
       }
     for (k=0; k<ptrhits2->GetEntries(); k++)
       {
	Hit *tst2=(Hit*)ptrhits2->At(k);
	cout<<"Hit "<<k<<") Z2, Phi2, etichetta = "<<tst2->GetZ()<<setw(10)<<tst2->GetPhi()<<setw(5)<<tst2->GetEtichetta()<<endl;
       }
      //fine Debug
      */
     

      tree->Fill();
      ptrhits1->Clear();
      ptrhits2->Clear();

  }//Fine loop eventi
  //___________________________________________________________________________________________________
  //___________________________________________________________________________________________________


  // Si salvano tutti gli oggetti del TTree nel file hfile  
  hfile.Write();
  hfile.Close();
  // Si cambia file per realizzare istogramma scarti
  file_scarti.cd();
  hist_scarti->Write(); 
  file_scarti.Close();

  //Stop timer
  timer.Stop();
  timer.Print();
  

}//fine funzione GeneraTree










