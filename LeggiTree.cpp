#include <Riostream.h>
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Punto.h"
#include "Retta.h" 
#include "Cilindro.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "MyRandom3.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "Hit.h"
#include "TH1F.h"
#include <vector>
#include <utility>
#include "CalcoliRicostruzione.h"
#include "TStopwatch.h"

using namespace std;

void LeggiTree()
 { 
   //Si monitora il tempo di esecuzione   
   TStopwatch timer;
   timer.Start();

   //DEFINIZIONE VARIABILI
   //variabili per riempire pair
   Double_t PhiL1;
   Double_t PhiL2;
   Double_t ZL1;
   Double_t ZL2;
   //variabili per ricostruzione del vertice
   Double_t Z_ric;
   Bool_t ricostruzione=kFALSE;
   Int_t Count_EventiNonRicostruiti=0;
   Int_t Count_EventiRicostruiti=0;
   //variabili per risultati finali
   Int_t NTot;
   Double_t Efficienza;
   Double_t Risoluzione;


   //Si definisce un oggetto di CalcoliRicostruzione
   CalcoliRicostruzione pluto=CalcoliRicostruzione();
   //Si calcola dall'istogramma realizzato con verita' MC il valore di deltaphi-di-accettanza
   Double_t FWHM=pluto.CalcoloFWHM();
   

   
   // DEFINIZIONE TTREE
   //_________________________________________________________

   //Definizione struct
   typedef struct {
      Double_t X,Y,Z;
      Int_t mult;} POINT;
   static POINT point;

  //Dichiarazione TClonesArray
  TClonesArray *hits1 = new TClonesArray("Hit",300);
  TClonesArray *hits2 = new TClonesArray("Hit",300);

  //Apertura file di input
  TFile hfile("htree.root");

  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1=tree->GetBranch("VertMult");
  TBranch *b2=tree->GetBranch("Hits1");
  TBranch *b3=tree->GetBranch("Hits2");

  //Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
  b2->SetAddress(&hits1);
  b3->SetAddress(&hits2);

  //___________________________________________________________________
 
  //Si utilizzeranno dei vector per copiare i dati contenuti nel TTree
  //Per ogni evento il vector sarà riempito con un pair <Phi, Z> in modo indipendente per il layer1 (L1) e per il layer2 (L2)

  //definizione degli iterator per lettura vector
  vector <pair<Double_t, Double_t> >::iterator it1;
  vector <pair<Double_t, Double_t> >::iterator it2;

  //vector per lettura info da TTree
  vector <pair<Double_t, Double_t> > vL1;
  vector <pair<Double_t, Double_t> > vL2;
  vector <UInt_t> vMult_Rec;
  vector <UInt_t> vMult_Tot;  
  //vector per ricostruzione
  vector <Double_t> vCandidatiVertice;
  vector <Double_t> vZVertice_MC;
  vector <Double_t> vZVertice_Ricostruita;
 

  //Creazione file "file_grafici" dove verranno stampati tutti i grafici
  TFile file_grafico("grafici.root","recreate");

  //Istogramma in cui raccogliere le intersezioni tra tracklets e asse del fascio
  TH1F *hist_vertice= new TH1F("hist_vertice","Candidati vertice",200.,-14.,14.);
  hist_vertice->GetXaxis()->SetTitle("Z [cm]");
  
  
  //_________________________________________________________
  // LETTURA TTREE
  //_________________________________________________________

  NTot=tree->GetEntries();

  //FOR LOOP sugli ingressi nel TTree
  for(Int_t ev=0; ev<NTot; ev++)
    {
      //Si svuotano l'istogramma e il vector con le z candidate vertice
      hist_vertice->Reset();
      vCandidatiVertice.clear();

      //Si seleziona l'evento da leggere
      tree->GetEvent(ev);

      //////////////////////////
      // lettura PRIMO BRANCH //
      //////////////////////////
      vMult_Tot.push_back(point.mult);
            

     
      ////////////////////////////
      // lettura SECONDO BRANCH //
      ////////////////////////////
      //Numero elementi nel TClonesArray1
      Int_t num1=hits1->GetEntries();
      
      //Loop per riempimento vector vL1
      //nei vector vL1 e vL2 vengono salvate solo i pair <Phi, Z> ma non le etichette (info utile per eventuale debug)
      for (Int_t j=0; j<num1; j++)
      {	Hit *tst1=(Hit*)hits1->At(j);	
	PhiL1=tst1->GetPhi();
	ZL1=tst1->GetZ();        
	vL1.push_back(make_pair(PhiL1,ZL1));       
      }

      //////////////////////////
      // lettura TERZO BRANCH //
      //////////////////////////
      //Numero elementi nel TClonesArray2
      Int_t num2=hits2->GetEntries();
      
      //Loop per riempimento vector vL1
      for (Int_t j=0; j<num2; j++)
	{ Hit *tst2=(Hit*)hits2->At(j);
	  PhiL2=tst2->GetPhi();
	  ZL2=tst2->GetZ();
	  vL2.push_back(make_pair(PhiL2,ZL2));
	}
        
    
      /*
      //debug per controllare che vL1 e vL2 siano stati riempiti correttamente
      //cout<<"Numero di elementi nel TClonesArray1: "<<num1<<" e nel TClonesArray2"<<num2<<endl;
      for(it2=vL2.begin(); it2!=vL2.end(); it2++)
      {	cout<<"Phi= "<<it2->first<<endl; 
        cout<<"Z= "<<it2->second<<endl; 
        cout<<"con FWHM giusta: "<<FWHM<<endl;
      }
      */

      //_______________________________________________________
      //_______________________________________________________
      //             R I C O S T R U Z I O N E
      //_______________________________________________________
      //_______________________________________________________
      //In questa fase si ricostruisce la coordinata z del vertice (criterio di "match": Phi di accettanza)
      //L'angolo di accettanza (chiamato FWHM) deriva dalla distribuzione Phi1-Phi2 della verita' MonteCarlo.

      
      //Si crea l'istogramma con i valori candidati di z del vertice
      pluto.GeneraIstoCandidatiVertice(FWHM, hist_vertice, vL1, vL2, it1, it2, vCandidatiVertice);
      
      //Calcolo della posizione del vertice come posizione piu' probabile di int delle tracklets
      // *se ricostruzione=kTRUE, ricostruzione avvenuta
      // *se ricostruzione=kFALSE, contatore Count_EventiNonRicostruiti++
      ricostruzione=pluto.GetZRicostruita_Hist(Z_ric, hist_vertice, Count_EventiNonRicostruiti, vCandidatiVertice);

      if(ricostruzione==kTRUE)  
	{ 
	  //Riempimento vector con coordinate z del vertice di collisione e quello con la verita' MC
	  vZVertice_MC.push_back(point.Z); 
	  vZVertice_Ricostruita.push_back(Z_ric); 
	  //Riempimento vector vMult
	  vMult_Rec.push_back(point.mult); 
	  
	} 
           
      //Si svuotano i vector vL1 e vL2
      vL1.clear();    
      vL2.clear();   

}//Fine for loop sugli ingressi del TTree

  //Si libera la memoria di hist_vertice
  delete  hist_vertice;
  hist_vertice=NULL;




  //_________________________________
  //GRAFICI
  //_________________________________

  pluto.GraficoRisoluzioneVSMolteplicita(vZVertice_Ricostruita, vZVertice_MC, vMult_Rec);
  pluto.GraficoRisoluzioneVSZGenerata(vZVertice_Ricostruita, vZVertice_MC);
  pluto.GraficoEfficienzaVSMolteplicita(vMult_Tot, vMult_Rec);
  
  
  //Istogramma per il calcolo della risoluzione (tutti gli eventi)
  TH1F *hist_risoluzione= new TH1F("hist_risoluzione","Risoluzione totale",200.,-0.3,0.3);
  hist_risoluzione->GetXaxis()->SetTitle("Zric - Zvero [cm]");  
  pluto.GeneraIstoRisoluzione(hist_risoluzione,vZVertice_Ricostruita, vZVertice_MC);    
  Risoluzione=hist_risoluzione->GetRMS();
   
  //Chiusura file che contiene tutti gli istogrammi
  file_grafico.Close();



  //__________________________________________________
  //RISULTATI FINALI  
  //_________________________________________________
  
  //Efficienza del programma di ricostruzione
  Count_EventiRicostruiti=NTot-Count_EventiNonRicostruiti;
  Efficienza=100.*(Count_EventiRicostruiti)/NTot;
  cout<<"***********************************************"<<endl;
  cout<<"Numero eventi simulati: "<<NTot<<endl;
  cout<<"Numero eventi ricostruiti: "<<Count_EventiRicostruiti<<endl;
  cout<<"Numero eventi non ricostruiti: "<<Count_EventiNonRicostruiti<<endl;
  cout<<"Efficienza di ricostruzione: "<<Efficienza<<" %"<<endl;
  cout<<"Risoluzione: "<<Risoluzione<<" cm "<<endl; 
  cout<<"***********************************************"<<endl;

  timer.Stop();
  cout<<"Tempo esecuzione LeggiTree"<<endl;
  timer.Print();
  
 
}//fine LeggiTree




