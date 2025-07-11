#include <Riostream.h>
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Hit.h"
#include "CalcoliRicostruzione.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TCanvas.h"
#include <vector>
#include <utility>
#include <algorithm>
using namespace std;

ClassImp(CalcoliRicostruzione)


//__________________________________________________________________
CalcoliRicostruzione::CalcoliRicostruzione():TObject() 
  {
    // Default constructor
  }

//___________________________________________________________________
CalcoliRicostruzione::~CalcoliRicostruzione()
 {
    // Default destructor
 }




//_______________________________________________________
Double_t CalcoliRicostruzione::CalcoloFWHM()
{
   TFile F("file_scarti.root");
   TH1F *hist=dynamic_cast<TH1F*>(F.Get("hist_scarti")); 
   Double_t FWHMstima=4*hist->GetRMS();    
   F.Close();
   return FWHMstima;
}


//_____________________________________________________
//intersezione tra tracklet e asse del fascio, restituisce z candidato
Double_t CalcoliRicostruzione::CalcoloCandidatiVertice(Double_t z1, Double_t z2)
{
  Double_t ZVertice = (7.*z1 - 4.*z2) / 3.;
  return ZVertice;
}


//_____________________________________________________
void CalcoliRicostruzione::GeneraIstoCandidatiVertice(Double_t FWHM, TH1F *hist_vertice, vector <pair<Double_t,Double_t> > vL1, vector <pair<Double_t,Double_t> > vL2 , vector <pair<Double_t, Double_t> >::iterator it1, vector <pair<Double_t, Double_t> >::iterator it2, vector<Double_t> &vCandidatiVertice)
{
  Int_t contatore1=0;
  Int_t contatore2=0;
  Double_t DeltaPhi=0.;
  Double_t ZCandidato=0.;

    //loop sugli elementi di vL1
    for(it1=vL1.begin(); it1!=vL1.end(); it1++)
      {
	//loop sugli elementi di vL2
	for(it2=vL2.begin(); it2!=vL2.end(); it2++)
	  {
	    //si definisce DeltaPhi
	    DeltaPhi=TMath::Abs(it1->first - it2->first);	    
	   
	    //_______________________________________________________________
	    //COMPATIBILITA'
	    if(DeltaPhi<FWHM)
	      {
		contatore1++; //debug		      	
	
		//calcolo intersezioni tra tracklet e asse del fascio
		ZCandidato=CalcoloCandidatiVertice(it1->second,it2->second);
		//riempimento istogramma con i candidati vertici
		hist_vertice->Fill(ZCandidato);
		//riempimento vector vCandidatiVertice
		vCandidatiVertice.push_back(ZCandidato);
	      }
	    //NO COMPATIBILITA'
	    else
	      {		
		contatore2++;                
	      }
	    //________________________________________________________________


	  }//fine for su elementi di vL2
      }//fine for su elementi di vL1         
   //!hist_vertice->Write(); debug
}//fine funzione GeneraIstoCandidatiVertice



//________________________________________________________________________________
Bool_t CalcoliRicostruzione::GetZRicostruita_Hist(Double_t &z, TH1F *hist, Int_t &Count_NonRec, vector<Double_t> vCandidatiVertice)
{
  //ampiezza finestra "zoom"
  Double_t ampiezzafinestra=0.5; //5 mm

  //inizializzazione variabili
  Int_t contatore=0; 
  Double_t Somma=0.;
  Int_t ContatoreSomma=0;

  //si ordina il vector in ordine crescente
  std::sort (vCandidatiVertice.begin(), vCandidatiVertice.end());

  //informazioni su hist
  TAxis *xa=hist->GetXaxis(); 
  Int_t celle = hist->GetSize();
  Int_t bin_maxconteggi=hist->GetMaximumBin();  
  Int_t count_max=hist->GetBinContent(bin_maxconteggi);
  Double_t  z_maxconteggi=xa->GetBinCenter(bin_maxconteggi);  
  
 
  //ciclo FOR per correggere imprecisione di GetMaximumBin
  //si conta quante volte il massimo ha un numero di counts effettivamente maggiore rispetto ai counts delle altre celle
  for(Int_t i=0; i<celle;i++)
      {	       
	if( count_max - hist->GetBinContent(i) > 0)
	    { contatore++; }
      }

  

  //RICOSTRUZIONE
  //____________________________________________________________
  //per poter ricostruire ci si assicura che il massimo sia tale
  if(contatore == celle-1)
    {
      //FASE 1
      //si determinano estremi della finestra di ricostruzione
      Double_t A=z_maxconteggi-ampiezzafinestra/2.;
      Double_t B=z_maxconteggi+ampiezzafinestra/2.;
      
      for(UInt_t i=0; i<vCandidatiVertice.size();i++)
	{
	  //FASE 2
	  //solo i candidati presenti nella finestra contribuiscono al calcolo del vertice ricostruito
	  if(vCandidatiVertice[i]>A && vCandidatiVertice[i]<B)
	    {
	      Somma=Somma+vCandidatiVertice[i];
	      ContatoreSomma++;
	    }//fine if
	}//fine for

      //si ricostruisce la z come media dei valori presenti nella finestra
      z=Somma/(static_cast<Double_t>(ContatoreSomma));
           
      return kTRUE;
    }//fine if esterno
 
  else
    {
      Count_NonRec++;
      return kFALSE;
    }
  //fine algoritmo ricostruzione
  //_____________________________________________________________


}//fine funzione GetZRicostruita_Hist




//____________________________________________________________________
void CalcoliRicostruzione::GeneraIstoRisoluzione(TH1F *hist, vector<Double_t> vZVertice_Ricostruita, vector<Double_t>  vZVertice_MC)
{
   Double_t scarto;
   for(UInt_t i=0; i<vZVertice_Ricostruita.size();i++)
      {
	scarto=vZVertice_Ricostruita[i]-vZVertice_MC[i];
	hist->Fill(scarto);
       }  
   hist->Write();
}//fine funzione GeneraIstoRisoluzione







//////////////////////////////////
///////////  GRAFICI   ///////////
//////////////////////////////////

//____________________________________________________________________
void CalcoliRicostruzione::GraficoRisoluzioneVSMolteplicita(vector<Double_t> vZVertice_Ricostruita, vector<Double_t>  vZVertice_MC, vector<UInt_t> vMult)
{ 
  //variabili
  UInt_t M;
  Double_t VVero;
  Double_t VRec;
  //variabili per TGraph (numero punti n, ascisse x, ordinate y, errori)
  Int_t n=8;
  Double_t x[8]={1.5,4.5,8.,12.5,20.,30.,40.,50.};
  Double_t y[8];
  Double_t Ex[8];
  Double_t Ey[8];
  for(Int_t k=0; k<8;k++)
    { Ex[k]=0.;
      Ey[k]=0.; }  

  //RisoluzioneVSMolteplicita'
  //vettore che conterra' gli istogrammi
  TH1F *vRis_Mol[8];

  //creazione istogrammi  
  vRis_Mol[0]=new TH1F("Int0","Risoluzione con molteplicita in [0, 3]",200,-0.3,0.3);
  vRis_Mol[1]=new TH1F("Int1","Risoluzione con molteplicita in [3, 6]", 200,-0.3,0.3);
  vRis_Mol[2]=new TH1F("Int2","Risoluzione con molteplicita in [6, 10]",200,-0.3,0.3);
  vRis_Mol[3]=new TH1F("Int3","Risoluzione con molteplicita in [10, 15]",200,-0.3,0.3);
  vRis_Mol[4]=new TH1F("Int4","Risoluzione con molteplicita in [15, 25]",200,-0.3,0.3);
  vRis_Mol[5]=new TH1F("Int5","Risoluzione con molteplicita in [25, 35]",200,-0.3,0.3);
  vRis_Mol[6]=new TH1F("Int6","Risoluzione con molteplicita in [35, 45]",200,-0.3,0.3);
  vRis_Mol[7]=new TH1F("Int7","Risoluzione con molteplicita in [45, 58]",200,-0.3,0.3);
  
  //riempimento istogrammi
  for(UInt_t i=0; i<vMult.size();i++)
    {
      M=vMult[i];
      VVero=vZVertice_MC[i];
      VRec=vZVertice_Ricostruita[i];

      if(M<=3)           vRis_Mol[0]->Fill(VVero-VRec);
      if(M>3 && M<=6)    vRis_Mol[1]->Fill(VVero-VRec);
      if(M>6 && M<=10)   vRis_Mol[2]->Fill(VVero-VRec);
      if(M>10 && M<=15)  vRis_Mol[3]->Fill(VVero-VRec);
      if(M>15 && M<=25)  vRis_Mol[4]->Fill(VVero-VRec);
      if(M>25 && M<=35)  vRis_Mol[5]->Fill(VVero-VRec);
      if(M>35 && M<=45)  vRis_Mol[6]->Fill(VVero-VRec);
      if(M>45)           vRis_Mol[7]->Fill(VVero-VRec);    
          
    }
 
  //riempimento arrays per TGraph
  for(UInt_t i=0; i<8;i++)
    { 
      y[i]=vRis_Mol[i]->GetRMS();
      Ex[i]=0.;
      Ey[i]=vRis_Mol[i]->GetRMSError();          
      vRis_Mol[i]->Write();
    }

  //opzioni grafiche: Canvas e TGraph
  TCanvas *c1=new TCanvas("c1","Risoluzione VS Molteplicita'");     
  TGraphErrors *grafico1 = new TGraphErrors(n,x,y,Ex,Ey);

  grafico1->SetTitle("Risoluzione VS Molteplicita'");
  grafico1->GetXaxis()->SetTitle("Molteplicita'");
  grafico1->GetYaxis()->SetTitle("Risoluzione [cm]");
  grafico1->Draw("AL*");
  grafico1->Write("Risoluzione VS Molteplicita'");
  c1->SaveAs("c1.gif");

  //si libera la memoria
  for(UInt_t i=0; i<8;i++)
   { delete vRis_Mol[i];
     vRis_Mol[i]=NULL;
   }
}//fine funzione GraficoRisoluzioneVSMolteplicita


//_____________________________________________________________
void CalcoliRicostruzione::GraficoRisoluzioneVSZGenerata(vector<Double_t> vZVertice_Ricostruita, vector<Double_t>  vZVertice_MC)

{  
  //variabili strumentali
  Double_t Z;
  Double_t ZRec;
  //variabili per TGraph
  Int_t n=9;
  Double_t x[9]={-13.5,-9.,-4.5, -1.9, 0., 1.9, 4.5, 9.,13.5};
  Double_t ex[9];
  Double_t ey[9];
  Double_t y[9];

  //RisoluzioneVSZGenerata 
  //vettore che conterra' gli istogrammi
  TH1F *vRis_Z[9];

  //creazione istogrammi
  vRis_Z[0]=new TH1F("1","Risoluzione per Z generata tra -13.5 e -12.5 ",200,-0.3,0.3);
  vRis_Z[1]=new TH1F("2","Risoluzione per Z generata tra -12.5 e -6",200,-0.3,0.3);
  vRis_Z[2]=new TH1F("3","Risoluzione per Z generata tra -6 e -3 ",200,-0.3,0.3);
  vRis_Z[3]=new TH1F("4","Risoluzione per Z generata tra -3 e -0.7 ",200,-0.3,0.3);
  vRis_Z[4]=new TH1F("5","Risoluzione per Z generata tra -0.7 e 0.7 ",200,-0.3,0.3);
  vRis_Z[5]=new TH1F("6","Risoluzione per Z generata tra 0.7 e 3",200,-0.3,0.3);
  vRis_Z[6]=new TH1F("7","Risoluzione per Z generata tra 3 e 6",200,-0.3,0.3);
  vRis_Z[7]=new TH1F("8","Risoluzione per Z generata tra 6 e 12.5",200,-0.3,0.3);
  vRis_Z[8]=new TH1F("9","Risoluzione per Z generata tra 12.5 e 13.5",200,-0.3,0.3);
  

  //riempimento istogrammi
  for(UInt_t i=0; i<vZVertice_MC.size();i++)
    {
      Z=vZVertice_MC[i];
      ZRec=vZVertice_Ricostruita[i];

      if(Z<-12.5)              vRis_Z[0]->Fill(Z-ZRec);  
      if(Z>-12.5 && Z<=-6)     vRis_Z[1]->Fill(Z-ZRec);     
      if(Z>-6 && Z<=-3.)       vRis_Z[2]->Fill(Z-ZRec);
      if(Z>-3. && Z<=-0.7)     vRis_Z[3]->Fill(Z-ZRec);     
      if(Z>-0.7 && Z<=0.7)     vRis_Z[4]->Fill(Z-ZRec);  //centrale          
      if(Z>0.7 && Z<=3.)       vRis_Z[5]->Fill(Z-ZRec);
      if(Z>3. && Z<=6)         vRis_Z[6]->Fill(Z-ZRec);      
      if(Z>6 && Z<=12.5)       vRis_Z[7]->Fill(Z-ZRec);
      if(Z>12.5)               vRis_Z[8]->Fill(Z-ZRec);
    }
  
  //definizione variabili per TGraph
  for(UInt_t i=0; i<9;i++)
    {  
      y[i]=vRis_Z[i]->GetRMS(); 
      ex[i]=0.;
      ey[i]=vRis_Z[i]->GetRMSError();  
      vRis_Z[i]->Write();
    }

  //TGraph e Canvas
  TCanvas *c2=new TCanvas("c2","Risoluzione VS ZGenerata");
  TGraphErrors *grafico2 = new TGraphErrors(n,x,y,ex,ey);
  grafico2->SetTitle("Risoluzione VS ZGenerata");
  grafico2->GetXaxis()->SetTitle("Z Generata");
  grafico2->GetYaxis()->SetTitle("Risoluzione [cm]");
  grafico2->Draw("AL*");
  grafico2->Write("Risoluzione VS ZGenerata");
  c2->SaveAs("c2.gif");

  //si libera la memoria
  for(UInt_t i=0; i<9;i++)
    {
      delete vRis_Z[i];
      vRis_Z[i]=NULL;
    }

}//fine funzione GraficoRisoluzioneVSZGenerata

//___________________________________________________________________________
 void CalcoliRicostruzione::GraficoEfficienzaVSMolteplicita(vector <UInt_t> vMult_Tot, vector <UInt_t> vMult_Rec)
 {   
   //istogrammi
   TH1F *Numeratore= new TH1F("Numeratore","eventi ricostruiti",10, 0., 50.);
   TH1F *Denominatore= new TH1F("Denominatore","totali",10, 0., 50.);
   TH1F *Divisione= new TH1F("Divisione", "divisione", 10, 0., 50.);
   TAxis *xdivisione=Divisione->GetXaxis();
   //variabili
   Int_t CelleDivisione=10;    
   Double_t EventiCella=0;
   //arrays per TGraph
   Double_t Eff[10];
   Double_t X[10]; 
   Double_t XErr[10];
   Double_t YErr[10];
   for(Int_t i=0; i<10; i++)
	{
	  Eff[i]=0.;
	  X[i]=0.;
	  XErr[i]=0.;
	  YErr[i]=0.;
	}  

  
  //Riempimento Numeratore
  for(UInt_t j=0; j<vMult_Rec.size();j++)
    { Numeratore->Fill(vMult_Rec[j]); }
  //Riempimento Denominatore
  for(UInt_t i=0; i<vMult_Tot.size();i++)
    { Denominatore->Fill(vMult_Tot[i]); }
  //Funzione Divide per Divisione
   Divisione->Divide(Numeratore, Denominatore, 1., 1., "B");
   Divisione->Write();
         


  //Calcolo efficienza
  for(Int_t i=0; i<CelleDivisione;i++)
  { 
    EventiCella=Denominatore->GetBinContent(i+1);
    X[i]=xdivisione->GetBinCenter(i+1);
    Eff[i]=Divisione->GetBinContent(i+1);  
    YErr[i]=TMath::Sqrt((Eff[i]*(1.-Eff[i]))/EventiCella);    
  }

    
  //TGraph e TCanvas
  TCanvas *c3=new TCanvas("c3","Efficienza VS Molteplicita'");
  TGraphErrors *grafico3 = new TGraphErrors(10, X, Eff, XErr, YErr);
  grafico3->SetTitle("Efficienza VS Molteplicita'");
  grafico3->GetXaxis()->SetTitle("Molteplicita'");
  grafico3->GetYaxis()->SetTitle("Efficienza");
  grafico3->Draw("AL*");
  grafico3->Write("Efficienza VS Molteplicita'");
  c3->SaveAs("c3.gif");
  

}//fine funzione GraficoEfficienzaVSMolteplicita

 
