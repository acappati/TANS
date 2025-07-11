//compilazione veloce

#include "TObject.h"
#include "TString.h"
#include "TSystem.h"

void Compilaclasse(TString myopt="fast"){
  TString opt;

  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("Punto.cpp",opt.Data());
  gSystem->CompileMacro("Retta.cpp",opt.Data());
  gSystem->CompileMacro("Cilindro.cpp",opt.Data());
  gSystem->CompileMacro("MyRandom3.cpp",opt.Data());
  gSystem->CompileMacro("Hit.cpp",opt.Data());
}
