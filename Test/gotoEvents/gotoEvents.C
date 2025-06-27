
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;

void gotoEvents(){

    clas12reader c12("/lustre19/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo",{0});
    
    auto good= c12.grabEvent(100003);
    if(good){
      std::cout<<"event : rows  "<<c12.event()->getRows()<<" runconfig : rows "<<c12.runconfig()->getRows()<<" number of particles "<< c12.getDetParticles().size()<<endl;
      
      cout<<" start time "<<c12.event()->getStartTime()<<" DAQ event number "<<c12.runconfig()->getEvent()<<endl;;
  }
 }   
    
  
