// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLG4WAVEFORM_C
#define MACRO_FUN4ALLG4WAVEFORM_C

#include <GlobalVariables.C>

#include <G4_CEmc_Spacal.C>
#include <G4_HcalIn_ref.C>
#include <G4_HcalOut_ref.C>
#include <G4_Input.C>
#include <G4_Production.C>
#include <G4_Global.C>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include<Calo_Calib.C>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllUtils.h>
//#include <fun4allutils/TimerStats.h>
#include <fun4all/Fun4AllServer.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4allutils/TimerStats.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libCaloWaveformSim.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allutils.so)

#include <caloana/CaloAna.h>
R__LOAD_LIBRARY(libcaloana.so)

void Fun4All_G4_Waveform(
    const int nevents = 1,
    const string &inputFile0 = "dst_calo_cluster.list",
    const string &inputFile1 = "g4hits.list",
    const string &outdir = ".",
    int iter = 2,
    const string &cdbtag = "MDC2_ana.418")

{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);  // set it to 1 if you want event printouts

  recoConsts *rc = recoConsts::instance();

  ifstream file(inputFile0);
  string first_file;
  getline(file, first_file);
  //===============
  // conditions DB flags
  //===============
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  CDBInterface::instance()->Verbosity(1);

  //--------------
  // Set up Input Manager
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_CALO_CLUSTER");
  Fun4AllInputManager *intruth = new Fun4AllDstInputManager("G4Hits");
  cout << "add listfiles to input manager" << endl;
  in->AddListFile(inputFile0,1);
  intruth->AddListFile(inputFile1,1);
  cout << "files added" << endl;
  se->registerInputManager(in);
  se->registerInputManager(intruth);
  cout << "input manager registered" << endl;

  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTHIST_iter_%s", filename.c_str());

  /*
  // re-cluster
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  ClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC2");
  se->registerSubsystem(ClusterBuilder);
  //*/
  //global vertex reco
  //Global_Reco();

  //--------------Calibrating EMCal
  Process_Calo_Calib();
  ///////////////////
  // analysis modules
  if (iter > 1)
  {
    CaloAna *ca = new CaloAna("calomodulename", OutFile);
    ca->set_timing_cut_width(16);
    ca->apply_vertex_cut(false);
    ca->set_vertex_cut(20.);
    se->registerSubsystem(ca);
    std::cout << "Subsystems registered" << std::endl;
    //
  }


  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;
  std::cout << "All done!" << std::endl;

  TFile *f_done_signal = new TFile("DONE.root", "recreate");
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);
}

#endif  // MACRO_FUN4ALLG4WAVEFORM_C
