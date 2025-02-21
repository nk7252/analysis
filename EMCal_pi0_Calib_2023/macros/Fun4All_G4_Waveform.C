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
#include <G4_Mbd.C>
#include <G4_Global.C>
#include <g4mbd/MbdDigitization.h>
#include <g4mbd/MbdVertexFastSimReco.h>
#include <mbd/MbdReco.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>

#include <calowaveformsim/CaloWaveformSim.h>

#include<Calo_Calib.C>

//jet includes
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <G4_Jets.C>
#include <HIJetReco.C>

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
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libg4mbd.so)

#include <caloana/CaloAna.h>
R__LOAD_LIBRARY(libcaloana.so)

void Fun4All_G4_Waveform(
    const int nevents = 1,
    const string &inputFile0 = "/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/dst_calo_cluster.list",
    const string &inputFile1 = "/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/listfiles/single/run24/pi0/g4hits.list",
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
  //std::cout << "Timestamp : " << rc->get_uint64Flag("TIMESTAMP") << std::endl;
  CDBInterface::instance()->Verbosity(1);

  //--------------
  // Set up Input Manager
  cout << "Setting up input manager" << endl;
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_CALO_CLUSTER");
  Fun4AllInputManager *intruth = new Fun4AllDstInputManager("DST_TRUTH");


  cout << "add listfiles to input manager" << endl;
  in->AddListFile(inputFile0,1);
  intruth->AddListFile(inputFile1,1);
  cout << "files added" << endl;

  cout << "register input manager" << endl;
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
  Enable::MBDFAKE = true;// for single particle
  //Enable::MBDRECO = true; //for full pythia 
  Mbd_Reco();
  Global_Reco();
  //--------------Calibrating EMCal
  Process_Calo_Calib();
  ///////////////////
  //std::cout << "Running HIJetReco" << std::endl;
  GlobalVertex::VTXTYPE vertex_type = GlobalVertex::MBD;

  JetReco *truthjetreco = new JetReco();
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(1);  // (1) for pythia simulations, (2) for pythia embedding into hijing
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjetreco->set_algo_node("ANTIKT");
  truthjetreco->set_input_node("TRUTH");
  truthjetreco->Verbosity(0);
  se->registerSubsystem(truthjetreco);

  /*  
  //Jet_Reco();
  //HIJetReco();
  std::string jetreco_input_prefix = "TOWERINFO_CALIB";
  std::string jetreco_wvm_input_prefix = "WAVEFORM";
  
  
  JetReco *truthjetreco = new JetReco();
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjetreco->set_algo_node("ANTIKT");
  truthjetreco->set_input_node("TRUTH");
  truthjetreco->Verbosity(0);
  se->registerSubsystem(truthjetreco);
  RetowerCEMC *rcemc = new RetowerCEMC(); 
  rcemc->Verbosity(0); 
  rcemc->set_towerinfo(true);
  rcemc->set_frac_cut(0.5); //fraction of retower that must be masked to mask the full retower
  rcemc->set_towerNodePrefix(jetreco_input_prefix);
  se->registerSubsystem(rcemc);

  JetReco *towerjetrecounsub = new JetReco("TOWERJETRECO");
  TowerJetInput *incemc = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, jetreco_input_prefix);
  TowerJetInput *inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO, jetreco_input_prefix);
  TowerJetInput *inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO, jetreco_input_prefix);
  incemc->set_GlobalVertexType(HIJETS::vertex_type);
  inihcal->set_GlobalVertexType(HIJETS::vertex_type);
  inohcal->set_GlobalVertexType(HIJETS::vertex_type);
  towerjetrecounsub->add_input(incemc);//Jet::CEMC_TOWER
  towerjetrecounsub->add_input(inihcal);//Jet::HCALIN_TOWER
  towerjetrecounsub->add_input(inohcal);//Jet::HCALOUT_TOWER
  towerjetrecounsub->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_r04");
  towerjetrecounsub->set_algo_node("ANTIKT");
  towerjetrecounsub->set_input_node("TOWER");
  towerjetrecounsub->Verbosity(0);
  se->registerSubsystem(towerjetrecounsub);
  //*/

  // analysis modules
  if (iter > 1)
  {
    CaloAna *ca = new CaloAna("calomodulename", OutFile);
    //ca->set_timing_cut_width(16);
    ca->Set_Debug(0);//set to 1 for debug output

    ca->set_SPMC_bools(false, false, false);//SPMC, Pythia_weight, eta_weight
    ca->set_general_bools(false, false, false, false);// debug, Cluster_Debug, Cluster_Debug2, etabyeta
    ca->set_cut_bools(false, true, false, true);// eTCutbool, etaCutbool, clusterprobcut, zvtxcut_bool
    //ca->set_clusprob_cut(0.1);if false above no need to set this
    ca->set_cluschi2_cut(10);
    ca->set_eta_cut(0.6);
    ca->set_zvtx_cut(30.);
    ca->set_pythiajets(true);// set to true if you want to cut on fully efficient range for jets
    ca->set_EfficiencyRange(0, 3000);// efficiency range for jets. effectively off for now using 0, 3000
    //14,30 for jet 10 // 30,3000 for jet 30 // 0,14 for MB
    ca->set_cluspt_cut(0.6, 1.0);

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
