/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

#include <QA.C>
#include <fun4all/Fun4AllUtils.h>
#include <GlobalVariables.C>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <tpcqa/TpcLaserQA.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

//#include <tpc/LaserEventIdentifier.h>

#include <tpccalib/TpcCentralMembraneMatching.h>

#include <stdio.h>
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libtpccalib.so)
R__LOAD_LIBRARY(libtpcqa.so)
void Fun4All_LaserQA(
    const int nEvents = 10,
    const std::string tpcdir = "",
    const std::string tpcfilename = "/sphenix/lustre01/sphnxpro/commissioning/slurp/tpcbeam/run_00041900_00042000/",
    const std::string outdir = "/sphenig/tg/tg01/jets/bkimelman/")
{

  auto se = Fun4AllServer::instance();
  se->Verbosity(2);

  std::string inputFile = tpcdir + tpcfilename;


  
  std::pair<int, int>
    runseg = Fun4AllUtils::GetRunSegment(tpcfilename);
  int runnumber = runseg.first;
  int segment = runseg.second;
  
  auto rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_IntFlag("RUNNUMBER",runnumber);
  rc->set_uint64Flag("TIMESTAMP",runnumber);

  auto hitsin = new Fun4AllDstInputManager("InputManager");
  hitsin->fileopen(inputFile);
  se->registerInputManager(hitsin);
  

  TRACKING::tpc_zero_supp = true;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  //  G4TPC::tpc_drift_velocity_reco = (8.0 / 1000) * 107.0 / 105.0;
  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS=false;


  G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  G4TPC::average_correction_filename = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS=false;
  G4TPC::average_correction_interpolate = false;
  
  //G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;

  //TRACKING::reco_tpc_maxtime_sample = 1023;

  TrackingInit();


  auto LaserQA = new TpcLaserQA;
  se->registerSubsystem(LaserQA);
  
  
  se->run(nEvents);
  //se->run(-1);
  se->End();
  se->PrintTimer();
  
  TString qaname = Form("%s/HIST-%08d-%05d.root",outdir.c_str(),runnumber,segment);
  //TString qaname = "HIST_" + outfilename;
  std::string qaOutputFileName(qaname.Data());
  QAHistManagerDef::saveQARootFile(qaOutputFileName);
  
  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
