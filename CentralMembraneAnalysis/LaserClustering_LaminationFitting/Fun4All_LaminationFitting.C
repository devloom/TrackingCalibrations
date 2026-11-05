/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

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

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

//#include <tpc/LaserEventIdentifier.h>

#include <tpccalib/TpcLaminationFitting.h>

#include <stdio.h>
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libtpccalib.so)
void Fun4All_LaminationFitting(
    const int nEvents = 10,
    const std::string filelist = "/sphenix/lustre01/sphnxpro/commissioning/slurp/tpcbeam/run_00041900_00042000/",
    const std::string outdir = "/sphenig/tg/tg01/jets/bkimelman/",
    const std::string prod = "ana441_2024p007",
    const std::string type = "run2pp")

{

  auto se = Fun4AllServer::instance();
  se->Verbosity(1);


  TRACKING::tpc_zero_supp = true;
  G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING = true;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  
  std::ifstream ifs(filelist);
  std::string filepath;

  auto rc = recoConsts::instance();

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  int i=0;
  int runnumber = 0;
  int segment = 0;
  while(std::getline(ifs,filepath))
  {
    std::cout << "Adding DST with filepath: " << filepath << std::endl;
    if(i==0)
    {
      std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(filepath);
      runnumber = runseg.first;
      segment = runseg.second;
      rc->set_IntFlag("RUNNUMBER",runnumber);
      rc->set_uint64Flag("TIMESTAMP",runnumber);
    }

    //std::string inputname = "InputManager" + std::to_string(i);
    //auto hitsin = new Fun4AllDstInputManager(inputname);
    //hitsin->fileopen(filepath);
    //hitsin->AddFile(filepath);
    //se->registerInputManager(hitsin);
    i++;
  }

  auto hitsin = new Fun4AllDstInputManager("InputManager");
  hitsin->AddListFile(filelist);
  se->registerInputManager(hitsin);


  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  //  G4TPC::tpc_drift_velocity_reco = (8.0 / 1000) * 107.0 / 105.0;
  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  // G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = false;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  G4TPC::static_correction_filename = "DistortionMaps_Static/static_only_inverted_10-new_hugoupdate.root";

  if(runnumber == 52077 || runnumber == 53098 || (runnumber >= 75103 && runnumber <= 75114) || (runnumber >= 75555 && runnumber <= 75574) || runnumber == 78242 || runnumber == 81000 || runnumber == 81001 || runnumber == 84235 || runnumber == 84236 || (runnumber >= 84274 && runnumber <= 84277) || (runnumber >= 84279 && runnumber <= 84282) || (runnumber >= 84284 && runnumber <= 84287) || (runnumber >= 84289 && runnumber <= 84292) || (runnumber >= 84294 && runnumber <= 84297) || (runnumber >= 84299 && runnumber <= 84302) || (runnumber >= 84304 && runnumber <= 84307) || (runnumber >= 84309 && runnumber <= 84312) || (runnumber >= 84314 && runnumber <= 84317) || (runnumber >= 84319 && runnumber <= 84322) || (runnumber >= 84324 && runnumber <= 84327) || (runnumber >= 84329 && runnumber <= 84283))
  {
    G4TPC::ENABLE_STATIC_CORRECTIONS = false;
    //if(runnumber != 52077) G4TPC::BFieldOff = true;
    G4TPC::BFieldOff = true;
  }

  if(runnumber >= 75073 && runnumber <= 75090)
  {
    G4TPC::ENABLE_STATIC_CORRECTIONS = true;
    G4TPC::BFieldOff = true;
  }    

  // static corrections but fit laminations with straight lines
  if(runnumber >= 83313 && runnumber <= 83319)
  {
    G4TPC::ENABLE_STATIC_CORRECTIONS = true;
    G4TPC::BFieldOff = true;
  } 
  
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS=false;


  G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
  //G4TPC::average_correction_filename = std::string(Form("/sphenix/tg/tg01/jets/bkimelman/BenProduction/Feb25_2025/Laminations_run2pp_ana466_2024p012_v001-%08d.root",runnumber));
  //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS=true;
  //G4TPC::average_correction_interpolate = false;
  
  //G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;

  //TRACKING::reco_tpc_maxtime_sample = 1023;

  G4TPC::laser_adc_threshold = 100;

  //set to false if you want AuAu or true for pp
  if(type.find("pp") != std::string::npos) TRACKING::pp_mode = true;
  else TRACKING::pp_mode = false;


  
  TrackingInit();


  G4TPC::LaminationOutputName = Form("%s/Laminations_%s_%s-%08d.root",outdir.c_str(),type.c_str(),prod.c_str(),runnumber);
  G4TPC::LaminationQAName = Form("%s/LaminationQA_%s_%s-%08d.pdf",outdir.c_str(),type.c_str(),prod.c_str(),runnumber);
  TPC_LaminationFitting();
  
  se->run(nEvents);
  //se->run(-1);
  se->End();
  se->PrintTimer();

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
