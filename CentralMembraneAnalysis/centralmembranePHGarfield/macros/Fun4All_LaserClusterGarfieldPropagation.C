#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllServer.h>
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
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phgarfield/PHGarfield.h>
#include <phool/recoConsts.h>

#include <centralmembranephgarfield/TpcLaserClusterGarfieldPropagator.h>

#include <fstream>
#include <iostream>
#include <string>

#include <TString.h>
#include <TSystem.h>
#include <TVector3.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4detectors_io.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libtrack_io.so)
R__LOAD_LIBRARY(libtrackbase_historic_io.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libPHGarfield.so)
R__LOAD_LIBRARY(libcentralmembraneanalysisv2.so)
R__LOAD_LIBRARY(libcentralmembranephgarfield.so)



void Fun4All_LaserClusterGarfieldPropagation(
    const int nEvents = 1000,
    const std::string filelist = "/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/laserclusters_lists/run79507.list",
    const std::string outdir = "/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/",
    const std::string fieldOffBeamOffStripeFile = "/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/78242_fieldoffbeamoff/idealtpccoordinates_updatedfitstripereference/Laminations_run3auau_new_nocdbtag_v001-00078242.root")
{
  auto se = Fun4AllServer::instance();
  se->Verbosity(0);

  TRACKING::tpc_zero_supp = true;
  G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING = true;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  auto rc = recoConsts::instance();
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  std::ifstream inputList(filelist);
  std::string firstInputFile;
  std::getline(inputList, firstInputFile);
  if (firstInputFile.empty())
  {
    std::cout << "Fun4All_LaserClusterGarfieldPropagation - empty input file list: " << filelist << std::endl;
    gSystem->Exit(1);
  }

  const auto runSegment = Fun4AllUtils::GetRunSegment(firstInputFile);
  const int runnumber = runSegment.first;
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  auto hitsin = new Fun4AllDstInputManager("InputManager");
  hitsin->AddListFile(filelist);
  se->registerInputManager(hitsin);

  const std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  auto ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
  G4TPC::laser_adc_threshold = 100;

  TrackingInit();
  ACTSGEOM::ActsGeomInit();

  const std::string electricFieldMap = "garfield_fields/sphenix_rossegger_garfield_field.root";
  auto* phg = new PHGarfield("PHGarfield", electricFieldMap, 1.0, 1.0);
  TVector3 Northxyz(-0.001, -0.001, 1123.109);  // mm
  TVector3 Southxyz(-3.354, -0.673, -1137.382);  // mm
  TVector3 center = 0.5 * (Northxyz + Southxyz);
  center *= 0.1;  // mm to cm
  phg->MoveTpc(center.X(), center.Y(), center.Z());
  phg->RotateTpc(0, 0.001485, 0);
  phg->RotateTpc(0.000298, 0, 0);
  phg->SetCMVoltageDefault(380.0);
  phg->SetZeroField(false);
  se->registerSubsystem(phg);

  auto propagator = new TpcLaserClusterGarfieldPropagator;
  propagator->setPHGarfield(phg);
  propagator->setDeltaMapBinning(48, 48);
  propagator->setStripeComparisonBranchWeights(1.0, 1.0, 1.0);
  // const std::string outputfile = Form("%s/LaserClusterGarfieldPropagation-%08d_zerofieldv4.root", outdir.c_str(), runnumber);
  const std::string outputfile = Form("%s/LaserClusterGarfieldPropagation-%08d_testzerofieldfalse.root", outdir.c_str(), runnumber);
  propagator->setOutputfile(outputfile);
  propagator->setFieldOffBeamOffStripeFile(fieldOffBeamOffStripeFile);
  se->registerSubsystem(propagator);

  se->run(nEvents);
  se->End();
  se->PrintTimer();

  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
