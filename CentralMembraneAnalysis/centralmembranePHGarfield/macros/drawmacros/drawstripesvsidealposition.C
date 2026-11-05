#include <cdbobjects/CDBTTree.h>

void drawstripesvsidealposition()
{
    TFile *fstripes = new TFile("../output/LaserClusterGarfieldPropagation-00078242_zerofield.root");
    TFile *fstripeswradialdrift = new TFile("../output/LaserClusterGarfieldPropagation-00078242_zerofield_radialdriftcomponent.root");

    TFile *fieldonbeamoffstripeswgarfield = new TFile("../output/LaserClusterGarfieldPropagation-00075073.root");
    TFile *fieldonbeamoffstripeswgarfieldwradialdrift = new TFile("../output/LaserClusterGarfieldPropagation-00075073_radialdriftcomponent.root");

    
    TH2D *hFieldOffBeamOff_negz = (TH2D*)fstripes->Get("hFieldOffBeamOffStripes_phir_negz");
    TH2D *hFieldOffBeamOff_posz = (TH2D*)fstripes->Get("hFieldOffBeamOffStripes_phir_posz");

    TH2D *hFieldOnBeamOff_negz = (TH2D*)fieldonbeamoffstripeswgarfield->Get("hLaserCluster_garfield_negz");
    TH2D *hFieldOnBeamOff_posz = (TH2D*)fieldonbeamoffstripeswgarfield->Get("hLaserCluster_garfield_posz");

    TH2D *hFieldOffBeamOffWithRadialDrift_negz = (TH2D*)fstripeswradialdrift->Get("hLaserCluster_garfield_negz");
    TH2D *hFieldOffBeamOffWithRadialDrift_posz = (TH2D*)fstripeswradialdrift->Get("hLaserCluster_garfield_posz");
    hFieldOffBeamOffWithRadialDrift_negz->SetName("hFieldOffBeamOffWithRadialDrift_negz");
    hFieldOffBeamOffWithRadialDrift_posz->SetName("hFieldOffBeamOffWithRadialDrift_posz");

    TH2D *hFieldOnBeamOffWithRadialDrift_negz = (TH2D*)fieldonbeamoffstripeswgarfieldwradialdrift->Get("hLaserCluster_garfield_negz");
    TH2D *hFieldOnBeamOffWithRadialDrift_posz = (TH2D*)fieldonbeamoffstripeswgarfieldwradialdrift->Get("hLaserCluster_garfield_posz");
    hFieldOnBeamOffWithRadialDrift_negz->SetName("hFieldOnBeamOffWithRadialDrift_negz");
    hFieldOnBeamOffWithRadialDrift_posz->SetName("hFieldOnBeamOffWithRadialDrift_posz");


    std::string sideal = "../../../../Distortions/StripePatterns/CMStripePattern_ideal.root";

    // Read ideal stripe positions
    CDBTTree *cdbttree = new CDBTTree(sideal);
    cdbttree->LoadCalibrations();
    auto cdbMap = cdbttree->GetDoubleEntryMap();
    
    // TH2D *hPetal_idealStripes[2];
    // hPetal_idealStripes[0] = new TH2D("hPetal_idealStripes_negz", "hPetal_idealStripes_negz", 360, 0, 2*M_PI, 100, 20.0, 80.0);
    // hPetal_idealStripes[1] = new TH2D("hPetal_idealStripes_posz", "hPetal_idealStripes_posz", 360, 0, 2*M_PI, 100, 20.0, 80.0);

    TGraph *gPetal_idealStripes[2];
    gPetal_idealStripes[0] = new TGraph();
    gPetal_idealStripes[0]->SetName("gPetal_idealStripes_negz");
    gPetal_idealStripes[1] = new TGraph();
    gPetal_idealStripes[1]->SetName("gPetal_idealStripes_posz");


    for (const auto &entry : cdbMap)
    {
      const auto index = entry.first;

      double phi, r;

      phi = cdbttree->GetDoubleValue(index, "truthPhi");
      r = cdbttree->GetDoubleValue(index, "truthR");

      if (phi < 0.0)
      {
        phi += 2 * M_PI;
      }
      if (index > 180000)
      {
        // hPetal_idealStripes[0]->Fill(cdbttree->GetDoubleValue(index, "truthPhi"), cdbttree->GetDoubleValue(index, "truthR"));
        gPetal_idealStripes[0]->SetPoint(gPetal_idealStripes[0]->GetN(), phi, r);
      }
      else
      {
        // hPetal_idealStripes[1]->Fill(cdbttree->GetDoubleValue(index, "truthPhi"), cdbttree->GetDoubleValue(index, "truthR"));
        gPetal_idealStripes[1]->SetPoint(gPetal_idealStripes[1]->GetN(), phi, r);
      }

    }

    TFile *fout = new TFile("stripesvsidealposition.root", "RECREATE");
    fout->cd();
    hFieldOffBeamOff_negz->Write();
    hFieldOffBeamOff_posz->Write();
    hFieldOffBeamOffWithRadialDrift_negz->Write();
    hFieldOffBeamOffWithRadialDrift_posz->Write();
    hFieldOnBeamOff_negz->Write();
    hFieldOnBeamOff_posz->Write();
    hFieldOnBeamOffWithRadialDrift_negz->Write();
    hFieldOnBeamOffWithRadialDrift_posz->Write();
    gPetal_idealStripes[0]->Write();
    gPetal_idealStripes[1]->Write();
    fout->Close();
}