#include "/sphenix/user/dloomis/Tools/sPhenixStyle.C"

void drawLaserClusters()
{


    gROOT->SetBatch(kTRUE);
    SetsPhenixStyle();
    TFile *f75073 = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00075073.root");
    TFile *f78242 = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00078242_zerofield.root");

    TFile *f75073radialdriftcomponent = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00075073_radialdriftcomponent.root");
    TFile *f78242radialdriftcomponent = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00078242_zerofield_radialdriftcomponent.root");


    TGraph *g75073negz = (TGraph*)f75073->Get("gGarfieldStripePositions_xy_negz");
    TGraph *g78242negz = (TGraph*)f78242->Get("gGarfieldStripePositions_xy_negz");
    TGraph *g75073posz = (TGraph*)f75073->Get("gGarfieldStripePositions_xy_posz");
    TGraph *g78242posz = (TGraph*)f78242->Get("gGarfieldStripePositions_xy_posz");

    TGraph *g75073negz_radialdriftcomponent = (TGraph*)f75073radialdriftcomponent->Get("gGarfieldStripePositions_xy_negz");
    TGraph *g78242negz_radialdriftcomponent = (TGraph*)f78242radialdriftcomponent->Get("gGarfieldStripePositions_xy_negz");
    TGraph *g75073posz_radialdriftcomponent = (TGraph*)f75073radialdriftcomponent->Get("gGarfieldStripePositions_xy_posz");
    TGraph *g78242posz_radialdriftcomponent = (TGraph*)f78242radialdriftcomponent->Get("gGarfieldStripePositions_xy_posz");


    TCanvas *c1 = new TCanvas("c1", "Laser Cluster Positions", 1200, 600);
    c1->Divide(2, 1);
    c1->cd(1);
    g75073negz->SetTitle("");
    g75073negz->SetMarkerStyle(20);
    g75073negz->SetMarkerSize(0.5);
    g75073negz->SetMarkerColor(kBlue);
    g78242negz->SetMarkerStyle(20);
    g78242negz->SetMarkerSize(0.5);
    g78242negz->SetMarkerColor(kRed);
    g75073negz->Draw("AP");
    g78242negz->Draw("P SAME");
    g75073negz->GetXaxis()->SetTitle("x [cm]");
    g75073negz->GetYaxis()->SetTitle("y [cm]");

    c1->cd(2);
    g75073posz->SetTitle("");
    g75073posz->SetMarkerStyle(20);
    g75073posz->SetMarkerSize(0.5);
    g75073posz->SetMarkerColor(kBlue);
    g78242posz->SetMarkerStyle(20);
    g78242posz->SetMarkerSize(0.5);
    g78242posz->SetMarkerColor(kRed);
    g75073posz->Draw("AP");
    g78242posz->Draw("P SAME");
    g75073posz->GetXaxis()->SetTitle("x [cm]");
    g75073posz->GetYaxis()->SetTitle("y [cm]");

    TCanvas *c2 = new TCanvas("c2", "Laser Cluster Radial Drift Component", 1200, 600);
    c2->Divide(2, 1);
    c2->cd(1);
    g75073negz_radialdriftcomponent->SetTitle("");
    g75073negz_radialdriftcomponent->SetMarkerStyle(20);
    g75073negz_radialdriftcomponent->SetMarkerSize(0.5);
    g75073negz_radialdriftcomponent->SetMarkerColor(kBlue);
    g78242negz_radialdriftcomponent->SetMarkerStyle(20);
    g78242negz_radialdriftcomponent->SetMarkerSize(0.5);
    g78242negz_radialdriftcomponent->SetMarkerColor(kRed);
    g75073negz_radialdriftcomponent->Draw("AP");
    g78242negz_radialdriftcomponent->Draw("P SAME");
    g75073negz_radialdriftcomponent->GetXaxis()->SetTitle("x [cm]");
    g75073negz_radialdriftcomponent->GetYaxis()->SetTitle("y [cm]");

    c2->cd(2);
    g75073posz_radialdriftcomponent->SetTitle("");
    g75073posz_radialdriftcomponent->SetMarkerStyle(20);
    g75073posz_radialdriftcomponent->SetMarkerSize(0.5);
    g75073posz_radialdriftcomponent->SetMarkerColor(kBlue);
    g78242posz_radialdriftcomponent->SetMarkerStyle(20);
    g78242posz_radialdriftcomponent->SetMarkerSize(0.5);
    g78242posz_radialdriftcomponent->SetMarkerColor(kRed);
    g75073posz_radialdriftcomponent->Draw("AP");
    g78242posz_radialdriftcomponent->Draw("P SAME");
    g75073posz_radialdriftcomponent->GetXaxis()->SetTitle("x [cm]");
    g75073posz_radialdriftcomponent->GetYaxis()->SetTitle("y [cm]");


    c1->SaveAs("LaserClusterPositions.png");
    c2->SaveAs("LaserClusterRadialDriftComponent.png");


    


}