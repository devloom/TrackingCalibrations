#include "/sphenix/user/dloomis/Tools/sPhenixStyle.C"

void drawLaminations()
{
    gROOT->SetBatch(kTRUE);
    SetsPhenixStyle();
    TFile *f79507 = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00079507_strength3.root");
    TFile *f78242 = TFile::Open("/sphenix/user/dloomis/DistortionsAnalysis/centralmembranePHGarfield/macros/output/LaserClusterGarfieldPropagation-00078242_zerofield.root");


    TH1D *g79507negz = (TH1D*)f79507->Get("hLaserCluster_garfield_negz");
    TH1D *g78242negz = (TH1D*)f78242->Get("hLaserCluster_garfield_negz");
    TH1D *g79507posz = (TH1D*)f79507->Get("hLaserCluster_garfield_negz");
    TH1D *g78242posz = (TH1D*)f78242->Get("hLaserCluster_garfield_negz");


    TExec *pal1 = new TExec("pal1", "gStyle->SetPalette(kBird);");
    TExec *pal2 = new TExec("pal2", "gStyle->SetPalette(56);");
    // TExec *pal3 = new TExec("pal3", "gStyle->SetPalette(kBird);");

    TCanvas *c1 = new TCanvas("c1", "Laser Cluster Positions", 1200, 600);
    c1->Divide(6, 3);

    for (int i = 0; i < 18; i++)
    {
        c1->cd(i + 1);

        float philow = ((i * 6.28 / 18) + ((i + 1) * 6.28 / 18)) / 2 - 0.05;
        float phihigh = ((i * 6.28 / 18) + ((i + 1) * 6.28 / 18)) / 2 + 0.05;
        // float philow = ((i * 6.28 / 18) + ((i + 1) * 6.28 / 18)) / 2 - 0.05 + 0.2;
        // float phihigh = ((i * 6.28 / 18) + ((i + 1) * 6.28 / 18)) / 2 + 0.05 + 0.2;

        gStyle->SetPalette(56);
        auto frame = gPad->DrawFrame(philow, 20, phihigh, 80);
        gPad->SetLogz();
        
        frame->SetTitle(";#phi;z");

        auto h1 = (TH1D*)g79507posz->Clone(Form("g79507posz_%d", i));
        auto h2 = (TH1D*)g78242posz->Clone(Form("g78242posz_%d", i));

        h1->GetXaxis()->SetRangeUser(philow, phihigh);
        h2->GetXaxis()->SetRangeUser(philow, phihigh);

        pal2->Draw();
        h2->Draw("colz same");

        pal1->Draw();
        h1->Draw("colz same");
    }

    c1->SaveAs("laminations_withspacecharge.png");

}