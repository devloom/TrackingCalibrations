void plotStripeRadialPositions()
{
    gROOT->SetBatch(kTRUE);

    const int nRadii{8};
    const int nStripes[4]{6,6,8,12};
    const int nPads[4]{96,96,128,192};
    const double RValues[4][8] = {
        {22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.63283566},
        {31.7648082, 32.89678074, 34.02875328, 35.16072582, 36.29269836, 37.4246709, 38.55664344, 39.68861597},
        {42.1705532, 44.2119258, 46.2532984, 48.29467608, 50.336069, 52.3774416, 54.4188015, 56.4601868},
        {59.46048725, 61.6545823, 63.84867738, 66.04277246, 68.23686754, 70.43096262, 72.6250577, 74.81915277}
    };

    const double RValuesFieldOnBeamOffSouth[4][8] = {
        {22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.4878},
        {31.6291, 32.7735, 33.9054, 35.0515, 36.177, 37.3157, 38.2845, 39.68861597},
        {42.0764, 44.1097, 46.1379, 48.1788, 50.2374, 52.3009, 54.3209, 56.4601868},
        {59.3461, 61.5651, 63.7675, 65.9972, 68.2256, 70.407, 72.5807, 74.81915277}
    };
    const double RValuesFieldOnBeamOffNorth[4][8] = {
        {22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.6064},
        {31.763, 32.8551, 34.036, 35.1597, 36.2798, 37.4069, 38.4812, 39.68861597},
        {42.0888, 44.1621, 46.2614, 48.3438, 50.4445, 52.487, 54.5199, 56.4601868},
        {59.5431, 61.7629, 63.946, 66.1569, 68.3603, 70.596, 72.7802, 74.81915277}
    };

    const double RValuesFieldOffBeamOffSouth[4][8] = {
        {22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.4878},
        {31.77, 32.5227, 33.706, 34.8955, 36.0802, 37.2196, 38.3897, 39.5012},
        {42.2528, 44.2081, 46.3044, 48.3487, 50.4203, 52.4585, 54.5596, 56.5377},
        {59.7544, 61.9264, 64.101, 66.3012, 68.5286, 70.7237, 72.9088, 74.9394}
    };
    const double RValuesFieldOffBeamOffNorth[4][8] = {
        {22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 31.5793},
        {32.5079, 33.6058, 34.757, 35.8817, 36.9847, 38.191, 39.2024, 39.68861597},
        {42.8751, 44.8936, 46.9494, 48.9999, 51.0452, 53.0825, 55.1411, 56.95},
        {60.3371, 62.5349, 64.728, 66.9481, 69.1665, 71.3437, 73.5473, 74.81915277}
    };

    TFile *f78242 = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/78242_fieldoffbeamoff/idealtpccoordinates_updatedfitstripereference/Laminations_run3auau_new_nocdbtag_v001-00078242.root");
    TFile *f83319 = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/83319_fieldonbeamoff/idealtpccoordinates_hugoupdatedstaticdistortionmap/Laminations_run3pp_new_nocdbtag_v001-00083319.root");

    TH2 *h78242South = (TH2*)f78242->Get("hPetal_South");
    TH2 *h78242North = (TH2*)f78242->Get("hPetal_North");
    TH2 *h83319South = (TH2*)f83319->Get("hPetal_South");
    TH2 *h83319North = (TH2*)f83319->Get("hPetal_North");

    TH2 *rResidualSouth = new TH2D("TPC South","TPC South;R [cm];R_{83319} - R_{78242} [cm]",50,30,80,50,-1.5,1.5);
    rResidualSouth->GetYaxis()->SetTitleSize(0.05);
    rResidualSouth->GetYaxis()->SetTitleOffset(0.75);
    TH2 *rResidualNorth = new TH2D("TPC North","TPC North;R [cm];R_{83319} - R_{78242} [cm]",50,30,80,50,-1.5,1.5);
    rResidualNorth->GetYaxis()->SetTitleSize(0.05);
    rResidualNorth->GetYaxis()->SetTitleOffset(0.75);

    double window = 0.5;
    double minWeight = 1.0;

    // ---------- SOUTH ----------
    for (int region = 0; region < 4; region++)
    {
      for (int ir = 0; ir < nRadii; ir++)
      {
        double truthR = RValuesFieldOffBeamOffSouth[region][ir];

        for (int ix = 1; ix <= h78242South->GetNbinsX(); ix++)
        {
          int ybinmin1 = h78242South->GetYaxis()->FindBin(truthR - window);
          int ybinmax1 = h78242South->GetYaxis()->FindBin(truthR + window);

          int ybinmin2 = h83319South->GetYaxis()->FindBin(truthR - window);
          int ybinmax2 = h83319South->GetYaxis()->FindBin(truthR + window);

          if (ybinmin1 < 1) ybinmin1 = 1;
          if (ybinmax1 > h78242South->GetNbinsY()) ybinmax1 = h78242South->GetNbinsY();

          if (ybinmin2 < 1) ybinmin2 = 1;
          if (ybinmax2 > h83319South->GetNbinsY()) ybinmax2 = h83319South->GetNbinsY();

          double sumW1 = 0.0, sumRW1 = 0.0;
          double sumW2 = 0.0, sumRW2 = 0.0;

          for (int iy = ybinmin1; iy <= ybinmax1; iy++)
          {
            double R = h78242South->GetYaxis()->GetBinCenter(iy);
            double w = h78242South->GetBinContent(ix, iy);
            sumW1 += w;
            sumRW1 += w * R;
          }

          for (int iy = ybinmin2; iy <= ybinmax2; iy++)
          {
            double R = h83319South->GetYaxis()->GetBinCenter(iy);
            double w = h83319South->GetBinContent(ix, iy);
            sumW2 += w;
            sumRW2 += w * R;
          }

          if (sumW1 < minWeight || sumW2 < minWeight) continue;

          double meanR1 = sumRW1 / sumW1;
          double meanR2 = sumRW2 / sumW2;
          double dR = meanR2 - meanR1;

          rResidualSouth->Fill(truthR, dR);
          // alternatively:
          // rResidualSouth->Fill(meanR1, dR);
        }
      }
    }

    // ---------- NORTH ----------
    for (int region = 0; region < 4; region++)
    {
      for (int ir = 0; ir < nRadii; ir++)
      {
        double truthR1 = RValuesFieldOffBeamOffNorth[region][ir];
        double truthR2 = RValuesFieldOnBeamOffNorth[region][ir];

        for (int ix = 1; ix <= h78242North->GetNbinsX(); ix++)
        {
          int ybinmin1 = h78242North->GetYaxis()->FindBin(truthR1 - window);
          int ybinmax1 = h78242North->GetYaxis()->FindBin(truthR1 + window);

          int ybinmin2 = h83319North->GetYaxis()->FindBin(truthR2 - window);
          int ybinmax2 = h83319North->GetYaxis()->FindBin(truthR2 + window);

          if (ybinmin1 < 1) ybinmin1 = 1;
          if (ybinmax1 > h78242North->GetNbinsY()) ybinmax1 = h78242North->GetNbinsY();

          if (ybinmin2 < 1) ybinmin2 = 1;
          if (ybinmax2 > h83319North->GetNbinsY()) ybinmax2 = h83319North->GetNbinsY();

          double sumW1 = 0.0, sumRW1 = 0.0;
          double sumW2 = 0.0, sumRW2 = 0.0;

          for (int iy = ybinmin1; iy <= ybinmax1; iy++)
          {
            double R = h78242North->GetYaxis()->GetBinCenter(iy);
            double w = h78242North->GetBinContent(ix, iy);
            sumW1 += w;
            sumRW1 += w * R;
          }

          for (int iy = ybinmin2; iy <= ybinmax2; iy++)
          {
            double R = h83319North->GetYaxis()->GetBinCenter(iy);
            double w = h83319North->GetBinContent(ix, iy);
            sumW2 += w;
            sumRW2 += w * R;
          }

          if (sumW1 < minWeight || sumW2 < minWeight) continue;

          double meanR1 = sumRW1 / sumW1;
          double meanR2 = sumRW2 / sumW2;
          double dR = meanR2 - meanR1;

          rResidualNorth->Fill(truthR1, dR);
          // or rResidualNorth->Fill(meanR1, dR);
        }
      }
    }

    // Profile overlays
    TProfile *pSouth = rResidualSouth->ProfileX("pSouth");
    TProfile *pNorth = rResidualNorth->ProfileX("pNorth");

    pSouth->SetDirectory(nullptr);
    pNorth->SetDirectory(nullptr);

    pSouth->SetMarkerStyle(20);
    pSouth->SetMarkerSize(0.9);
    pSouth->SetMarkerColor(kRed);
    pSouth->SetLineColor(kRed);
    pSouth->SetLineWidth(2);

    pNorth->SetMarkerStyle(20);
    pNorth->SetMarkerSize(0.9);
    pNorth->SetMarkerColor(kRed);
    pNorth->SetLineColor(kRed);
    pNorth->SetLineWidth(2);

    TCanvas *c1 = new TCanvas("c1","c1",700,1200);
    c1->Divide(1,2);

    c1->cd(1);
    gPad->SetRightMargin(0.12);
    rResidualSouth->Draw("COLZ");
    pSouth->Draw("SAME");

    c1->cd(2);
    gPad->SetRightMargin(0.12);
    rResidualNorth->Draw("COLZ");
    pNorth->Draw("SAME");

    c1->SaveAs("stripeRadialResiduals.pdf");

    delete pSouth;
    delete pNorth;
    delete c1;
    delete rResidualSouth;
    delete rResidualNorth;
    f78242->Close();
    f83319->Close();
    delete f78242;
    delete f83319;
}





// void plotStripeRadialPositions()
// {
//     gROOT->SetBatch(kTRUE);

//     const int nRadii{8};
//     const int nStripes[4]{6,6,8,12};
//     const int nPads[4]{96,96,128,192};
//     const double RValues[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.63283566},{31.7648082, 32.89678074, 34.02875328, 35.16072582, 36.29269836, 37.4246709, 38.55664344, 39.68861597},{42.1705532, 44.2119258, 46.2532984, 48.29467608, 50.336069, 52.3774416, 54.4188015, 56.4601868},{59.46048725, 61.6545823, 63.84867738, 66.04277246, 68.23686754, 70.43096262, 72.6250577, 74.81915277}};
    

//     const double RValuesFieldOnBeamOffSouth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.4878},{31.6291, 32.7735, 33.9054, 35.0515, 36.177, 37.3157, 38.2845, 39.68861597},{42.0764, 44.1097, 46.1379, 48.1788, 50.2374, 52.3009, 54.3209, 56.4601868},{59.3461, 61.5651, 63.7675, 65.9972, 68.2256, 70.407, 72.5807, 74.81915277}};
//     const double RValuesFieldOnBeamOffNorth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.6064},{31.763, 32.8551, 34.036, 35.1597, 36.2798, 37.4069, 38.4812, 39.68861597},{42.0888, 44.1621, 46.2614, 48.3438, 50.4445, 52.487, 54.5199, 56.4601868},{59.5431, 61.7629, 63.946, 66.1569, 68.3603, 70.596, 72.7802, 74.81915277}};

//     const double RValuesFieldOffBeamOffSouth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.4878},{31.77, 32.5227, 33.706, 34.8955, 36.0802, 37.2196, 38.3897, 39.5012},{42.2528, 44.2081, 46.3044, 48.3487, 50.4203, 52.4585, 54.5596, 56.5377},{59.7544, 61.9264, 64.101, 66.3012, 68.5286, 70.7237, 72.9088, 74.9394}};
//     const double RValuesFieldOffBeamOffNorth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 31.5793},{32.5079, 33.6058, 34.757, 35.8817, 36.9847, 38.191, 39.2024, 39.68861597},{42.8751, 44.8936, 46.9494, 48.9999, 51.0452, 53.0825, 55.1411, 56.95},{60.3371, 62.5349, 64.728, 66.9481, 69.1665, 71.3437, 73.5473, 74.81915277}};


//     TFile *f78242 = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/78242_fieldoffbeamoff/idealtpccoordinates_updatedfitstripereference/Laminations_run3auau_new_nocdbtag_v001-00078242.root");
//     TFile *f83319 = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/83319_fieldonbeamoff/idealtpccoordinates_hugoupdatedstaticdistortionmap/Laminations_run3pp_new_nocdbtag_v001-00083319.root");

//     TH2 *h78242South = (TH2*)f78242->Get("hPetal_South");
//     TH2 *h78242North = (TH2*)f78242->Get("hPetal_North");
//     TH2 *h83319South = (TH2*)f83319->Get("hPetal_South");
//     TH2 *h83319North = (TH2*)f83319->Get("hPetal_North");


//     TH2 *rResidualSouth = new TH2D("rResidualSouth","rResidualSouth",100,20,80,100,-0.5,0.5);
//     TH2 *rResidualNorth = new TH2D("rResidualNorth","rResidualNorth",100,20,80,100,-0.5,0.5);

//     for(int i = 1; i <= h78242South->GetNbinsX()/80; i++)
//     {
//       std::cout << i << " " << h78242South->GetNbinsX() << std::endl;
//       for(int j = 1; j <= h78242South->GetNbinsY(); j++)
//       {
//         double R1 = h78242South->GetYaxis()->GetBinCenter(j);

//         for (int k = 1; k <= h83319South->GetNbinsX()/80; k++)
//         {
//           int lmin = h83319South->GetYaxis()->FindBin(R1 - 1.0);
//           int lmax = h83319South->GetYaxis()->FindBin(R1 + 1.0);

//           if (lmin < 1) lmin = 1;
//           if (lmax > h83319South->GetNbinsY()) lmax = h83319South->GetNbinsY();

//           for (int l = lmin; l <= lmax; l++)
//           {
//             double R2 = h83319South->GetYaxis()->GetBinCenter(l);
//             rResidualSouth->Fill(R1, R2 - R1);
//           }
//         }
//       }
//     }


//     TCanvas *c1 = new TCanvas("c1","c1",800,600);
//     c1->cd();
//     gStyle->SetOptStat(0);
//     rResidualSouth->SetTitle("South Stripe Radial Position Residuals;R [cm]; Y [cm]");
//     rResidualSouth->Draw("COLZ");
//     c1->SaveAs("SouthStripeRadialPositionResiduals_FieldOnBeamOff.pdf");

    
//     // // TGraph *gNorthGeometry = new TGraph();
//     // // TGraph *gSouthGeometry = new TGraph();
//     // TGraph *gNorthFieldOnBeamOff = new TGraph();
//     // TGraph *gSouthFieldOnBeamOff = new TGraph();
    

//     // for(int module=1; module<4; module++)
//     // {
//     //   for(int j=0; j<nRadii-1; j++)
//     //   {
//     //     // gNorthGeometry->SetPoint(gNorthGeometry->GetN(), RValues[module][j], 0.0);
//     //     // gSouthGeometry->SetPoint(gSouthGeometry->GetN(), RValues[module][j], 0.0);
//     //     gNorthFieldOnBeamOff->SetPoint(gNorthFieldOnBeamOff->GetN(), RValues[module][j], RValuesFieldOnBeamOffNorth[module][j]-RValuesFieldOffBeamOffNorth[module][j]);
//     //     gSouthFieldOnBeamOff->SetPoint(gSouthFieldOnBeamOff->GetN(), RValues[module][j], RValuesFieldOnBeamOffSouth[module][j]-RValuesFieldOffBeamOffSouth[module][j]);
//     //   }
//     // }


//     // TCanvas *c1 = new TCanvas("c1","c1",800,600);
//     // c1->cd();
//     // gStyle->SetOptStat(0);
//     // gNorthFieldOnBeamOff->SetTitle("North Stripe Radial Positions;R [cm]; #Delta R [cm]");
//     // gNorthFieldOnBeamOff->SetMarkerStyle(20);
//     // gNorthFieldOnBeamOff->SetMarkerColor(kBlue);
//     // gNorthFieldOnBeamOff->Draw("AP");

//     // // TLine *line = new TLine(30,30,80,80);
//     // // line->SetLineColor(kBlack);
//     // // line->SetLineStyle(2);
//     // // line->Draw("same");

//     // c1->SaveAs("NorthStripeRadialPositions_FieldOnBeamOff.pdf");

//     // TCanvas *c2 = new TCanvas("c2","c2",800,600);
//     // c2->cd();
//     // gStyle->SetOptStat(0);
//     // gSouthFieldOnBeamOff->SetTitle("South Stripe Radial Positions;R [cm]; #Delta R [cm]");
//     // gSouthFieldOnBeamOff->SetMarkerStyle(20);
//     // gSouthFieldOnBeamOff->SetMarkerColor(kRed);
//     // gSouthFieldOnBeamOff->Draw("AP");

//     // // TLine *line2 = new TLine(30,30,80,80);
//     // // line2->SetLineColor(kBlack);
//     // // line2->SetLineStyle(2);
//     // // line2->Draw("same");

//     // c2->SaveAs("SouthStripeRadialPositions_FieldOnBeamOff.pdf");



// }