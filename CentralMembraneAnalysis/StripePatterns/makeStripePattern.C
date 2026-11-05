#include <cdbobjects/CDBTTree.h>

R__LOAD_LIBRARY(libcdbobjects.so);

void makeStripePattern()
{
  gROOT->SetBatch(kTRUE);

  const int nRadii{8};
  const int nStripes[4]{6,6,8,12};
  const int nPads[4]{96,96,128,192};
  const double RValues[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.63283566},{31.7648082, 32.89678074, 34.02875328, 35.16072582, 36.29269836, 37.4246709, 38.55664344, 39.68861597},{42.1705532, 44.2119258, 46.2532984, 48.29467608, 50.336069, 52.3774416, 54.4188015, 56.4601868},{59.46048725, 61.6545823, 63.84867738, 66.04277246, 68.23686754, 70.43096262, 72.6250577, 74.81915277}};

  const double RValuesFieldOnBeamOffSouth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.4878},{31.6291, 32.7735, 33.9054, 35.0515, 36.177, 37.3157, 38.2845, 39.68861597},{42.0764, 44.1097, 46.1379, 48.1788, 50.2374, 52.3009, 54.3209, 56.4601868},{59.3461, 61.5651, 63.7675, 65.9972, 68.2256, 70.407, 72.5807, 74.81915277}};
  const double RValuesFieldOnBeamOffNorth[4][8] = {{22.70902789, 23.84100043, 24.97297296, 26.1049455, 27.23691804, 28.36889058, 29.50086312, 30.6064},{31.763, 32.8551, 34.036, 35.1597, 36.2798, 37.4069, 38.4812, 39.68861597},{42.0888, 44.1621, 46.2614, 48.3438, 50.4445, 52.487, 54.5199, 56.4601868},{59.5431, 61.7629, 63.946, 66.1569, 68.3603, 70.596, 72.7802, 74.81915277}};

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

  const int keepThisAndAfter[8]{1,0,1,0,1,0,1,0};
  const int keepUntil[4][8]{{4,4,5,4,5,5,5,5},{5,5,6,5,6,5,6,5},{7,7,8,7,8,8,8,8},{11,10,11,11,11,11,12,11}};

  const double phi_petal = TMath::Pi()/6.0;
  const int pr_mult = 3;
  const int dw_mult = 8;
  const double diffwidth = 0.06;
  const double adjust = 0.015;

  std::vector<double> truthR;
  std::vector<double> truthPhi;
  std::vector<double> truthPhiOriginal;
  std::vector<int> truthIndex;

  vector<TVector3> pos[4][nRadii];
  vector<TVector3> posFieldOnBeamOffSouth[4][nRadii];
  vector<TVector3> posFieldOnBeamOffNorth[4][nRadii];

  // TFile *flaminationreferencepositions = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationreferencepositions/fieldOffPhiOffsets_GPR_updated.root");
  // TFile *flaminationreferencepositions = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationreferencepositions/fieldOffPhiOffsets_GPR_75074.root");
  TFile *flaminationreferencepositions = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationreferencepositions/fieldOffPhiOffsets_GPR_78242.root"); // field off beam off in ideal coordinates (78242)
  // TFile *flaminationreferencepositions = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationreferencepositions/fieldOnPhiOffsets_GPR_83319.root"); // field on beam off with static distortion corrections in ideal coordinates (83319)
  // TFile *flaminationreferencepositions = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationreferencepositions/fieldOnPhiOffsets_GPR_83313_alignedcoordinates.root"); // field on beam off with static distortion corrections in aligned coordinates (83313)
  TGraph *gPhiOffsetsSouth = (TGraph*)flaminationreferencepositions->Get("gOffsetsSouth_GPmean");
  TGraph *gPhiOffsetsNorth = (TGraph*)flaminationreferencepositions->Get("gOffsetsNorth_GPmean");

  for(int module=0; module<4; module++)
  {
    double spacing[nRadii];
    for(int j=0; j<nRadii; j++)
    {
      spacing[j] = 2.0 * ((dw_mult * diffwidth / RValues[module][j]) + (pr_mult * phi_petal / nPads[module]));
    }
    double phiTmp = 0.0;
    for(int j=0; j<nRadii; j++)
    {
      for(int k=keepThisAndAfter[j]; k<keepUntil[module][j]; k++)
      {
        if(j % 2 == 0)
        {
          phiTmp = k * spacing[j] + (spacing[j] / 2.0) - adjust;
        }
        else
        {
          phiTmp = (k + 1) * spacing[j] - adjust;
        }
        pos[module][j].emplace_back(RValues[module][j]*cos(phiTmp), RValues[module][j]*sin(phiTmp), 0.0);
      }
    }
  }

  

  for(int module=0; module<4; module++)
  {
    for(int j=0; j<nRadii; j++)
    {
      for(int i=0; i<(int)pos[module][j].size(); i++)
      {
        for(int k=0; k<18; k++)
        {
          if(k>0) pos[module][j][i].RotateZ(TMath::Pi() / 9);


          TVector3 NorthPos = pos[module][j][i];
          NorthPos.SetX(-NorthPos.X());
          //NorthPos.RotateZ(k * TMath::Pi() / 9);

          // truthR.push_back((NorthPos.Perp()+0.75)/(1.0 + 0.001)); // this is Ben's field off ideal coordinates 
          // truthR.push_back((NorthPos.Perp()+0.500)/(1.0 - 0.006)); // this is from my field off aligned coordinates
          // truthR.push_back((NorthPos.Perp()+0.420)/(1.0 - 0.02225)); // this is from my field on beam off aligned coordinates (without static corrections)
          // truthR.push_back((NorthPos.Perp()-0.10523)/(1.0 - 0.00260366)); // this is field on beam off ideal coordinate from 83319 (with static corrections)
          // truthR.push_back((NorthPos.Perp()-0.168367)/(1.0 - 0.00446341)); // this is field on beam off aligned coordinates from 83313 (with static corrections)
          truthR.push_back(RValuesFieldOffBeamOffNorth[module][j]); // this is with individual fits to each stripe row without being restricted to a linear fit across all R
          // truthR.push_back(RValuesFieldOnBeamOffNorth[module][j]); // this is with individual fits to each stripe row without being restricted to a linear fit across all R


          //truthR.push_back(NorthPos.Perp());
          double phiTmp = NorthPos.Phi();
          //double phi = phiTmp;
          // double phi = phiTmp - 0.00303345 + 0.0010828 * cos(phiTmp - 1.03718); // this is Ben's field off ideal coordinates
          double phi = phiTmp + gPhiOffsetsNorth->Eval(phiTmp); // this is from GPR

          if(phi < -TMath::Pi()) phi += 2*TMath::Pi();
          if(phi > TMath::Pi()) phi -= 2*TMath::Pi();
          truthPhi.push_back(phi);
          if(phiTmp < 0.0) phiTmp += 2*TMath::Pi();
          truthPhiOriginal.push_back(phiTmp);
          if(k<=8)
          {
            truthIndex.push_back((8-k)*10000 + (j+(module*nRadii))*100 + ((int)pos[module][j].size()-i-1));
            //truthIndex.push_back((26-k)*10000 + (j+(module*nRadii))*100 + ((int)pos[module][j].size()-i-1));
          }
          else
          {
            truthIndex.push_back((26-k)*10000 + (j+(module*nRadii))*100 + ((int)pos[module][j].size()-i-1));
            //truthIndex.push_back((44-k)*10000 + (j+(module*nRadii))*100 + ((int)pos[module][j].size()-i-1));
          }

          TVector3 SouthPos = pos[module][j][i];
          SouthPos.RotateZ(-TMath::Pi() / 18);
          //SouthPos.RotateZ(k * TMath::Pi() / 9);

          // truthR.push_back((SouthPos.Perp()-0.875)/(1.0 - 0.02)); // this is Ben's field off ideal coordinates 
          // truthR.push_back((SouthPos.Perp()-0.500)/(1.0 - 0.011)); // this is from my field off aligned coordinates
          // truthR.push_back((SouthPos.Perp()+0.670)/(1.0 - 0.02225)); // this is from my field on beam off aligned coordinates (without static corrections)
          // truthR.push_back((SouthPos.Perp()-0.315689)/(1.0 - 0.00409146)); // this is field on beam off ideal coordinates from 83319 (with static corrections)
          // truthR.push_back((SouthPos.Perp()-0.336735)/(1.0 - 0.00446341)); // this is field on beam off aligned coordinates from 83313 (with static corrections)
          truthR.push_back(RValuesFieldOffBeamOffSouth[module][j]); // this is with individual fits to each stripe row without being restricted to a linear fit across all R
          // truthR.push_back(RValuesFieldOnBeamOffSouth[module][j]); // this is with individual fits to each stripe row without being restricted to a linear fit across all R

          //truthR.push_back(SouthPos.Perp());
          phiTmp = SouthPos.Phi();
          // phi = phiTmp - 0.00148465 + 0.00219335 * cos(phiTmp - 1.24219); // this is Ben's field off ideal coordinates
          phi = phiTmp + gPhiOffsetsSouth->Eval(phiTmp); // this is from GPR

          //phi = phiTmp;
          if(phi < -TMath::Pi()) phi += 2*TMath::Pi();
          if(phi > TMath::Pi()) phi -= 2*TMath::Pi();
          truthPhi.push_back(phi);
          if(phiTmp < 0.0) phiTmp += 2*TMath::Pi();
          truthPhiOriginal.push_back(phiTmp);
          truthIndex.push_back((18+k)*10000 + (j+(module*nRadii))*100 + ((int)pos[module][j].size()-i-1));

        }
      }
    }
  }

  vector<double> phiNorth;
  vector<double> rNorth;
  vector<double> xNorth;
  vector<double> yNorth;
  vector<double> xNorthOriginal;
  vector<double> yNorthOriginal;
  vector<int> indexNorth;
  vector<double> phiSouth;
  vector<double> rSouth;
  vector<double> xSouth;
  vector<double> ySouth;
  vector<double> xSouthOriginal;
  vector<double> ySouthOriginal;
  vector<int> indexSouth;

  for(int i=0; i<(int)truthIndex.size(); i++)
  {
    if(truthIndex[i] / 10000 >= 18)
    {
      phiNorth.push_back(truthPhi[i]);
      rNorth.push_back(truthR[i]);
      xNorth.push_back(truthR[i] * cos(truthPhi[i]));
      yNorth.push_back(truthR[i] * sin(truthPhi[i]));
      xNorthOriginal.push_back(truthR[i] * cos(truthPhiOriginal[i]));
      yNorthOriginal.push_back(truthR[i] * sin(truthPhiOriginal[i]));
      indexNorth.push_back(truthIndex[i]);
    }
    else if(truthIndex[i] / 10000 < 18)
    {
      phiSouth.push_back(truthPhi[i]);
      rSouth.push_back(truthR[i]);
      xSouth.push_back(truthR[i] * cos(truthPhi[i]));
      ySouth.push_back(truthR[i] * sin(truthPhi[i]));
      xSouthOriginal.push_back(truthR[i] * cos(truthPhiOriginal[i]));
      ySouthOriginal.push_back(truthR[i] * sin(truthPhiOriginal[i]));
      indexSouth.push_back(truthIndex[i]);
    }
  }

  TGraph *patternNorth = new TGraph(xNorth.size(), &xNorth[0], &yNorth[0]);
  TGraph *patternSouth = new TGraph(xSouth.size(), &xSouth[0], &ySouth[0]);

  TCanvas *c1 = new TCanvas("c1","c1",800*4,600*4);
  gStyle->SetOptStat(0);
  patternNorth->SetTitle("Stripe Pattern North;X [cm];Y [cm]");
  patternNorth->SetMarkerStyle(20);
  patternNorth->SetMarkerColor(kRed);
  patternNorth->SetMarkerSize(0.25);
  patternNorth->Draw("AP");
  for(int i=0; i<(int)indexNorth.size(); i++)
  {
    TText *t = new TText(xNorth[i], yNorth[i], Form("%d", indexNorth[i]));
    t->SetTextSize(0.0015);
    t->Draw();
  }
  c1->SaveAs("stripePatternNorth.pdf");

  c1->Clear();
  patternSouth->SetTitle("Stripe Pattern South;X [cm];Y [cm]");
  patternSouth->SetMarkerStyle(20);
  patternSouth->SetMarkerColor(kRed);
  patternSouth->SetMarkerSize(0.25);
  patternSouth->Draw("AP");
  for(int i=0; i<(int)indexSouth.size(); i++)
  {
    TText *t = new TText(xSouth[i], ySouth[i], Form("%d", indexSouth[i]));
    t->SetTextSize(0.0015);
    t->Draw();
  }
  c1->SaveAs("stripePatternSouth.pdf");


  CDBTTree cdbttree("CMStripePattern_Apr15_2026_idealcoordinates_78242fieldoffbeamoff_updatedstripefits.root");

  for(int i=0; i<(int)truthIndex.size(); i++)
  {
    cdbttree.SetDoubleValue(truthIndex[i], "truthR", truthR[i]);
    cdbttree.SetDoubleValue(truthIndex[i], "truthPhi", truthPhi[i]);
  }

  cdbttree.Commit();
  cdbttree.WriteCDBTTree();

}
