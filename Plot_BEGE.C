////////////////////////////////////
/////  N.Becerici-Schmidt   ////////
//////// 10 October 2012 /////////// 
////////////////////////////////////

#include <TROOT.h>
#include <TUnixSystem.h>

#include <iostream>  
#include <iomanip>

#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

#include <RA_PoiStat.h>
#include <fstream>

using namespace std;

void Plot_BEGE( string filename, string outfilename )
{

// style settings
//.............................................................................//
  gROOT->SetStyle("Plain");
  gStyle->SetLabelFont(42,"xy");
  gStyle->SetTitleFont(42,"xy");
  gStyle->SetOptFit(111111);
  gStyle->SetOptStat(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLabelSize(0.07,"xy");
  gStyle->SetTitleSize(0.07,"xy");
  gStyle->SetLabelOffset(0.02,"x");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetTitleOffset(1.1,"x");
  //   gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleOffset(0.72,"y");

  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.11);

  gROOT->ForceStyle();
//.............................................................................//

  TFile* fitfile = new TFile( filename.c_str() ); 


  // read in the histograms

  TH1D* hdata=(TH1D*)fitfile->Get("henergySum");
//   hdata->SetMarkerStyle(20);
//   hdata->SetMarkerSize(0.8);
  hdata->SetFillColor(kGray);
  hdata->SetLineColor(kGray);
  hdata->GetYaxis()->SetTitle("Counts/(50 keV)");
  hdata->GetXaxis()->SetTitle("E (keV)");


  TH1D* hMC=(TH1D*)fitfile->Get("hMC");
  hMC->SetLineColor(1);
  hMC->SetLineWidth(2);
  hMC->SetLineStyle(1);

  TH1D* hPo210=(TH1D*)fitfile->Get("h_Po210_pPlus_dl600nm");
  hPo210->SetLineColor(kCyan);
  hPo210->SetLineStyle(1);
  hPo210->SetLineWidth(2);

  TH1D* hRa226pPlus=(TH1D*)fitfile->Get("h_Ra226_pPlus_dl600nm");
  TH1D* hRn222pPlus600=(TH1D*)fitfile->Get("h_Rn222_pPlus_dl600nm");
  hRa226pPlus->Add( hRn222pPlus600 );

//  TH1D* hPo218pPlus=(TH1D*)fitfile->Get("h_Po218_pPlus");
//  hRa226pPlus->Add(hPo218pPlus);
//  TH1D* hPo214pPlus=(TH1D*)fitfile->Get("h_Po214_pPlus");
//  hRa226pPlus->Add(hPo214pPlus);

  hRa226pPlus->SetLineColor(kViolet-3);
  hRa226pPlus->SetLineStyle(2);
  hRa226pPlus->SetLineWidth(3);

  TH1D* hRa226LAr=(TH1D*)fitfile->Get("h_Ra226_inLArBH");
  TH1D* hRn222LAr=(TH1D*)fitfile->Get("h_Rn222_inLArBH");
  hRa226LAr->Add(hRn222LAr);
//  TH1D* hPo218LAr=(TH1D*)fitfile->Get("h_Po218_inLArBH");
//  hRa226LAr->Add(hPo218LAr);
//  TH1D* hPo214LAr=(TH1D*)fitfile->Get("h_Po214_inLArBH");
//  hRa226LAr->Add(hPo214LAr);

  hRa226LAr->SetLineColor(kMagenta);
  hRa226LAr->SetLineStyle(6);
  hRa226LAr->SetLineWidth(3);


  TCanvas * c1= new TCanvas("c1", "c1", 5, 5, 600, 500);
  TPad *pad1 = new TPad("pad1",
  "The pad with the function",0.02,0.40,0.98,0.98);
  TPad *pad2 = new TPad("pad2",
  "The pad with the histogram",0.02,0.02,0.98,0.40);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetTickx(1);
  pad1->SetLogy(1);
  pad1->SetTopMargin(0.05);
  pad1->SetBottomMargin(0);


  hdata->GetXaxis()->SetRangeUser(3500.,5300.);
  hdata->GetYaxis()->SetRangeUser(0.01,200.);

  hdata->Draw("");
  hMC->Draw("same");
  hPo210->Draw("same");
  hRa226pPlus->Draw("same");
  hRa226LAr->Draw("same");


  TLegend * myLeg1 = new TLegend(0.31,0.75,0.42,0.91,"");
  myLeg1->SetBorderSize(0);
  myLeg1->SetFillColor(kWhite);
  myLeg1->SetTextSize(0.07);
  myLeg1->SetTextFont(42); 
  myLeg1->AddEntry(hdata,"data","f");
  myLeg1->AddEntry(hMC,"model","L");
  myLeg1->Draw();

  TLegend * myLeg2 = new TLegend(0.51,0.67,0.69,0.91,"");
  myLeg2->SetBorderSize(0);
  myLeg2->SetFillColor(kWhite);
  myLeg2->SetTextSize(0.07);
  myLeg2->SetTextFont(42); 
  myLeg2->AddEntry(hPo210,"^{210}Po on surface","L");
  myLeg2->AddEntry(hRa226pPlus,"^{226}Ra & daughters on surface ","L");
  myLeg2->AddEntry(hRa226LAr,"^{226}Ra & daughters in LAr ","L");
  myLeg2->Draw();

  TLatex tl;
  tl.SetTextSize(0.06);
  tl.SetTextFont(22);
//   tl.DrawLatex(620., 3764, "golden data set");
  tl.DrawLatex(1850., 50., "BEGe sum data set");


// draw color bands
  pad2->cd();
  pad2->SetTickx(1);
  pad2->SetBottomMargin(0.25);
  pad2->SetTopMargin(0);

  TH1D * hmodel1 = new TH1D();
  hMC->Copy(*hmodel1);
  hmodel1->GetXaxis()->SetRangeUser(3500., 5300.); 
  TH1D * hdata1 = new TH1D();
  hdata->Copy(*hdata1);
  hdata1->GetXaxis()->SetRangeUser(3500., 5300.); 
  hdata1->GetYaxis()->SetTitle("data/model ratio ");
  hdata1->GetXaxis()->SetTitle("E (keV)");
  hdata1->SetTitleOffset(0.50,"y");
  hdata1->SetLabelSize(0.11,"xy");
  hdata1->SetTitleSize(0.11,"xy");
  hdata1->GetYaxis()->SetDecimals(1);
  hdata1->GetYaxis()->SetNdivisions(505);

  hdata1->GetXaxis()->SetNdivisions(105);
  hdata1->SetMarkerStyle(24);
  hdata1->SetMarkerSize(0.7);
 

  TLegend * myLeg11 = new TLegend(0.15,0.83,0.25,0.96,"");
  myLeg11->SetBorderSize(0);
  myLeg11->SetFillColor(kWhite);
  myLeg11->SetTextSize(0.095);
  myLeg11->SetTextFont(42);
  myLeg11->AddEntry(hdata1,"data/model","P");

  TLegend * myLeg12 = new TLegend(0.38,0.83,0.48,0.96,"");
  myLeg12->SetBorderSize(0);
  myLeg12->SetFillColor(kWhite);
  myLeg12->SetTextSize(0.095);
  myLeg12->SetTextFont(42);
  TH1D*h1=new TH1D();
  h1->SetFillColor(kGreen-4);
  myLeg12->AddEntry(h1,"68%","F");

  TLegend * myLeg13 = new TLegend(0.38,0.71,0.48,0.84,"");
  myLeg13->SetBorderSize(0);
  myLeg13->SetFillColor(kWhite);
  myLeg13->SetTextSize(0.095);
  myLeg13->SetTextFont(42);
  TH1D*h2=new TH1D();
  h2->SetFillColor(kYellow-4);
  myLeg13->AddEntry(h2,"95%","F");

  TLegend * myLeg14 = new TLegend(0.38,0.59,0.48,0.72,"");
  myLeg14->SetBorderSize(0);
  myLeg14->SetFillColor(kWhite);
  myLeg14->SetTextSize(0.095);
  myLeg14->SetTextFont(42);
  TH1D*h3=new TH1D();
  h3->SetFillColor(kRed-7);
  myLeg14->AddEntry(h3,"99.9%","F");

  RA_PoiStat * h2_=  new RA_PoiStat();
  h2_->RA_PoiStat::Plot_w3ProbLines_ratio(hmodel1,hdata1,1.,0.,"Smallest",0.68,0.95,0.999,"infinite");
  myLeg11->Draw();
  myLeg12->Draw();
  myLeg13->Draw();
  myLeg14->Draw();

  TFile* file=new TFile( outfilename.c_str() , "RECREATE");

  c1->Write();

  file->Close();

  outfilename.replace( outfilename.end()-4, outfilename.end(), "svg" );
  c1->SaveAs( outfilename.c_str() );
}
