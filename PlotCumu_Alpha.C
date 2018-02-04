/*
 * Author: Katharina von Sturm 4.2.2018
 *
 *
 */

 // C/C++
 #include <iostream>
 #include <iomanip>
 #include <fstream>

// root cern
#include "TROOT.h"
#include "TUnixSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include "RA_PoiStat.h"

using namespace std;

void load_style();

void PlotCumu_Alpha( string filename, string outfilename )
{
    load_style();

    RA_PoiStat * h2_ = new RA_PoiStat();

    TFile * fitfile = new TFile( filename.c_str() );

    TH1D* h_data = (TH1D*) fitfile -> Get("hSum");
    h_data->GetXaxis()->SetRangeUser( 3500., 5300. );
    h_data->GetYaxis()->SetRangeUser( 0.01, 200. );

    TH1D* h_mc   = (TH1D*) fitfile -> Get("hMC");

    TCanvas * c1= new TCanvas("c1", "c1", 5, 5, 600, 750);
    TPad *pad1 = new TPad("pad1", "The pad with the data",0.02,0.50,0.98,0.98);
    TPad *pad2 = new TPad("pad2", "The pad with the brazilian plot",0.02,0.25,0.98,0.50);
    TPad *pad3 = new TPad("pad3", "The pad with the cumulative brazilian plot",0.02,0.02,0.98,0.25);
    pad1->Draw(); pad2->Draw(); pad3->Draw();

    pad1->cd(); pad1->SetTickx(1); pad1->SetLogy(1); pad1->SetTopMargin(0.05); pad1->SetBottomMargin(0);

    h_data -> Draw("p");
    h_mc   -> Draw("histsame");

    TLegend * myLeg1 = new TLegend(0.31,0.75,0.42,0.91,"");
    myLeg1->SetBorderSize(0);
    myLeg1->SetFillColor(kWhite);
    myLeg1->SetTextSize(0.07);
    myLeg1->SetTextFont(42);
    myLeg1->AddEntry(h_data,"data","f");
    myLeg1->AddEntry(h_mc,"model","L");
    myLeg1->Draw();

  //   TLatex tl;
  //   tl.SetTextSize(0.06);
  //   tl.SetTextFont(22);
  // //   tl.DrawLatex(620., 3764, "golden data set");
  //   tl.DrawLatex(1850., 50., "BEGe sum data set");

    pad2->cd(); pad2->SetTickx(1); pad2->SetBottomMargin(0); pad2->SetTopMargin(0);
    h2_ -> Plot_w3ProbLines( h_mc, h_data, 1., 0., "Smallest", 0.683, 0.954, 0.997, "infinite" );

    pad3->cd(); pad3->SetTickx(1); pad3->SetBottomMargin(0.25); pad3->SetTopMargin(0);
    h2_ -> Plot_Cumulative( h_mc, h_data, 1., 0., "Smallest", 0.683, 0.954, 0.997, "infinite" );

    // Draw Legend
    TLegend * myLeg11 = new TLegend(0.15,0.83,0.25,0.96,"");
    myLeg11->SetBorderSize(0);
    myLeg11->SetFillColor(kWhite);
    myLeg11->SetTextSize(0.095);
    myLeg11->SetTextFont(42);
    myLeg11->AddEntry(h_data,"data/model","P");

    TLegend * myLeg12 = new TLegend(0.38,0.83,0.48,0.96,"");
    myLeg12->SetBorderSize(0);
    myLeg12->SetFillColor(kWhite);
    myLeg12->SetTextSize(0.095);
    myLeg12->SetTextFont(42);
    TH1D * h1 = new TH1D(); h1->SetFillColor(kGreen-7);
    myLeg12->AddEntry(h1,"1#sigma","F");

    TLegend * myLeg13 = new TLegend(0.38,0.71,0.48,0.84,"");
    myLeg13->SetBorderSize(0);
    myLeg13->SetFillColor(kWhite);
    myLeg13->SetTextSize(0.095);
    myLeg13->SetTextFont(42);
    TH1D * h2 = new TH1D(); h2->SetFillColor(kYellow-7);
    myLeg13->AddEntry(h2,"2#sigma","F");

    TLegend * myLeg14 = new TLegend(0.38,0.59,0.48,0.72,"");
    myLeg14->SetBorderSize(0);
    myLeg14->SetFillColor(kWhite);
    myLeg14->SetTextSize(0.095);
    myLeg14->SetTextFont(42);
    TH1D * h3 = new TH1D(); h3->SetFillColor(kRed-7);
    myLeg14->AddEntry(h3,"3#sigma","F");

//    myLeg11->Draw();
//    myLeg12->Draw();
//    myLeg13->Draw();
//    myLeg14->Draw();

    TFile * file = new TFile( outfilename.c_str() , "RECREATE");
    c1->Write();
    file->Close();

    outfilename.replace( outfilename.end()-4, outfilename.end(), "png" );
    c1->SaveAs( outfilename.c_str() );

    return;
}


void load_style()
{
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
}
