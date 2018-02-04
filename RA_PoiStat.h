// ************************************
// R. Aggarwal, A.Caldwell
// See: https://link.springer.com/article/10.1140/epjp/i2012-12024-0
// Class to find
// a) smallest/central error intervals
// b) to plot histograms with error bands with (limit)% probability
// c) MCopt = "infinite" MC expectations are known with very small error (see paper)
//          = "" for finite MC cases
// d) use opt = "" for now ... to be kept for future changes
// ************************************

#ifndef RA_POISTAT_H
#define RA_POISTAT_H

// C C++
#include <iostream>
#include <iomanip>

// root cern
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"
#include "Rtypes.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"

using namespace std;

class RA_PoiStat : public TObject
{

 public:

  // factorial
  double ln_faci(int x);

  // constructor and destructor
  RA_PoiStat()
  {
      cout << "***Object of RA_PoiStat created***"<< endl;
  };
  ~RA_PoiStat()
  {
      cout << "***Object of RA_PoiStat destroyed***"<< endl;
  };

  // poissonian probability no uncertainty on lumi value
  double cal_P_obs_nu(int obs, double nu, double s, TString MCopt);
  // poissonian probability case with uncertainty on lumi value (s)
  double cal_P_obs_nu_s_sigmaS(int obs, double n, double s, double error);

  // to find n for whose prob >= limit
  int * cal_n_gauss(double limit,double nu);    // gaussian
  int * cal_n_gaussErf(double limit,double nu); // gaussian
  double cal_gaussP_atx(double x,double nu);    // gaussian P(x|nu)
  double cal_gaussCF_atx(double x,double nu);   // gaussian CF P(x|nu)
  int cal_n(double limit,double nu, double s, TString MCopt); // poisson with s(>=0)
  int cal_n_obs_nu_s_sigmaS(double limit,double nu,double s,double error);//poisson with error on MC

  // Katharina: Cumulative Brazilian Plot
  int cal_cumu( double obs, double nu, double s, double error, TString MCopt );
  void Plot_Cumulative(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale , Double_t Percent_error, TString ProbSet, Double_t Prob1, Double_t Prob2, Double_t Prob3, TString MCopt);

  // finding the Probability intervals
  int * Central_ProbSet(double nu, double prob, double s, double error, TString MCopt);
  int * Smallest_ProbSet(double nu, double alpha, double s, double error, TString MCopt);

  // an example to use above functions and draw probablity bands
  void Plot_w1ProbLine(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error, TString ProbSet, Double_t Prob, TString MCopt);
  void Plot_w3ProbLines(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error, TString ProbSet , Double_t Prob1, Double_t Prob2,  Double_t Prob3, TString MCopt);
  void Plot_w3ProbLines_ratio(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error, TString ProbSet , Double_t Prob1, Double_t Prob2,  Double_t Prob3, TString MCopt);
  void Plot_w3ProbLines_lin(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error, TString ProbSet , Double_t Prob1, Double_t Prob2,  Double_t Prob3, TString opt, TString MCopt);

  // this has to be the last line of the class definition
  ClassDef(RA_PoiStat,1);

 };

#endif

#ifdef RA_PoiStat_CXX

#endif
