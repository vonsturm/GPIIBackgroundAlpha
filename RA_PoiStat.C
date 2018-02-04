#define RA_POISTAT_CXX

#include "RA_PoiStat.h"


// log factorial
double RA_PoiStat::ln_faci( int x )
{
  double c = 0.;

  if( x == 0 || x == 1 ) c = 0.;
  else if( x > 1 )
  {
    for( int i = 1; i < x+1; i++ ) c += log(i);
  }

  return c;
}


// number of events just above cumulative poissonian probability limit
int RA_PoiStat::cal_n( double limit, double nu, double s, TString MCopt )
{
  int c = 1;

  int i_limit = int(nu) * 100;
  if( nu < 1. ) i_limit = 20;

  double f_i = 0.;

  for( int i = 0; i < i_limit; i++ )
  {
    f_i += cal_P_obs_nu( i, nu, s, MCopt );

    if( f_i >= limit )
    {
      c = i;
      break;
    }
  }

  return c;
}


// Katharina
// cumulative value at obs
double RA_PoiStat::cal_cumu( double obs, double nu, double s, double error, TString MCopt )
{
  double f_i = 0.;
  int i_limit = int(obs);

  if( error <= 0. )
  {
      for( int i = 0; i <= i_limit; i++ )
        f_i += cal_P_obs_nu( i, nu, s, MCopt );
  }
  else
  {
      for( int i = 0; i <= i_limit; i++ )
        f_i += cal_P_obs_nu_s_sigmaS( i, nu, s, error );
  }

  return f_i;
}


// n above gauss cumulative
int * RA_PoiStat::cal_n_gauss( double limit, double nu )
{
  double c = 1;

  int o_star = int (nu);
  double step_size  = 0.001;
  int i_limit = (int) nu / step_size;

  double f_i = 0.;

  double  ii = o_star;			// - step_size/2.;
  for( int i = 0; i < i_limit; i++ )
  {
    ii = o_star + step_size * i;
    f_i += (  exp(-(ii - o_star)*(ii - o_star)/(2.*nu) )  )
                           / (sqrt(2.*TMath::Pi()*nu)) * step_size;

    if( f_i >= limit/2 )
    {
      c = ii - o_star;
      cout << "actual up value (n+0.5) " << ii << endl;
      break;
    }
  }

  int * i_set = new int[2];

  i_set[0] = (int) (o_star - c); // see eq. in paper

  if(ii > int(ii)) i_set[1] = (int) (ii + 1); // see eq. in paper
  else             i_set[1] = (int)  ii;

  return i_set;
}


int * RA_PoiStat::cal_n_gaussErf( double limit, double nu )
{
  int c = 1;

  double f_i = 0.;
  int i_limit = int(nu) * 100;

  for( int i = 0; i < i_limit; i++ )
  {
    double x = i/sqrt(2.*nu);
    f_i = 0.5 * TMath::Erf(x);

    if( f_i >= limit/2 )
    {
      c = i;
      break;
    }
  }

  int * i_set = new int[2];
  int o_star = int (nu);

  i_set[0] = (int) (o_star - c); // see eq. in paper
  i_set[1] = (int) (o_star + c);// see eq. in paper

  return i_set;
}


double RA_PoiStat::cal_gaussP_atx( double x, double nu )
{
  int o_star = int (nu);
  double f_i = 0.;

  f_i = (  exp(-(x - o_star)*(x - o_star)/(2.*nu) )  )
    		/ (sqrt(2.*TMath::Pi()*nu)) ;
  cout << "P(x|nu)gauss" <<  f_i  <<  endl;
  return f_i;
}


double RA_PoiStat::cal_gaussCF_atx( double x, double nu )
{
  int o_star = int (nu);
  double f_i = 0.;

  for( int i = o_star; i < x+1; i++ )
  {
    f_i += ( exp(-(i - o_star)*(i - o_star)/(2.*nu) ) )
                           / (sqrt(2.*TMath::Pi()*nu));

    cout << "i=" << i << "f_i = " << f_i << endl;
  }

    cout << "CF --- > P(x|nu)gauss" <<  f_i  <<  endl;
    return f_i;
}

// poissonian probability
double RA_PoiStat::cal_P_obs_nu( int obs, double nu, double s, TString MCopt )
{
  double f_i = 0.;

  if( MCopt == "infinite" && s == 1. )
  {
    f_i += exp( - nu + ( obs * log(nu) ) - ln_faci(obs) );
  }
  else
  {
    f_i += exp( ( nu + 0.5 ) * log(s)
			- (nu+obs + 0.5)*log(s+1) - (obs*log(4.))
		     	+ ln_faci( 2.*(obs + nu))
		     	+ ln_faci(nu)
		     	- ln_faci(obs)
		     	- ln_faci(2.*nu)
		     	- ln_faci(obs+nu));
  }

  return f_i;
}

double  RA_PoiStat::cal_P_obs_nu_s_sigmaS( int Obs, double n, double s, double error )
{
  // calculates P(obs|n,s,error)
  // double for N as it is usually sum of weights in MC.
  // give error in %

  double P_case = 1.;

  int N = int(n);
  // calculating the factorial factor
  double log_x = N * log(4) + ln_faci(N) - ln_faci(2.*N);
  double x = exp(log_x);
  //----------------------------------------

  double Sigma_s = error*s / 100.; // error is in %

  double s_min = max( 0.1, s - 5*Sigma_s);
  double s_max = s + 5. * Sigma_s;
  int Ssteps = 1000;

  double ds = (s_max - s_min)/double(Ssteps);

  // take care if the guassian is truncated at 0.
  //---------------------------------------------
  double s_i = s_min - 0.5*ds ;
  double gint = 0.;
  for (int i_s = 0; i_s < Ssteps ; i_s++)
  {
      s_i += ds;
      gint += ds * exp( -pow(s_i-s,2)/(2.*Sigma_s*Sigma_s) )
                   / sqrt( 2.*TMath::Pi() ) / Sigma_s;
  }
  //---------------------------------------------

  // nu integ variables

  double nu_max_limits[21]={18., 22., 25., 27., 29., 32., 33., 35., 37., 39., 41., 43., 44., 46., 48., 50., 52., 53., 54., 56., 58};
  // keeping in mind that the senstivity of poisson factor in integral
  // is less than 10-8

  double nu_min = 0.; double nu_max = 0.; //. * Obs;
  if( Obs <= 20 )
  {
    nu_max = nu_max_limits[Obs];
  }
  else
  {
    nu_min = max(0.,Obs-5.*sqrt(Obs));
    nu_max = Obs+5.*sqrt(Obs);
  }

  int Nusteps = 1000;

  double dnu = (nu_max - nu_min)/double(Nusteps);

  // integration started ......................
  //===========================================
  double Val_nu  = 0.;
  double nu = nu_min - 0.5 * dnu;

  for (int i_nu = 0; i_nu < Nusteps ; i_nu++)
  {
    nu += dnu;
    double Val_nu_i = TMath::Exp(-nu)
                    * TMath::Power(nu,Obs)/TMath::Factorial(Obs)
                    * dnu;

    double Val_s = 0.;
    s_i = s_min - 0.5 * ds;

    for (int i_s = 0; i_s < Ssteps ; i_s++)
    {
      s_i += ds;
      double Val_s_i = s_i
                     * ( TMath::Power(s_i*nu,N-0.5) * TMath::Exp(-s_i*nu))
                     * x / sqrt(TMath::Pi())
                     * ( 1./sqrt( 2.*TMath::Pi() ) / Sigma_s)
                     * ( TMath::Exp( -(s_i - s)*(s_i - s)/(2.*Sigma_s*Sigma_s) ) )
                     * ds;

      Val_s += Val_s_i;
    } // s integ

    Val_nu += Val_nu_i * Val_s / gint;
  } // nu integ

  P_case = Val_nu;

  return Val_nu;
}


int RA_PoiStat::cal_n_obs_nu_s_sigmaS(double limit,double nu,double s,double error)
{
  int c = 0;

  int i_limit = int(nu) * 100;
  if( nu < 1. ) i_limit = 20;

  double f_i = 0.;

  for(int i = 0; i < i_limit; i++ )
  {
    f_i += cal_P_obs_nu_s_sigmaS(i,nu,s,error);

    if(f_i >= limit)
    {
      c = i;
      break;
    }
  }

  return c;
}



int * RA_PoiStat::Central_ProbSet( double nu, double alpha, double s, double error, TString MCopt )
{
  // to have a set which has alpha % central probability
  double alphaby2  = (1-alpha) / 2.;
  double prob_down = alphaby2;
  double prob_up   = 1. - alphaby2;

  int o_star     = (nu + 0.5)/s; // 2nd change here
  double nu_star = nu/s;

  int obs_low = int (nu_star);
  int obs_up  = int (nu_star);

  if(nu_star > 50)
  {
    cout << " o* = " << obs_low << " : applying Gaussian PDF " << endl;
    int * n = cal_n_gauss(alpha, nu_star);
    obs_low = n[0];
    obs_up  = n[1];
  }
  else
  {
    cout << " o* = " << o_star << " : applying Poisson PDF " << endl;
    if(error > 0.)
    {
      cout << "scale: " << s << " +- " << error << "%" <<endl;
      obs_low = cal_n_obs_nu_s_sigmaS( prob_down, nu, s, error );
      obs_up  = cal_n_obs_nu_s_sigmaS( prob_up,   nu, s, error );
    }
    else
    { // taking error == 0 case
      cout << "scale: " << s << " +- " << error << "%" <<endl;
      obs_low = cal_n( prob_down, nu, s, MCopt );
      obs_up  = cal_n( prob_up,   nu, s, MCopt );
    }
  }

  int * i_set = new int[2];

  i_set[0] = obs_low;
  i_set[1] = obs_up;

  cout << "=============================================="<< endl;
  cout << alpha << " % Central interval for "<< nu << " expectation.."<<  endl;
  cout << "[ " << i_set[0] << " : " << i_set[1] << " ]" << endl;
  cout << "=============================================="<< endl;

  return i_set;
}

// smallest probability interval
// successively add more probable values to interval until reaching the requested probability
int * RA_PoiStat::Smallest_ProbSet( double nu, double alpha, double s, double error, TString MCopt )
{
  int o_star = int (nu);
  if(MCopt != "infinite") o_star = int (nu + 0.5)/s; //changed here 3rd

  double nu_exp = nu/s; // for gaussian

  int obs_low = o_star;
  int obs_up  = o_star;

  int i_limit = int(nu) * 100;
  if(nu  < 1.) i_limit = 20;
  if(nu == 0.) i_limit = 0;

  double prob = alpha;

  if(o_star > 50)
  {
    int * n = cal_n_gauss(prob, nu_exp);
    obs_low = n[0];
    obs_up  = n[1];
  }
  else
  {
    double prob_obs = cal_P_obs_nu(obs_low, nu, s, MCopt);
    if(error > 0.) prob_obs = cal_P_obs_nu_s_sigmaS(obs_low, nu, s, error);

    if(prob_obs > prob)
    {
      obs_low = o_star;
      obs_up  = o_star;
    }
    else
    {
      double obs_low_i = obs_low;
      double obs_up_i  = obs_up;

      for(int i = 0; i < i_limit; i++)
      {
          if(obs_low > 0) obs_low_i = obs_low - 1;
          obs_up_i = obs_up + 1;

          double Prob_up_i  = 0.;
          double Prob_low_i = 0.;

          if(error > 0.)
          {
              Prob_up_i = cal_P_obs_nu_s_sigmaS(obs_up_i, nu, s, error);
              if(obs_low > 0) Prob_low_i = cal_P_obs_nu_s_sigmaS(obs_low_i, nu, s, error);
          }
          else
          { // taking the case error == 0
              Prob_up_i = cal_P_obs_nu(obs_up_i, nu, s, MCopt);
              if(obs_low > 0) Prob_low_i = cal_P_obs_nu(obs_low_i, nu, s, MCopt);
          }

          double Prob_i = Prob_up_i;

          if(Prob_up_i > Prob_low_i) obs_up = obs_up_i;
          else if(Prob_up_i < Prob_low_i)
          {
              obs_low = obs_low_i;
              Prob_i  = Prob_low_i;
          }
          else if(Prob_up_i == Prob_low_i)
          { //eg. when nu is an interger
              obs_up  = obs_up_i;
              obs_low = obs_low_i;
          }

          prob_obs += Prob_i;

          if(prob_obs >= prob) break;
      } // for loop
    } // go ahead of second loop
  }

  int * o_set = new int[2];
  o_set[0] = obs_low;
  o_set[1] = obs_up;

  cout << "=============================================="<< endl;
  cout << alpha << " % Smallest interval for " << nu << " expectation.."<<  endl;
  cout << "[ " << o_set[0] <<  " : " << o_set[1] << " ]" << endl;
  cout << "=============================================="<< endl;

  return o_set;
}


// Katharina
void  RA_PoiStat::Plot_Cumulative(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale , Double_t Percent_error,
		TString ProbSet, Double_t Prob1, Double_t Prob2, Double_t Prob3, TString MCopt)
{
    h_mc->Scale(1./Lumi_scale);

    Int_t myLineWidth=1;

    TH1D * hProb_1 = new TH1D(); h_mc->Copy(*hProb_1);
    TH1D * hProb_2 = new TH1D(); h_mc->Copy(*hProb_2);
    TH1D * hProb_3 = new TH1D(); h_mc->Copy(*hProb_3);

    TH1D * hProb_1l = new TH1D(); h_mc->Copy(*hProb_1l);
    TH1D * hProb_2l = new TH1D(); h_mc->Copy(*hProb_2l);
    TH1D * hProb_3l = new TH1D(); h_mc->Copy(*hProb_3l);

    TH1D * da_copy = new TH1D(); h_data->Copy(*da_copy);
    TH1D * mc_copy = new TH1D(); h_mc->Copy(*mc_copy);

    da_copy->SetMarkerStyle(23);
    da_copy->SetMarkerSize(1.5);
    da_copy->SetMarkerColor(kBlack);

    mc_copy->SetFillColor(kWhite);

    Int_t nBin = h_data->GetNbinsX();
    cout << " Nbins in this hist : " << nBin << endl;
    int * N_1i_prob = new int[2];
    int * N_2i_prob = new int[2];
    int * N_3i_prob = new int[2];

    for( Int_t i = 0; i < nBin; i++ )
    {
      double N_obs = h_data->GetBinContent(i+1);
      double N_exp = h_mc->GetBinContent(i+1);

      if(ProbSet == "Central")
      {
          N_1i_prob = Central_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
          N_2i_prob = Central_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
          N_3i_prob = Central_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
      }
      else if(ProbSet == "Smallest")
      {
          N_1i_prob = Smallest_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
          N_2i_prob = Smallest_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
          N_3i_prob = Smallest_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
      }

      // translate interval in cumulatvice space
      double * C_1i_prob = new double[2];
      double * C_2i_prob = new double[2];
      double * C_3i_prob = new double[2];

      C_1i_prob[0] = cal_cumu( N_1i_prob[0], N_exp, Lumi_scale, Percent_error, MCopt );
      C_2i_prob[0] = cal_cumu( N_2i_prob[0], N_exp, Lumi_scale, Percent_error, MCopt );
      C_3i_prob[0] = cal_cumu( N_3i_prob[0], N_exp, Lumi_scale, Percent_error, MCopt );

      C_1i_prob[1] = cal_cumu( N_1i_prob[1], N_exp, Lumi_scale, Percent_error, MCopt );
      C_3i_prob[1] = cal_cumu( N_3i_prob[1], N_exp, Lumi_scale, Percent_error, MCopt );
      C_2i_prob[1] = cal_cumu( N_2i_prob[1], N_exp, Lumi_scale, Percent_error, MCopt );

      double C_data   = cal_cumu( N_obs, N_exp, Lumi_scale, Percent_error, MCopt );
      double C_mc     = cal_cumu( N_exp, N_exp, Lumi_scale, Percent_error, MCopt );

      hProb_1->SetBinContent( i+1, C_1i_prob[1] );
      hProb_2->SetBinContent( i+1, C_2i_prob[1] );
      hProb_3->SetBinContent( i+1, C_3i_prob[1] );

      hProb_1l->SetBinContent( i+1, C_1i_prob[0] );
      hProb_2l->SetBinContent( i+1, C_2i_prob[0] );
      hProb_3l->SetBinContent( i+1, C_3i_prob[0] );

      da_copy->SetBinContent( i+1, C_data );
      mc_copy->SetBinContent( i+1, C_mc );

      da_copy->SetBinError( i+1, 0. );

      cout << setprecision(4) << "Nmc = " << N_exp << "\n"
           << "C_" << int(Prob1*100) << " = [ " << C_1i_prob[0] << " : " << C_1i_prob[1] << " ]" << "\n"
           << "C_" << int(Prob2*100) << " = [ " << C_2i_prob[0] << " : " << C_2i_prob[1] << " ]" << "\n"
           << "C_" << int(Prob3*100) << " = [ " << C_3i_prob[0] << " : " << C_3i_prob[1] << " ]" << "\n"
           << "C_mc(" << N_exp << ") = " << C_mc << "\n"
           << "C_data(" << N_obs << ") = " << C_data << endl;

    }

    da_copy->SetStats(0);
    da_copy->SetMarkerStyle(20);
    da_copy->SetMarkerSize(0.7);
    da_copy->SetMaximum(1.1);

    hProb_1->SetLineWidth(myLineWidth*2.0);
    hProb_1->SetLineColor(kGreen-7);
    hProb_2->SetLineColor(kYellow-7);
    hProb_3->SetLineColor(kRed-7);
    hProb_1l->SetLineColor(kGreen-7);
    hProb_2l->SetLineColor(kYellow-7);
    hProb_3l->SetLineColor(kRed-7);

    hProb_1->SetFillColor(kGreen-7);
    hProb_2->SetFillColor(kYellow-7);
    hProb_3->SetFillColor(kRed-7);
    hProb_1l->SetFillColor(kYellow-7);
    hProb_2l->SetFillColor(kYellow-7);
    hProb_3l->SetFillColor(kRed-7);

    mc_copy->SetLineColor(kBlack);
    mc_copy->SetLineWidth(myLineWidth*2.0);

    TString mcOpt = "same";

    da_copy->SetTitle("");
    da_copy->DrawCopy("p");

    hProb_3->DrawCopy(mcOpt+"hist");
    hProb_2->DrawCopy(mcOpt+"hist");
    hProb_1->DrawCopy(mcOpt+"hist");
    hProb_1l->DrawCopy(mcOpt+"hist");
    hProb_2l->DrawCopy(mcOpt+"hist");
    hProb_3l->DrawCopy(mcOpt+"hist");

//    da_copy->DrawCopy("Sameep");
//    da_copy->DrawCopy("SameAxis");
    mc_copy->DrawCopy(mcOpt+"hist");
    da_copy->DrawCopy("Same p");

    return;
}


void  RA_PoiStat::Plot_w1ProbLine(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale , Double_t Percent_error,
		TString ProbSet, Double_t Prob, TString MCopt)
{
  h_mc->Scale(1./Lumi_scale);

  Int_t myLineWidth=1;

  Double_t tmpMaxData = h_data->GetBinContent( h_data->GetMaximumBin() );
  Double_t tmpMaxMC   = h_mc->GetMaximum();
  Double_t tmpMax     = TMath::Max( (tmpMaxData+TMath::Sqrt(tmpMaxData)), tmpMaxMC );

  TH1D * hProb_1  = new TH1D(); h_mc->Copy(*hProb_1);
  TH1D * hProb_1l = new TH1D(); h_mc->Copy(*hProb_1l);

  TH1D * da_copy  = new TH1D(); h_data->Copy(*da_copy);
  TH1D * mc_copy  = new TH1D(); h_mc->Copy(*mc_copy);

  da_copy->SetMarkerStyle(23);
  da_copy->SetMarkerSize(1.5);
  da_copy->SetMarkerColor(kBlack);

  const Int_t nBin = h_data->GetNbinsX();
  int * N_i_prob = new int[2];
  cout << " Nbins in this hist : " << nBin << endl;

  for( Int_t i = 0; i < nBin; i++ )
  {
    double N_obs = h_data->GetBinContent(i+1);
    double N_exp = h_mc->GetBinContent(i+1);

    if(ProbSet == "Central")
        N_i_prob = Central_ProbSet( N_exp, Prob, Lumi_scale, Percent_error, MCopt );
    else if(ProbSet == "Smallest")
        N_i_prob = Smallest_ProbSet( N_exp, Prob, Lumi_scale, Percent_error, MCopt );

    double N_1Prob = N_i_prob[1] + 0.5;
    if(N_i_prob[1] == 0)N_1Prob = 0.1;

    double N_1lProb = N_i_prob[0] +0.5;
    if(N_i_prob[0] == 0)N_1lProb = 0.1;

    hProb_1->SetBinContent(i+1, N_1Prob);
    hProb_1l->SetBinContent(i+1, N_1lProb);

    h_data->SetBinError(i+1,0.);

    // trick to have some representation for 0 data events
    if(N_obs > 0) da_copy->SetBinContent(i+1,0.);
    else          da_copy->SetBinContent(i+1,0.1);
  }

  h_data->SetStats(0);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.7);
  h_data->SetMaximum(1.1*tmpMax);
  h_data->SetMinimum(0.48);

  hProb_1->SetLineWidth(myLineWidth*2.0);
  hProb_1->SetLineColor(kGray);
  hProb_1l->SetLineColor(kGray);

  hProb_1->SetFillColor(kGray);
  h_mc->SetFillColor(kGray);
  hProb_1l->SetFillColor(kYellow-12);

  h_mc->SetLineColor(kBlack);
  h_mc->SetLineWidth(myLineWidth *2.0);

  TString mcOpt = "same";

  h_data->Draw("p");
  hProb_1->DrawCopy(mcOpt+"hist");
  h_mc->DrawCopy(mcOpt+"hist");
  hProb_1l->DrawCopy(mcOpt+"hist");

  h_data->DrawCopy("Sameep");
  h_data->DrawCopy("SameAxis");
}


void  RA_PoiStat::Plot_w3ProbLines(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error,
		TString ProbSet , Double_t Prob1, Double_t Prob2,  Double_t Prob3, TString MCopt)
{
  h_mc->Scale(1./Lumi_scale);

  Int_t myLineWidth=1;

  Double_t tmpMaxData = h_data->GetBinContent( h_data->GetMaximumBin() );
  Double_t tmpMaxMC   = h_mc->GetMaximum();
  Double_t tmpMax     = TMath::Max( (tmpMaxData+TMath::Sqrt(tmpMaxData) ), tmpMaxMC);

  TH1D * hProb_1 = new TH1D(); h_mc->Copy(*hProb_1);
  TH1D * hProb_2 = new TH1D(); h_mc->Copy(*hProb_2);
  TH1D * hProb_3 = new TH1D(); h_mc->Copy(*hProb_3);

  TH1D * hProb_1l = new TH1D(); h_mc->Copy(*hProb_1l);
  TH1D * hProb_2l = new TH1D(); h_mc->Copy(*hProb_2l);
  TH1D * hProb_3l = new TH1D(); h_mc->Copy(*hProb_3l);

  TH1D * da_copy = new TH1D(); h_data->Copy(*da_copy);
  TH1D * mc_copy = new TH1D(); h_mc->Copy(*mc_copy);

  da_copy->SetMarkerStyle(23);
  da_copy->SetMarkerSize(1.5);
  da_copy->SetMarkerColor(kBlack);

  Int_t nBin = h_data->GetNbinsX();
  cout << " Nbins in this hist : " << nBin << endl;
  int * N_1i_prob = new int[2];
  int * N_2i_prob = new int[2];
  int * N_3i_prob = new int[2];

  for( Int_t i = 0; i < nBin; i++ )
  {
    double N_obs = h_data->GetBinContent(i+1);
    double N_exp = h_mc->GetBinContent(i+1);

    if(ProbSet == "Central")
    {
        N_1i_prob = Central_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Central_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Central_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }
    else if(ProbSet == "Smallest")
    {
        N_1i_prob = Smallest_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Smallest_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Smallest_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }

    double N_1Prob  = N_1i_prob[1] + 0.5;
    double N_1lProb = N_1i_prob[0] + 0.5;

    double N_2Prob  = N_2i_prob[1] + 0.5;
    double N_2lProb = N_2i_prob[0] + 0.5;

    double N_3Prob  = N_3i_prob[1] + 0.5;
    double N_3lProb = N_3i_prob[0] + 0.5;

    if( N_1i_prob[0] == 0) N_1lProb = 0.1;
    if( N_2i_prob[0] == 0) N_2lProb = 0.1;
    if( N_3i_prob[0] == 0) N_3lProb = 0.1;

    hProb_1->SetBinContent(i+1,N_1Prob);
    hProb_2->SetBinContent(i+1,N_2Prob);
    hProb_3->SetBinContent(i+1,N_3Prob);

    hProb_1l->SetBinContent(i+1,N_1lProb);
    hProb_2l->SetBinContent(i+1,N_2lProb);
    hProb_3l->SetBinContent(i+1,N_3lProb);

    da_copy->SetBinError( i+1, 0. );

    cout << "Nmc = "                         << N_exp
         << " N_" << int(Prob1*100) << " = " << N_1Prob
         << " N_" << int(Prob2*100) << " = " << N_2Prob
         << " N_" << int(Prob3*100) << " = " << N_3Prob << endl;
  }

  da_copy->SetStats(0);
  da_copy->SetMarkerStyle(20);
  da_copy->SetMarkerSize(0.7);
  da_copy->SetMaximum(1.5*tmpMaxData);

  da_copy->SetMinimum(0.01);

  hProb_1->SetLineWidth(myLineWidth*2.0);
  hProb_1->SetLineColor(kGreen-7);
  hProb_2->SetLineColor(kYellow-7);
  hProb_3->SetLineColor(kRed-7);
  hProb_1l->SetLineColor(kYellow-7);
  hProb_2l->SetLineColor(kRed-7);
  hProb_3l->SetLineColor(kWhite);

  hProb_1->SetFillColor(kGreen-7);
  hProb_2->SetFillColor(kYellow-7);
  hProb_3->SetFillColor(kRed-7);
  hProb_1l->SetFillColor(kYellow-7);
  hProb_2l->SetFillColor(kRed-7);
  hProb_3l->SetFillColor(kWhite);

  TString mcOpt = "same";

  da_copy->SetTitle("");
  da_copy->DrawCopy("p");

  hProb_3->DrawCopy(mcOpt+"hist");
  hProb_2->DrawCopy(mcOpt+"hist");
  hProb_1->DrawCopy(mcOpt+"hist");
  mc_copy->DrawCopy(mcOpt+"hist");
  hProb_1l->DrawCopy(mcOpt+"hist");
  hProb_2l->DrawCopy(mcOpt+"hist");
  hProb_3l->DrawCopy(mcOpt+"hist");

  da_copy->DrawCopy("Sameep");
  da_copy->DrawCopy("SameAxis");
}


// ratio plot
void  RA_PoiStat::Plot_w3ProbLines_ratio(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale,Double_t Percent_error,
    TString ProbSet, Double_t Prob1, Double_t Prob2, Double_t Prob3, TString MCopt)
{
  h_mc->Scale(1./Lumi_scale);

  Int_t myLineWidth=1;

  Double_t tmpMaxData = h_data->GetBinContent( h_data->GetMaximumBin() );
  Double_t tmpMaxMC   = h_mc->GetMaximum();
  Double_t tmpMax     = TMath::Max( (tmpMaxData+TMath::Sqrt(tmpMaxData)), tmpMaxMC );

  TH1D * hProb_1 = new TH1D(); h_mc->Copy(*hProb_1);
  TH1D * hProb_2 = new TH1D(); h_mc->Copy(*hProb_2);
  TH1D * hProb_3 = new TH1D(); h_mc->Copy(*hProb_3);

  TH1D * hProb_1l = new TH1D(); h_mc->Copy(*hProb_1l);
  TH1D * hProb_2l = new TH1D(); h_mc->Copy(*hProb_2l);
  TH1D * hProb_3l = new TH1D(); h_mc->Copy(*hProb_3l);

  Int_t nBin = h_data->GetNbinsX();
  cout << " Nbins in this hist : " << nBin << endl;
  int * N_1i_prob = new int[2];
  int * N_2i_prob = new int[2];
  int * N_3i_prob = new int[2];

  for( Int_t i = 0; i < nBin; i++ )
  {
    double N_obs = h_data->GetBinContent(i+1);
    double N_exp = h_mc->GetBinContent(i+1);

    if(ProbSet == "Central")
    {
        N_1i_prob = Central_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Central_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Central_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }
    else if(ProbSet == "Smallest")
    {
        N_1i_prob = Smallest_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Smallest_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Smallest_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }

    double N_1lProb = 0, N_2lProb = 0, N_3lProb = 0;

    double N_1Prob = N_1i_prob[1] + 0.5;
    if( N_1i_prob[0] != 0 ) N_1lProb = N_1i_prob[0] - 0.5;

    double N_2Prob = N_2i_prob[1] + 0.5;
    if(N_2i_prob[0] != 0) N_2lProb =  N_2i_prob[0] - 0.5;

    double N_3Prob = N_3i_prob[1] + 0.5;
    if(N_3i_prob[0] != 0) N_3lProb =N_3i_prob[0] - 0.5;

    hProb_1->SetBinContent( i+1, N_1Prob/N_exp );
    hProb_2->SetBinContent( i+1, N_2Prob/N_exp );
    hProb_3->SetBinContent( i+1, N_3Prob/N_exp );

    hProb_1l->SetBinContent( i+1, N_1lProb/N_exp );
    hProb_2l->SetBinContent( i+1, N_2lProb/N_exp );
    hProb_3l->SetBinContent( i+1, N_3lProb/N_exp );

    h_data->SetBinError(i+1,0.);
    h_data->SetBinContent(i+1,N_obs/N_exp);
    h_mc->SetBinContent(i+1,N_exp/N_exp);

    cout << "Nmc = "                         << N_exp
         << " N_" << int(Prob1*100) << " = " << N_1Prob
         << " N_" << int(Prob2*100) << " = " << N_2Prob
         << " N_" << int(Prob3*100) << " = " << N_3Prob << endl;

    cout << "Nmc = "                         << N_exp
         << " N_" << int(Prob1*100) << " = " << N_1lProb
         << " N_" << int(Prob2*100) << " = " << N_2lProb
         << " N_" << int(Prob3*100) << " = " << N_3lProb << endl;
  }

  tmpMaxData = h_data->GetBinContent( h_data->GetMaximumBin() );

  h_data->SetStats(0);
  h_data->SetMarkerStyle(24);
  h_data->SetMarkerSize(0.7);
  h_data->SetMaximum(0.5*tmpMaxData);

  h_data->SetMinimum(-0.009);

  hProb_1->SetLineWidth(myLineWidth*2.0);
  hProb_1->SetLineColor(kGreen-7);
  hProb_2->SetLineColor(kYellow-7);
  hProb_3->SetLineColor(kRed-7);
  hProb_1l->SetLineColor(kGreen-7);
  hProb_2l->SetLineColor(kYellow-7);
  hProb_3l->SetLineColor(kRed-7);

  hProb_1->SetFillColor(kGreen-7);
  hProb_2->SetFillColor(kYellow-7);
  hProb_3->SetFillColor(kRed-7);
  hProb_1l->SetFillColor(kGreen-7);
  hProb_2l->SetFillColor(kYellow-7);
  hProb_3l->SetFillColor(kRed-7);

  h_mc->SetLineColor(kBlack);
  h_mc->SetLineWidth(myLineWidth*2.0);

  TString mcOpt = "same";

  h_data->SetMaximum(25.);
  h_data->SetMinimum(-0.009);
  h_data->Draw("p");

  hProb_3->DrawCopy(mcOpt+"hist");
  hProb_2->DrawCopy(mcOpt+"hist");
  hProb_1->DrawCopy(mcOpt+"hist");
  h_mc->DrawCopy(mcOpt+"hist");
  hProb_1l->DrawCopy(mcOpt+"hist");
  hProb_2l->DrawCopy(mcOpt+"hist");
  hProb_3l->DrawCopy(mcOpt+"hist");

  h_data->DrawCopy("Sameep");
  h_data->DrawCopy("SameAxis");
}


// linear plot
void  RA_PoiStat::Plot_w3ProbLines_lin(TH1D* h_mc ,TH1D* h_data, Double_t Lumi_scale, Double_t Percent_error,
    TString ProbSet , Double_t Prob1, Double_t Prob2,  Double_t Prob3, TString opt, TString MCopt)
{
  h_mc->Scale(1./Lumi_scale);

  Int_t myLineWidth=1;

  Double_t tmpMaxData = h_data->GetBinContent( h_data->GetMaximumBin() );
  Double_t tmpMaxMC   = h_mc->GetMaximum();
  Double_t tmpMax     = TMath::Max( (tmpMaxData + TMath::Sqrt(tmpMaxData)), tmpMaxMC );

  TH1D * hProb_1 = new TH1D(); h_mc->Copy(*hProb_1);
  TH1D * hProb_2 = new TH1D(); h_mc->Copy(*hProb_2);
  TH1D * hProb_3 = new TH1D(); h_mc->Copy(*hProb_3);

  TH1D * hProb_1l = new TH1D(); h_mc->Copy(*hProb_1l);
  TH1D * hProb_2l = new TH1D(); h_mc->Copy(*hProb_2l);
  TH1D * hProb_3l = new TH1D(); h_mc->Copy(*hProb_3l);

  Int_t nBin = h_data->GetNbinsX();
  cout << " Nbins in this hist : " << nBin << endl;
  int * N_1i_prob = new int[2];
  int * N_2i_prob = new int[2];
  int * N_3i_prob = new int[2];

  for( Int_t i = 0; i < nBin; i++ )
  {
    double N_obs = h_data->GetBinContent(i+1);
    double N_exp = h_mc->GetBinContent(i+1);

    if(ProbSet == "Central")
    {
        N_1i_prob = Central_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Central_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Central_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }
    else if(ProbSet == "Smallest")
    {
        N_1i_prob = Smallest_ProbSet(N_exp, Prob1, Lumi_scale, Percent_error, MCopt);
        N_2i_prob = Smallest_ProbSet(N_exp, Prob2, Lumi_scale, Percent_error, MCopt);
        N_3i_prob = Smallest_ProbSet(N_exp, Prob3, Lumi_scale, Percent_error, MCopt);
    }

    // see paper plot half a bin more
    double N_1Prob  = N_1i_prob[1] + 0.5;
    double N_1lProb = N_1i_prob[0] + 0.5;

    double N_2Prob  = N_2i_prob[1] + 0.5;
    double N_2lProb = N_2i_prob[0] + 0.5;

    double N_3Prob  = N_3i_prob[1] + 0.5;
    double N_3lProb = N_3i_prob[0] + 0.5;

    if( N_1i_prob[0] == 0 ) N_1lProb = 0.1;
    if( N_2i_prob[0] == 0 ) N_2lProb = 0.1;
    if( N_3i_prob[0] == 0 ) N_3lProb = 0.1;

    hProb_1->SetBinContent( i+1, N_1Prob );
    hProb_2->SetBinContent( i+1, N_2Prob );
    hProb_3->SetBinContent( i+1, N_3Prob );

    hProb_1l->SetBinContent( i+1, N_1lProb );
    hProb_2l->SetBinContent( i+1, N_2lProb );
    hProb_3l->SetBinContent( i+1, N_3lProb );

    h_data->SetBinError( i+1, 0. );

    cout << "Nmc = "                         << N_exp
         << " N_" << int(Prob1*100) << " = " << N_1Prob
         << " N_" << int(Prob2*100) << " = " << N_2Prob
         << " N_" << int(Prob3*100) << " = " << N_3Prob << endl;
  }

  h_data->SetStats(0);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.9);
  h_data->SetMaximum(1.5*tmpMax);
  h_data->SetMinimum(-0.1);

  hProb_1->SetLineWidth(myLineWidth*2.0);
  hProb_1->SetLineColor(kGreen-7);
  hProb_2->SetLineColor(kYellow-7);
  hProb_3->SetLineColor(kRed-7);
  hProb_1l->SetLineColor(kGreen-7);
  hProb_2l->SetLineColor(kYellow-7);
  hProb_3l->SetLineColor(kRed-7);

  hProb_1->SetFillColor(kGreen-7);
  hProb_2->SetFillColor(kYellow-7);
  hProb_3->SetFillColor(kRed-7);
  hProb_1l->SetFillColor(kGreen-7);
  hProb_2l->SetFillColor(kYellow-7);
  hProb_3l->SetFillColor(kRed-7);

  h_mc->SetLineColor(kBlack);
  h_mc->SetLineWidth(myLineWidth*2.0);

  TString mcOpt = "same";

  h_data->Draw("p");

  hProb_3->DrawCopy(mcOpt+"hist");
  hProb_2->DrawCopy(mcOpt+"hist");
  hProb_1->DrawCopy(mcOpt+"hist");
  h_mc->DrawCopy(mcOpt+"hist");
  hProb_1l->DrawCopy(mcOpt+"hist");
  hProb_2l->DrawCopy(mcOpt+"hist");
  hProb_3l->DrawCopy(mcOpt+"hist");

  h_data->DrawCopy("Sameep");
  h_data->DrawCopy("SameAxis");
  h_mc->DrawCopy(mcOpt+"hist");
}
