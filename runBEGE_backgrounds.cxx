// ***************************************************************
// This file was created using the CreateProject.sh script
// for project 2nbb.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCEngineMCMC.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCParameter.h>

#include "TStyle.h"

#include "BEGE_backgrounds.h"

using namespace std;

int main()
{

  // --------- set input parameters ------

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();
  gStyle->SetOptStat(0);
  
  // open log file
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);
  
  // create new BEGE_backgrounds object
  BEGE_backgrounds* m = new BEGE_backgrounds();
   
  // create a new summary tool object
  BCSummaryTool * summary = new BCSummaryTool(m);

  BCLog::OutSummary("Model created");
  
  // ==== Set precision for fit ====
  
  //--->
  //---- SET HERE ----
  const int precisionValue = 3;
  //------------------

  string precisionString="";
  BCEngineMCMC::Precision precision;
  if(precisionValue==0)
    {
      precisionString="kLow";
      precision=BCEngineMCMC::kLow;
    }
  if(precisionValue==1)
    {
      precisionString="kMedium";
      precision=BCEngineMCMC::kMedium;
    }
  if(precisionValue==2)
    {
      precisionString="kHigh";
      precision=BCEngineMCMC::kHigh;
    }
  if(precisionValue==3)
    {
      precisionString="kVeryHigh";
      precision=BCEngineMCMC::kVeryHigh;
    }

  cout << "Precision: " << precisionString << endl;

  // set the precision for the fit
  m->MCMCSetPrecision(precision);

  // === Set the histograms for the data spectra ===

  //--->
  //---- SET HERE ----
  double hemin=3500.; //try also 2900
  double hemax=5300.;
  double binwidth=50.;
  //------------------

  int hnumbins=(int)((hemax-hemin)/binwidth);

  m->SetHistogramParameters(hnumbins, hemin, hemax);

  cout << endl;
  cout << "**********************************************" << endl;
  cout << "ENERGY RANGE AND BINWIDTH INFORMATION" << endl;
  cout << "Energy range: (" << hemin << " - " << hemax << ")keV" << endl;
  cout << "Binwidth: " << binwidth << endl;
  cout << "**********************************************" << endl;
  cout << endl;

  // --- define data input ---
  // read in the names of the files you want to use, 
  vector<int> runlist = { 53,54,55,56,57,58,59,60,61,62,63 };
  m->ReadDataEnrBEGe( runlist );

  // --- read in the MC histograms ---
  m->ReadMCAlpha();

  // --- define parameters ---
  m->DefineParameters();

  // to get better sensitivity on a certain parameter
//   m->SetNbins("inverseHalflife_BEGE_backgrounds",2000);
//  m->SetNbins("halflife_2nbb",2000);
  // Set parameter binning
  size_t nPars = m->GetNParameters();
  cout << nPars << " Parameters in the fit." << endl;
  for ( size_t i = 0; i < nPars; i++ )  
  {
	m->GetParameter(i)->SetNbins(500);
  	cout << "Binning par " << i << ": " << m->GetParameter(i)->GetNbins() << endl;
  }

  //m->SetNbins(2000);


//   size_t nPars = m->GetNParameters();
//   for (size_t i=0;i<nPars;i++)        
//     m->MCMCSetFlagFillHistograms(i,false); //disable histos for nuisance parameters

  // create new output object
  string OUTPUT_DIR = "/opt/exp_software/gerda/gerda_gpfs/vonsturm/BAT/BEGE_alphas/ModelOutput/";
  OUTPUT_DIR += precisionString;

  BCModelOutput* mout = new BCModelOutput( m, Form( "%s/BEGE_alphas_model.root", OUTPUT_DIR.c_str() ) );

  // switch writing of Markov Chains on
  mout->WriteMarkovChain(true);

  // normalize the posterior, i.e. integrate posterior
  // over the full parameter space
  m->Normalize();

  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters
  
  m->MarginalizeAll();

  // run mode finding; by default using Minuit
  //   m->FindMode();
  
  // if MCMC was run before (MarginalizeAll()) it is
  // possible to use the mode found by MCMC as
  // starting point of Minuit minimization

  m->FindMode( m->GetBestFitParameters() );
   
  // ----------------------------------
  // write out all results, plots, ...
  // ----------------------------------

  // draw all marginalized distributions into a PostScript file
  m->PrintAllMarginalized( Form( "%s/BEGE_alphas_plots.ps", OUTPUT_DIR.c_str() ) );

  // print all summary plots
//   summary->PrintParameterPlot(Form("/raid4/gerda/hemmer/BAT_BEGE_backgrounds_results_M1/%s/%dkeVBins/BEGE_backgrounds_parameters.eps",precisionString.c_str(),(int)binwidth));
//   summary->PrintCorrelationPlot(Form("/raid4/gerda/hemmer/BAT_BEGE_backgrounds_results_M1/%s/%dkeVBins/BEGE_backgrounds_correlation.eps",precisionString.c_str(),(int)binwidth));
//   summary->PrintKnowledgeUpdatePlots(Form("/raid4/gerda/hemmer/BAT_BEGE_backgrounds_results_M1/%s/%dkeVBins/BEGE_backgrounds_update.ps",precisionString.c_str(),(int)binwidth)); 
  // print results of the analysis into a text file
  m->PrintResults( Form( "%s/BEGE_alphas_results.txt", OUTPUT_DIR.c_str() ) );
  

  // calculate p-value
  double pvalue = m->EstimatePValue();
  BCLog::OutSummary(Form(" --> p-value : %.6g", pvalue));

  // -----------------------------------------------
  // write out user defined plots, information, ...
  // -----------------------------------------------

  cout << "------------------------------------------------------------" << endl;
  cout << "************************************************************" << endl;
  cout << "------------------------------------------------------------" << endl;

  for( int i = 0; i < (int)nPars; i++ )
  {
      vector<string> MCnames = m->GetMCParNames();

      //Get the parameters of interest
      BCH1D* output = m->GetMarginalized(Form("par_%d_%s",i,MCnames.at(i).c_str() ));
  
      double mode = output->GetMode();
      double xmin = 0., xmax = 0.;
      output->GetSmallestInterval(xmin, xmax);

      double quantile = output->GetQuantile(0.90);

      cout << "Parameter " << i << " " << MCnames.at(i) << ":" << endl;
      cout << "Mode: " << mode << " + " << xmax-mode << " - " << mode-xmin << endl;
      cout << "Interval: (" << xmin << " - " << xmax << ")" << endl;
      cout << "90% quantile: " << quantile << endl;
    }
  
  
  cout << "------------------------------------------------------------" << endl;
  cout << "************************************************************" << endl;
  cout << "------------------------------------------------------------" << endl;

  
//   //write additional info in the output file
//   mout->GetFile()->cd();

//   //Get the parameter of interest
//   output = m->GetMarginalized("halflife_2nbb");  
//   output->GetHistogram()->Write("histo_halflife_2nbb");

  // close output file
  mout->Close();

  // dump event information and plots for best fit parameters
  char* rootOutput = Form( "%s/BEGE_alphas_histograms.root", OUTPUT_DIR.c_str() );
  m->DumpHistosAndInfo(m->GetBestFitParameters(), rootOutput);

  delete m;
  delete summary;
  delete mout;
  
  BCLog::OutSummary("Program ran successfully");
  BCLog::OutSummary("Exiting");
  
  // close log file
  BCLog::CloseLog();
  
  return 0;

}

