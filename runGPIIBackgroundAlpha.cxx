// ***************************************************************
// This file was created using the CreateProject.sh script
// for project 2nbb.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

// BAT includes
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCEngineMCMC.h>
#include <BAT/BCH1D.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCParameter.h>

// root includes
#include "TStyle.h"

// own includes
#include "GPIIBackgroundAlpha.h"

using namespace std;

void UsageGPIIBackgroundAlpha();

int main( int argc, char* argv[] )
{
    // arguments
    string data_set = "enrBEGe";
    string precisionString = "kLow";
    string runlist;
    string detectorlist = "";
    bool useDetectorList = false;
    string fitselection;

    double hMin = 3500., hMax = 7500.;  // fit range in keV
    double hBinning = 30.;              // bin size in keV

    int choice = 0;

    static struct option long_options[] =
    {
        { "runlist",    required_argument, 0,   'R' },
        { "help",       no_argument,       0,   'h' },
        { "detectors",  required_argument, 0,   'D' },
        { "precision",  required_argument, 0,   'P' },
        { "dataset",    required_argument, 0,   'S'},
        { "binning",    required_argument, 0,   'B' },
        { "min",        required_argument, 0,   'm' },
        { "max",        required_argument, 0,   'M' },
        { "fit",        required_argument, 0,   'F' },
        { 0, 0, 0, 0 }
    };

    int option_index = 0;

    while( (choice = getopt_long(argc, argv, "hR:D:P:S:B:m:M:F:", long_options, &option_index)) != -1 )
    {
        switch (choice)
        {
            case 'h':
                UsageGPIIBackgroundAlpha();
                return 0;
            case 'R':
                runlist = optarg;
                break;
            case 'D':
                detectorlist = optarg;
                useDetectorList = true;
                break;
            case 'P':
                precisionString = optarg;
                break;
            case 'S':
                data_set = optarg;
                useDetectorList = false;
                break;
            case 'B':
                hBinning = atof( optarg );
                break;
            case 'm':
                hMin = atof( optarg );
                break;
            case 'M':
                hMax = atof( optarg );
                break;
            case 'F':
                fitselection = optarg;
                break;
            default:
                cout << "Unknown option: -" << choice << endl;
                UsageGPIIBackgroundAlpha();
                return -1;
        }
    }

    // --------- set input parameters ------

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();
    gStyle->SetOptStat(0);

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // create new GPIIBackgroundAlpha object
    GPIIBackgroundAlpha* m = new GPIIBackgroundAlpha();

    // create a new summary tool object
    BCSummaryTool * summary = new BCSummaryTool(m);

    BCLog::OutSummary("Model created");

    // Set the precision for the fit
    BCEngineMCMC::Precision precision;
    if( precisionString == "kLow" ) precision = BCEngineMCMC::kLow;
    else if( precisionString == "kMedium" ) precision = BCEngineMCMC::kMedium;
    else if( precisionString == "kHigh" ) precision = BCEngineMCMC::;kHigh
    else if( precisionString == "kVeryHigh" ) precision = BCEngineMCMC::kVeryHigh;
    m->MCMCSetPrecision(precision);

    BCLog::OutSummary( Form( "Precision: %s ", precisionString.c_str() ) );

    // Set the histograms for the data spectra
    // If hMax does not fit the bin width then lower hMax accordingly
    int hNbins = (int) ( (hMin - hMax) / hBinning );
    hMax = hMin + hNbins * hBinning;

    m->SetHistogramParameters(hNbins, hMin, hMax);

    BCLog::OutSummary( Form( "Fit Range: %.0f - %.0f keV", hMin, hMax) );
    BCLog::OutSummary( Form( "Binning: %.0f keV", hBinning ) );
    BCLog::OutSummary( Form( "Number of Bins: %i", hNbins ) );

    // read in the names of the files you want to use,
    m->ReadData( runlist, data_set, detectorlist, useDetectorList );

/*
    // --- read in the MC histograms ---
    m->ReadMCAlpha();

  // --- define parameters ---
  m->DefineParameters();

  // to get better sensitivity on a certain parameter
//   m->SetNbins("inverseHalflife_GPIIBackgroundAlpha",2000);
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
//   summary->PrintParameterPlot(Form("/raid4/gerda/hemmer/BAT_GPIIBackgroundAlpha_results_M1/%s/%dkeVBins/GPIIBackgroundAlpha_parameters.eps",precisionString.c_str(),(int)binwidth));
//   summary->PrintCorrelationPlot(Form("/raid4/gerda/hemmer/BAT_GPIIBackgroundAlpha_results_M1/%s/%dkeVBins/GPIIBackgroundAlpha_correlation.eps",precisionString.c_str(),(int)binwidth));
//   summary->PrintKnowledgeUpdatePlots(Form("/raid4/gerda/hemmer/BAT_GPIIBackgroundAlpha_results_M1/%s/%dkeVBins/GPIIBackgroundAlpha_update.ps",precisionString.c_str(),(int)binwidth));
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
*/
  return 0;

}



void UsageGPIIBackgroundAlpha()
{
    cout << "Usage: runGPIIBackgroundAlpha" << endl;
    cout << "\tOptions" << endl;
    cout << "\t\t-h : this help menu" << endl;
    cout << "\t\t-R : run selection from file" << endl;
    cout << "\t\t-D : detector selection from file" << endl;
    cout << "\t\t-P : precision [kLow(default)|kMedium|kHigh|kVeryHigh]" << endl;
    cout << "\t\t-S : dataset [enrBEGe(default)|enrCoax|natCoax]" << endl;
    cout << "\t\t-B : binning in keV" << endl;
    cout << "\t\t-m : lower histogram limit in keV" << endl;
    cout << "\t\t-M : upper histogram limit in keV" << endl;
    cout << "\t\t-F : pdf selection to fit [Po210]" << endl;
    return;
}
