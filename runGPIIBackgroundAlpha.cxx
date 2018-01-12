// ***************************************************************
// This file was created using the CreateProject.sh script
// for project 2nbb.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <getopt.h>

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
    string masterfileJSON = "./config/default-masterconf.json";

    int choice = 0;

    static struct option long_options[] =
    {
        { "help",       no_argument,       0,   'h' },
        { "master",     required_argument, 0,   'Z' }, // master json config file
        { 0, 0, 0, 0 }
    };

    int option_index = 0;

    while( (choice = getopt_long(argc, argv, "hZ:", long_options, &option_index)) != -1 )
    {
        switch (choice)
        {
            case 'h':
                UsageGPIIBackgroundAlpha();
                return 0;
            case 'Z':
                masterfileJSON = optarg;
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
    GPIIBackgroundAlpha * m = new GPIIBackgroundAlpha( masterfileJSON );

    // create a new summary tool object
    BCSummaryTool * summary = new BCSummaryTool(m);

    BCLog::OutSummary("Model created");

    // initialize data histograms and read data
    m->ReadData();

    // initialize MC histograms and read MC pdfs
    m->ReadMC();

    string OUTPUT_DIR = m->GetOutputDirectory();
    string OUT_FILE_BASE = OUTPUT_DIR; OUT_FILE_BASE += "/"; OUT_FILE_BASE += m->GetOutputFilenameBase();

    // create directory if it does not exist
    string cmd = "mkdir -p "; cmd += OUTPUT_DIR;
    system( cmd.c_str() );

    string mout_filename = OUT_FILE_BASE; mout_filename += "_model.root";
    BCModelOutput* mout = new BCModelOutput( m, mout_filename.c_str() );

    // write MCMC chain if requested
    mout->WriteMarkovChain( m->GetWriteMCMCChain() );

    // normalize distributions
    m->SetIntegrationMethod(BCIntegrate::kIntDefault);
    m->Normalize();

    // Marginalize all
    m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // find mode using Minuit using MCMC mode as starting point
    m->FindMode( m->GetBestFitParameters() );

    // ----------------------------------
    // write results
    // ----------------------------------

    string updatefilename = OUT_FILE_BASE; updatefilename += "_update.pdf";
	summary -> PrintKnowledgeUpdatePlots( updatefilename.c_str() );

    string paramfilename = OUT_FILE_BASE; paramfilename += "_param.pdf";
	summary -> PrintParameterPlot( paramfilename.c_str() );

    string corrfilename = OUT_FILE_BASE; corrfilename += "_correlations.pdf";
	summary -> PrintCorrelationPlot( corrfilename.c_str() );

    string corrmatrixfilename = OUT_FILE_BASE; corrmatrixfilename += "_corrmatrix.pdf";
	summary -> PrintCorrelationMatrix( corrmatrixfilename.c_str() );

//    string marg_filename = OUT_FILE_BASE; marg_filename += "_plots.pdf";
//    m->PrintAllMarginalized( marg_filename.c_str() );

    string resu_filename = OUT_FILE_BASE; resu_filename += "_results.txt";
    m->PrintResults( resu_filename.c_str() );

//    m->DumpHistosAndInfo( m->GetBestFitParameters() );

    // Caluculate p-Value
    double pvalue = m->EstimatePValue();
    BCLog::OutSummary( Form(" --> p-value : %.6g", pvalue) );

    BCLog::OutSummary("Program ran successfully");
    BCLog::OutSummary("Exiting");

    delete m;
    delete summary;
    delete mout;

    // close log file
    BCLog::CloseLog();

    return 0;
}


void UsageGPIIBackgroundAlpha()
{
    cout << "Usage: runGPIIBackgroundAlpha -Z masterconf.json" << endl;
    cout << "\tOptions" << endl;
    cout << "\t\t-h : print this help menu" << endl;
    cout << "\t\t-Z : set master json config file (./config/default-masterconf.json)" << endl;
    return;
}
