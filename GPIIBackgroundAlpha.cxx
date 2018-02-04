// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

// C/C++ includes
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>

// jsoncpp
#include "json/value.h"
#include "json/json.h"

// ROOT includes
#include "TChain.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TColor.h"
#include "TF1.h"
#include "TLegend.h"
#include "TParameter.h"
#include "TObject.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TLine.h"

// BAT includes
#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

#include "GETRunConfiguration.hh"
#include "GERunConfigurationManager.hh"
#include "GADatasetManager.h"

#include "ProgressBar.h"

// own includes
#include "GPIIBackgroundAlpha.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
GPIIBackgroundAlpha::GPIIBackgroundAlpha() : BCModel()
{
    f_verbosity = 0;
    f_obin = 0;
};

// ---------------------------------------------------------
GPIIBackgroundAlpha::GPIIBackgroundAlpha(const char * name) : BCModel(name)
{
    f_verbosity = 0;
    f_obin = 0;
};

// ---------------------------------------------------------
GPIIBackgroundAlpha::GPIIBackgroundAlpha(string masterconfname) : BCModel()
{
    f_obin = 0;
    SetMasterConf( masterconfname );
    UnwrapMasterConf();
    DefineParameters();
};

// ---------------------------------------------------------
GPIIBackgroundAlpha::~GPIIBackgroundAlpha()
{
    for( auto i : f_hdata ) delete i.second;
    f_hdata.clear();

    for( auto i : f_MC ) delete i;
    f_MC.clear();

    f_vdata.clear();
    f_lowerlimits.clear();
    f_upperlimits.clear();
    f_vMC.clear();
};

// ---------------------------------------------------------
void GPIIBackgroundAlpha::SetMasterConf( string masterconfname )
{
    f_masterconfname = masterconfname;

    f_j_masterconf  = GetJsonValueFromFile( f_masterconfname );
    SetParConf( f_j_masterconf["parconf"].asString() );
    SetDetConf( f_j_masterconf["detconf"].asString() );
    SetRunConf( f_j_masterconf["runconf"].asString() );
    SetEnvConf( f_j_masterconf["envconf"].asString() );

    return;
}

// ---------------------------------------------------------
void GPIIBackgroundAlpha::UnwrapMasterConf()
{
    SetVerbosity( f_j_masterconf["verbosity"].asInt() );
    SetPrecision( f_j_masterconf["precision"].asString() );
    if( FitOverflowBin() ) f_obin = 1;
    SetHistogramParameters(
        f_j_masterconf["histo"]["binning"].asDouble(),
        f_j_masterconf["histo"]["min"].asDouble(),
        f_j_masterconf["histo"]["max"].asDouble()
    );

    f_fitbins = f_hnumbins + f_obin;

    MCMCSetFlagFillHistograms( f_j_masterconf["MCMC-fill-histograms"].asBool() );

    return;
}

// ---------------------------------------------------------
void GPIIBackgroundAlpha::SetPrecision( string precisionString )
{
    BCEngineMCMC::Precision precision = BCEngineMCMC::kLow;
    if( precisionString == "kLow" ) precision = BCEngineMCMC::kLow;
    else if( precisionString == "kMedium" ) precision = BCEngineMCMC::kMedium;
    else if( precisionString == "kHigh" ) precision = BCEngineMCMC::kHigh;
    else if( precisionString == "kVeryHigh" ) precision = BCEngineMCMC::kVeryHigh;

    MCMCSetPrecision(precision);

    return;
}

void GPIIBackgroundAlpha::SetHistogramParameters(double hBinning, double hMin, double hMax)
{
    int hNbins = (int) ( (hMax - hMin) / hBinning );
    hMax = hMin + hNbins * hBinning;

    SetHistogramParameters(hNbins, hMin, hMax);

    BCLog::OutSummary( Form( "Fit Range: %.0f - %.0f keV", hMin, hMax) );
    BCLog::OutSummary( Form( "Binning: %.0f keV", hBinning ) );
    BCLog::OutSummary( Form( "Number of Bins: %i", hNbins ) );

    return;
}

// ---------------------------------------------------------
void GPIIBackgroundAlpha::DefineParameters()
{
    if( f_verbosity > 0 ) cout << "Defining parameters" << endl;

    int nParametersSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++ )
    {
        // skip parameter if requested
        if( SkipParameter(p) )
        {
            if( f_verbosity > 0 ) cout << "Parameter " << p << " skipped " << endl;
            nParametersSkipped++;
            continue;
        }

        string name = ParameterName(p);
        double min = f_j_parconf["parameters"][p]["min"].asDouble();
        double max = f_j_parconf["parameters"][p]["max"].asDouble();
        int nbins = f_j_parconf["parameters"][p]["nbins"].asInt();

        // add parameter, set range and binning
        AddParameter( name.c_str(), min, max );
        GetParameter( name.c_str() )->SetNbins( nbins );
        SetPriorConstant( p-nParametersSkipped );

        if( f_verbosity > 0 )
        {
            cout << "\t" << name << ": " << "[" << nbins << "|" << min << ":" << max << "]" << endl;
        }
    }

    return;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::InitializeDataHistograms()
{
    string name = "hSum";
    string name_fine = name; name_fine += "_fine";
    string name_all = name; name_all += "_all";
    int bins = (int)f_hemax-f_hemin;

    f_hdataSum = new TH1D( name.c_str(), name.c_str(), f_hnumbins, f_hemin, f_hemax);
    f_hdataSum_fine = new TH1D( name_fine.c_str(), name_fine.c_str(), bins, f_hemin, f_hemax);
    f_hdataSum_all = new TH1D( name_all.c_str(), name_all.c_str(), 7500, 0., 7500.);

    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();

        string name_single = Form( "hSingle_%s", det.c_str() );

        TH1D * henergy = new TH1D( name_single.c_str(), name_single.c_str(), 7500, 0., 7500.);

        f_hdata[det] = henergy;
        f_DetectorLiveTime[det] = 0.;
    }

    BCLog::OutSummary( "Data histograms initialized" );

    return 0;
}

// WRAPPER
// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadData()
{
    InitializeDataHistograms();

    int stat;

    string filename = f_j_masterconf["reprocess-data"]["filename"].asString();

    if( f_j_masterconf["reprocess-data"]["force"].asBool() )
        stat = ReadDataFromEvents( filename );
    else
        stat = ReadDataFromHistogram( filename );

    FillDataArray();

    return stat;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadDataFromHistogram( string infilename )
{
    // read * detector live times (from json file)
    //      * run livetimes (from json file),
    //      * histogram from file (from root file) with binning possibly changed
    TFile * file = new TFile( infilename.c_str(), "READ" );

    if( !file->IsOpen() )
    {
        BCLog::OutSummary( Form( "File not found: %s", infilename.c_str()) );
        BCLog::OutSummary( "Reading data from events" );
        int stat = ReadDataFromEvents( infilename );

    return stat;
    }

    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();

        TH1D * hh = (TH1D*) file->Get( Form("hSingle_%s", det.c_str()) );

        if( f_verbosity > 0 )
            cout << "Histo integral " << det << ": " << hh->Integral() << endl;

        f_hdata.at(det) -> Add( hh );
    }

    TH1D * hsum = (TH1D*)file->Get("hSum_all");

    if( f_verbosity > 0 )
        cout << "Histo integral sum: " <<  hsum -> Integral() << endl;

    f_hdataSum_all->Add( hsum );

    int nbins = f_hdataSum_all->GetNbinsX();

    for( int b = 1; b <= nbins; b++ )
    {
        double bincontent = f_hdataSum_all->GetBinContent( b );
        double bincenter = f_hdataSum_all->GetBinCenter( b );

        f_hdataSum -> Fill( bincenter, bincontent );
        f_hdataSum_fine -> Fill( bincenter, bincontent );
    }

    file->Close();

    if( f_verbosity > 0 )
    {
        cout << "Data histo integral: " << f_hdataSum -> Integral() << endl;
        cout << "Data fine histo integral: " << f_hdataSum_fine -> Integral() << endl;
        cout << "Data all histo integral: " << f_hdataSum_all -> Integral() << endl;
    }

    BCLog::OutSummary( Form( "Data histograms read from file: %s", infilename.c_str() ) );

    ifstream runLTconf( GetRunConfLTName() );
    Json::Value j_RunLiveTime;
    runLTconf >> j_RunLiveTime;
    runLTconf.close();

    for( unsigned int r = 0; r < f_nruns; r++ )
    {
        string run = to_string( f_j_runconf["runs"][r].asInt() );
        double rLT = j_RunLiveTime[run].asDouble();

        f_RunLiveTime.push_back(rLT);
    }

    BCLog::OutSummary( Form( "Run Livetimes read from file: %s", GetRunConfLTName().c_str() ) );

    ifstream detLTconf( GetDetConfLTName() );
    Json::Value j_DetectorLiveTime;
    detLTconf >> j_DetectorLiveTime;
    detLTconf.close();

    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();
        double dLT = j_DetectorLiveTime[det].asDouble();

        f_DetectorLiveTime[det] = dLT;
    }

    BCLog::OutSummary( Form( "Detector Livetimes read from file: %s", GetDetConfLTName().c_str() ) );

    return 0;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadDataFromEvents( string outfilename )
{
    // write * detector live times (to json file)
    //       * run livetimes (to json file),
    //       * histograms from file (to root file)

    // read analysis key lists of each run
    string GERDA_DATA_SETS = f_j_envconf["GERDA_DATA_SETS"].asString();
    vector<string> keylist;

    if( f_verbosity > 0 ) cout << "Reading run data: " << endl;

    for( unsigned int r = 0; r < f_nruns; r++ )
    {
        int run = f_j_runconf["runs"][r].asInt();

        string keylist = GERDA_DATA_SETS;
        keylist += "/run"; keylist += Form( "%04d", run );
        keylist += "-phy-analysis.txt";

        if( f_verbosity > 0 ) cout << "\t" << keylist << endl;

        ReadRunData( keylist );
    }

    if( f_verbosity > 0 )
    {
        cout << "Detector LiveTimes: " << endl;

        for( unsigned int d = 0; d < f_ndets; d++ )
        {
            string det = f_j_detconf["detectors"][d].asString();
            cout << "\t" << det << ": " << f_DetectorLiveTime[det] << endl;
        }
    }

    WriteDataToFileForFastAccess( outfilename );

    return 0;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::WriteDataToFileForFastAccess( string outfilename )
{
    TFile * outfile = new TFile( outfilename.c_str(), "RECREATE" );

    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();
        f_hdata[det]->Write();
    }

    f_hdataSum->Write();
    f_hdataSum_fine->Write();
    f_hdataSum_all->Write();

    outfile->Close();

    BCLog::OutSummary( Form( "Data histograms written to file: %s", outfilename.c_str() ) );

    // Write to run livetimes to json files
    Json::Value j_RunLiveTime;

    for( unsigned int r = 0; r < f_nruns; r++ )
    {
        int run = f_j_runconf["runs"][r].asInt();

        j_RunLiveTime[ to_string(run) ] = f_RunLiveTime[r];
    }

    ofstream runLTconf( GetRunConfLTName() );
    runLTconf << j_RunLiveTime;
    runLTconf.close();

    BCLog::OutSummary( Form( "Run Livetimes written to file: %s", GetRunConfLTName().c_str() ) );

    // Write to detector livetimes to json files
    Json::Value j_DetectorLiveTime;

    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();

        j_DetectorLiveTime[ det ] = f_DetectorLiveTime[det];
    }

    ofstream detLTconf( GetDetConfLTName() );
    detLTconf << j_DetectorLiveTime;
    detLTconf.close();

    BCLog::OutSummary( Form( "Detector Livetimes written to file: %s", GetDetConfLTName().c_str() ) );

    return 0;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadRunData( string keylist )
{
    string GERDA_PHASEII_DATA = f_j_envconf["GERDA_PHASEII_DATA"].asString();
    string MU_CAL = f_j_envconf["MU_CAL"].asString();

    if( GERDA_PHASEII_DATA.empty() )
    {
        cout << "Environment variable GERDA_PHASEII_DATA not set" << endl;
        exit(EXIT_FAILURE);
    }
    if( MU_CAL.empty() )
    {
        cout << "Environment variable MU_CAL not set" << endl;
        cout << "LNGS: /nfs/gerda5/gerda-data/blind/active/meta/config/_aux/geruncfg" << endl;
        cout << "MPIK: /lfs/l3/gerda/Daq/data-phaseII/blind/active/meta/config/_aux/geruncfg" << endl;
        exit(EXIT_FAILURE);
    }

    // initialize run configuration manager
    GERunConfigurationManager * RunConfManager = new GERunConfigurationManager();
    GETRunConfiguration * RunConf = 0;
    RunConfManager -> AllowRunConfigurationSwitch(true);
    RunConfManager -> SetVerbosity(1);

    GADatasetManager * theDatasetManager = new GADatasetManager();

    // Here the data loader could be included...
    gada::FileMap myMap;
    myMap.SetRootDir( GERDA_PHASEII_DATA );
    myMap.BuildFromListOfKeys( keylist );

    gada::DataLoader l;
    l.AddFileMap(&myMap);
    l.BuildTier3();
    bool hasTier4 = l.BuildTier4();

    TChain * chain = l.GetSharedMasterChain();
    int nentries = chain->GetEntries();

    cout << "Entries in Run: " << nentries << endl;

    // set chain branches
    // TIER3
    int eventChannelNumber;
    int multiplicity;
    int isTP;
    int isVetoedInTime;
    unsigned long long timestamp;
    vector<Double_t> * energy = new vector<double>(f_ndets);
    vector<Int_t> * firedFlag = new vector<int>(f_ndets);
    vector<Int_t> * failedFlag = new vector<int>(f_ndets);
    vector<Int_t> * failedFlag_isPhysical = new vector<int>(f_ndets);
    vector<Int_t> * failedFlag_isSaturated = new vector<int>(f_ndets);
    // TIER4
    bool hasPsdIsEval = true;
    vector<Int_t> * psdIsEval = new vector<int>(f_ndets);

    chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
    chain -> SetBranchAddress("multiplicity",&multiplicity);
    chain -> SetBranchAddress("isTP",&isTP);
    chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
    chain -> SetBranchAddress("timestamp",&timestamp);
    chain -> SetBranchAddress("rawEnergyGauss",&energy);
    chain -> SetBranchAddress("firedFlag", &firedFlag);
    chain -> SetBranchAddress("failedFlag",&failedFlag);
    chain -> SetBranchAddress("failedFlag_isPhysical",&failedFlag_isPhysical);
    chain -> SetBranchAddress("failedFlag_isSaturated",&failedFlag_isSaturated);
    if( hasTier4 )
    {
        if ( chain->GetBranchStatus("psdIsEval") )
        {
            chain->SetBranchAddress("psdIsEval",&psdIsEval);
        }
        else
        {
            hasPsdIsEval = false;
            cout << "WARNING: psdIsEval branch not found" << endl;
        }
    }
    else
    {
        cout << "WARNING: No tier4 file found. Assuming all BEGe data valid for PSD" << endl;
    }

    int nTP = 0;
    double frequencyTP = 1./20.;

    ProgressBar bar( nentries, '#', false );

    // loop over all events
    for (int e = 0; e < nentries; e++)
    {
        bar.Update();

        chain->GetEntry( e );

        RunConf = RunConfManager -> GetRunConfiguration( timestamp );

        if ( isTP )
        {
            nTP++;
            for( unsigned int d = 0; d < f_ndets; d++ )
            {
                string det = f_j_detconf["detectors"][d].asString();
                if( IsOn( RunConf, det ) ) f_DetectorLiveTime[det] += 1./(frequencyTP*60.*60.*24.);
            }
            continue;
        }
        if ( multiplicity > 1 ) continue;
        if ( isVetoedInTime )   continue;

        for( unsigned int d = 0; d < f_ndets; d++  )
        {
            string det = f_j_detconf["detectors"][d].asString();
            if( !IsOn( RunConf, det ) ) continue;

            // 0,"EnrBEGe"
            // 1,"EnrCoax"
            // 2,"Runs47-49"
            // 3,"Natural"
            // 4,"EnrBEGeNoPSD"
            int idataset  = theDatasetManager->FindDataset( timestamp, d, hasPsdIsEval ? psdIsEval->at(d) : false );
            if( idataset  < 0 ) continue; // not in any data set
            if( idataset == 4 ) continue; // in EnrBEGeNoPSD data set

            int c = GetChannel( RunConf, det );

            double en = energy->at(c);

            // fill energy spectra
            if( ( multiplicity == 1 && !failedFlag_isPhysical->at(c) ) ||
                !failedFlag_isSaturated->at(c) )
            {
                if( !failedFlag_isSaturated->at(c) ) en = 10000.;

                f_hdata[det] -> Fill( en );
                f_hdataSum -> Fill( en );
                f_hdataSum_fine -> Fill( en );
                f_hdataSum_all -> Fill( en );
            }
        }
    }

    double runLiveTimeInDays = nTP /( frequencyTP * 60. * 60. * 24. );

    // run live time in days
    f_RunLiveTime.push_back( runLiveTimeInDays );

    if( f_verbosity > 0 ) cout << "Run Livetime: " << runLiveTimeInDays << endl;

    return 0;
}

// ---------------------------------------------------------

int GPIIBackgroundAlpha::FillDataArray()
{
    for(int ibin = 1; ibin <= f_fitbins; ibin++)
    {
        double value = f_hdataSum->GetBinContent( ibin );
        f_vdata.push_back(value);
    }

    // fill also the arrays with upper and lower limits
    // of the bins in the data histograms
    // dont't know what they should be good for...
    for(int ibin = 1; ibin <= f_fitbins; ibin++)
    {
        double lowerlimit = f_hdataSum->GetBinLowEdge( ibin );
        double upperlimit = lowerlimit + f_hdataSum->GetBinWidth( ibin );

        f_lowerlimits.push_back(lowerlimit);
        f_upperlimits.push_back(upperlimit);
    }

    BCLog::OutSummary( "Data arrays filled" );

    return 0;
}


int GPIIBackgroundAlpha::InitializeMCHistograms()
{
//    int nParametersSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++ )
    {
//        if( SkipParameter(p) ) nParametersSkipped++; continue;

        string pname = ParameterName(p);
        int ncorrelations = NumberOfCorrelatedHistos(p);

        for( int c = 0; c < ncorrelations; c++ )
        {
            string hname = f_j_parconf["parameters"][p]["mc"][c]["histoname"].asString();
            string fname = f_j_parconf["parameters"][p]["mc"][c]["filename"].asString();

            string name = pname; name += "_n"; name += to_string(c);
            string name_fine = name; name_fine += "_fine";
            string name_all = name; name_all += "_all";
            string title = fname; title += " : "; title += hname;
            int bins = (int)f_hemax-f_hemin;

            TH1D * hmc = new TH1D( name.c_str(), title.c_str(), f_hnumbins, f_hemin, f_hemax);
            TH1D * hmc_fine = new TH1D( name_fine.c_str(), title.c_str(), bins, f_hemin, f_hemax);
            TH1D * hmc_all = new TH1D( name_all.c_str(), title.c_str(), 7500, 0., 7500.);

            f_MC.push_back( hmc );
            f_MC_fine.push_back( hmc_fine );
            f_MC_all.push_back( hmc_all );
        }
    }

    BCLog::OutSummary( "MC histograms initialized"  );

    return 0;
}


// FIX ME: LOW
// for now only single histogram MC pdfs can be read
// would be nice to construct them from single detector histograms scaled by detector live time
// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadMC()
{
    InitializeMCHistograms();

    BCLog::OutSummary( "Read MC pdfs" );

    string GERDA_MC_PDFS = f_j_envconf["GERDA_MC_PDFS"].asString();

    if( GERDA_MC_PDFS.empty() )
    {
        cout << "ERROR: environment variable GERDA_MC_PDFS not set!" << endl;
        return -1;
    }

    // loop over parameters
    int index = 0;
    int nParametersSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++ )
    {
        // skip parameter if requested
        if( SkipParameter(p) ) { nParametersSkipped++; continue; }

        int ncorr = NumberOfCorrelatedHistos(p);

        for( int c = 0; c < ncorr; c++ )
        {
            string histoname = f_j_parconf["parameters"][p]["mc"][c]["histoname"].asString();
            string filename = GERDA_MC_PDFS; filename += "/";
            filename += f_j_parconf["parameters"][p]["mc"][c]["filename"].asString();

            ReadSingleMC( p, c, index, histoname, filename );

            index++;
        }
    }
/*
    // Reading MC histograms with a binning of 1keV
    //-------------
    // Po210_pPlus
    //-------------
    f_MC_FileName += "/histograms_Po210_onPplusSurface.root";
    f_MC_FileName += "/histograms_Ra226_onPplusSurface.root";
    f_MC_FileName += "/histograms_Rn222_onPplusSurface.root";
    f_MC_FileName += "/Ra226chain_inLArBH_PhaseI.root";
//    AddMCSingle("Ra226_pPlus","h_Ra226_onPplusSurface");
//    AddMCSingle("Rn222_pPlus","h_Rn222_onPplusSurface");
//    AddMCSingle("Po214_pPlus","h_Po214_onPplusSurface");
//    AddMCSingle("Po218_pPlus","h_Po218_onPplusSurface");
*/
    FillMCArrays();

    return 0;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::ReadSingleMC( int par_index, int local_index, int global_index, string histoname, string filename )
{
    BCLog::OutSummary( "--------------------------------------------" );
    BCLog::OutSummary(  Form( "parameter: %d", par_index  ) );
    BCLog::OutSummary(  Form( "local index: %d", local_index ) );
    BCLog::OutSummary(  Form( "global index: %d", global_index ) );
    BCLog::OutSummary(  Form( "File: %s", filename.c_str() ) );
    BCLog::OutSummary(  Form( "Histogram: %s", histoname.c_str() ) );

    // open file
    TFile * MCfile = new TFile( filename.c_str() );

    // get histogram [7500|0:7500]
    TH1D * MChisto = (TH1D*) MCfile->Get( histoname.c_str() );

    int nbins = MChisto->GetNbinsX();

    for( int b = 1; b <= nbins; b++ )
    {
        double bincontent = MChisto->GetBinContent( b );
        double bincenter = MChisto->GetBinCenter( b );

        f_MC[global_index]      -> Fill( bincenter, bincontent );
        f_MC_fine[global_index] -> Fill( bincenter, bincontent );
        f_MC_all[global_index]  -> Fill( bincenter, bincontent );
    }

    MCfile->Close();

    // scale histograms with number of primaries
    double primaries = GeneratedPrimaries(par_index,local_index);

    BCLog::OutSummary( Form( "Primaries: %.0f", primaries ) );

    f_MC[global_index]      -> Scale(1./primaries);
    f_MC_fine[global_index] -> Scale(1./primaries);
    f_MC_all[global_index]  -> Scale(1./primaries);

    if( f_verbosity > 0 )
    {
        cout << "fMC scaled integral: " << f_MC[global_index] -> Integral() << endl;
        cout << "fMC_fine scaled integral: " << f_MC_fine[global_index] -> Integral() << endl;
        cout << "fMC_all scaled integral: " << f_MC_all[global_index] -> Integral() << endl;
    }

    BCLog::OutSummary( "--------------------------------------------" );

    return 0;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::FillMCArrays()
{
    int nMC = f_MC.size();

      for( int iMC = 0; iMC < nMC; iMC++ )
      {
        for( int ibin = 1; ibin <= f_fitbins; ibin++ )
          {
              double value = f_MC[iMC]->GetBinContent( ibin );
              f_vMC.push_back(value);
          }
      }

    BCLog::OutSummary( "MC arrays filled" );

      return 0;
}

// ---------------------------------------------------------
double GPIIBackgroundAlpha::LogLikelihood(const vector <double> & parameters)
{
    // This methods returns the logarithm of the conditional probability
    // p(data|parameters). This is where you have to define your model.
    //
    // For this method to not be alpha specific we have to exclude the blinding region from the likelihood

    double logprob = 0.;

    for( int ibin = 0; ibin < f_fitbins; ibin++)
    {
        double lambda = 0.;
        int nHistosRead = 0;
        int nParametersSkipped = 0;

        for( unsigned int p = 0; p < f_npars; p++ )
        {
            // skip parameter if requested
            if( SkipParameter(p) ) { nParametersSkipped++; continue; }

            int index = p-nParametersSkipped;

            int ncorrelations = NumberOfCorrelatedHistos(p);

            for( int c = 0; c < ncorrelations; c++ )
            {
                double weight = WeightOfHistogram(p,c);
                lambda += parameters[index] * weight * f_vMC[ nHistosRead * f_fitbins + ibin ];
                nHistosRead++;
            }
      }

      int bincontent = (int)f_vdata[ibin];

      double sum = bincontent*log(lambda) - lambda - BCMath::LogFact(bincontent);

      logprob += sum;
  }

  return logprob;
}

/*
// ONLY NEEDED IF THE PRIOR PROBABILITY IS MORE COMPLICATED
// ---------------------------------------------------------
double GPIIBackgroundAlpha::LogAPrioriProbability(const vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

  double logprob = 0.;

  for( unsigned int p = 0; p < f_npars; p++ )
  {
      // skip parameter if requested
      if( SkipParameter(p) ) continue;


      double range = GetParameter(p-nParametersSkipped)->GetRangeWidth();
      // flat prior for all contributions
      logprob += log( 1. / range );
  }

  return logprob;
}
*/

// ---------------------------------------------------------
double GPIIBackgroundAlpha::EstimatePValue()
{
    //Allen's routine for the evaluation of p-value
    //This is derived from PRD 83 (2011) 012004, appendix
    //taken from Luciano
    double logp0 = LogLikelihood( GetBestFitParameters() ); //likelihood at the mode

    /*
      Now we initialize the Markov Chain by setting the number of entries in a bin to
      the integer
      part of the mean in that bin.  This will give the maximum possible likelihood.
    */
    double sumlog = 0;


    /* mean is the array where you have the expected mean in each bin (calculated
       from the parameter values at the mode. Nom is the nominal value of entries in
       each bin (integer part)
    */

    vector<double> mean( f_fitbins, 0 );
    vector<int> nom( f_fitbins, 0 );

    vector<double> parameters = GetBestFitParameters();

    for( int ibin = 0; ibin < f_fitbins; ibin++ )
    {
        double lambda = 0.;
        int nHistosRead = 0;
        int nParametersSkipped = 0;

        for( unsigned int p = 0; p < f_npars; p++)
        {
            // skip parameter if requested
            if( SkipParameter(p) ) { nParametersSkipped++; continue; }

            int index = p-nParametersSkipped;

            int ncorrelations = NumberOfCorrelatedHistos(p);

            for( int c = 0; c < ncorrelations; c++ )
            {
                double weight = WeightOfHistogram(p,c);
                lambda += parameters[index] * weight * f_vMC[ nHistosRead * f_fitbins + ibin ];
                nHistosRead++;
            }
        }

        int counter = ibin;

        mean[counter] = max( lambda, 1e-8 );
        nom[counter] = int( mean[counter] );
        sumlog += BCMath::LogPoisson( nom[counter], mean[counter] );
    }

    cout << "Logprob for best: " << sumlog << endl;

    /*
    Now we run the Markov chain to generate new configurations.  One iteration means
    a loop over all the bins, with an attempt to vary each bin up or down by one unit.  We
    accept/reject at each step  and compare the data logprob to the simulated at the end of each iteration.
    */
    const int nloops = 100000;
    int Pgood = 0;

    for( int iloop = 0; iloop < nloops; iloop++ )
    {
        for( int ibin = 0; ibin < f_fitbins; ibin++)
        {
            int counter = ibin;

            if ( rand() > RAND_MAX/2 ) // Try to increase the bin content by 1
            {
                double r = mean[counter]/(nom[counter]+1);
                double rtest = double(rand())/RAND_MAX;
                if( rtest < r ) //Accept
                {
                    nom[counter] = nom[counter]+1;
                    sumlog += log(r);
                }
            }
            else // Try to decrease the bin content by 1
            {
                double r = nom[counter]/mean[counter];
                double rtest = double(rand())/RAND_MAX;
                if ( rtest < r ) //Accept
                {
                    nom[counter] = nom[counter]-1;
                    sumlog += log(r);
                }
            }
        }
        if ( sumlog < logp0 ) Pgood++;
    }

    double pvalue = double(Pgood)/double(nloops);

    if( f_verbosity > 0) cout << "p-value is " << pvalue << endl;

    return pvalue;
}

// ---------------------------------------------------------
// FIX ME overflowbin plot
void GPIIBackgroundAlpha::DumpHistosAndInfo( string rootfilename )
{
    TFile* rootOut = new TFile( rootfilename.c_str(), "RECREATE" );

    if(!rootOut->IsOpen()) cout<<"No rootfile opened!"<<endl;

    rootOut->cd();

    const vector<double> parameters = GetBestFitParameters();
    const vector<double> parameters_error = GetBestFitParameterErrors();

    int nParametersSkipped = 0;

    BCLog::OutSummary( "Best fit parameters:" );
    for( unsigned int p = 0; p < f_npars; p++ )
    {
        if( SkipParameter(p) ){ nParametersSkipped++; continue; }

        string pname = ParameterName(p);
        int i = p-nParametersSkipped;

        BCLog::OutSummary( Form("Par%d: %.2f+-%.2f (%s)",
            i, parameters[i], parameters_error[i], pname.c_str()) );
    }

    int bins = int( f_hemax - f_hemin );
    TH1D* hMC = new TH1D("hMC", "hMC", f_hnumbins, f_hemin, f_hemax);
    TH1D* hMC_fine = new TH1D("hMC_fine", "hMC_fine", bins, f_hemin, f_hemax);
    TH1D* hMC_all = new TH1D("hMC_all", "model", 7500, 0, 7500.);

    TH1D* hresiduals = new TH1D("hresiduals", "residuals", f_hnumbins, f_hemin, f_hemax);

    vector<TH1D*> p_MC;
    vector<TH1D*> p_MC_fine;
    vector<TH1D*> p_MC_all;

    vector<double> eventsMC;
    vector<double> eventsMC_all;

    int nHistosRead = 0;
    nParametersSkipped = 0;

    // prepare singel MC pdfs
    for( unsigned int p = 0; p < f_npars; p++ )
    {
        // skip parameter if requested
        if( SkipParameter(p) ) { nParametersSkipped++; continue; }

        int index = p-nParametersSkipped;

        double scale = parameters[index];

        int ncorrelations = NumberOfCorrelatedHistos(p);

        for( int c = 0; c < ncorrelations; c++ )
        {
            double weight = WeightOfHistogram(p,c);

            int iMC = nHistosRead+c;

            f_MC[iMC]       -> Scale( scale*weight );
            f_MC_fine[iMC]   -> Scale( scale*weight );
            f_MC_all[iMC]    -> Scale( scale*weight );

            hMC         -> Add( f_MC[iMC] );
            hMC_fine    -> Add( f_MC_fine[iMC] );
            hMC_all     -> Add( f_MC_all[iMC] );

            eventsMC.push_back( f_MC[iMC]->Integral() );
            eventsMC_all.push_back( f_MC_all[iMC]->Integral() );

            if( c == 0 )
            {
                p_MC.push_back( f_MC[iMC] );
                p_MC_fine.push_back( f_MC_fine[iMC] );
                p_MC_all.push_back( f_MC_all[iMC] );
            }
            else
            {
                p_MC.back()->Add( f_MC[iMC] );
                p_MC_fine.back()->Add( f_MC_fine[iMC] );
                p_MC_all.back()->Add( f_MC_all[iMC] );
            }
        }

        nHistosRead += ncorrelations;
    }

    cout << "Parameters: " << f_npars << endl;
    cout << "Of which skipped: " << nParametersSkipped << endl;
    cout << "Combined MC pdfs: " << p_MC.size() << endl;

    if( (f_npars-nParametersSkipped) != p_MC.size() )
    {
        cout << "Something went wrong in preparing the MC histograms" << endl;
        exit(EXIT_FAILURE);
    }

    int binning = (int)f_j_masterconf["histo"]["binning"].asDouble();

    // write the single detector data spectra
    for( unsigned int d = 0; d < f_ndets; d++ )
    {
        string det = f_j_detconf["detectors"][d].asString();

        f_hdata.at(det)->SetLineWidth(2);
        f_hdata.at(det)->GetXaxis()->SetTitle("energy (keV)");
        f_hdata.at(det)->GetYaxis()->SetTitle(Form("cts/(%d keV)",binning));
        f_hdata.at(det)->Write();
    }

    // write the finer binning data histogram
    f_hdataSum_fine->SetLineWidth(2);
    f_hdataSum_fine->GetXaxis()->SetTitle("energy (keV)");
    f_hdataSum_fine->GetYaxis()->SetTitle("cts/keV");
    f_hdataSum_fine->Write();

    // write the finer binning data histogram
    f_hdataSum_all->SetLineWidth(2);
    f_hdataSum_all->GetXaxis()->SetTitle("energy (keV)");
    f_hdataSum_all->GetYaxis()->SetTitle("cts/keV");
    f_hdataSum_all->Write();

    // write the finer binning MC histogram
    hMC_fine->SetLineWidth(2);
    hMC_fine->GetXaxis()->SetTitle("energy (keV)");
    hMC_fine->GetYaxis()->SetTitle("cts/keV");
    hMC_fine->Write();

    // write the complete spectrum MC histogram
    hMC_all->SetLineWidth(2);
    hMC_all->GetXaxis()->SetTitle("energy (keV)");
    hMC_all->GetYaxis()->SetTitle("cts/keV");
    hMC_all->Write();

    // write the fit results to a canvas
    TCanvas* canvas = new TCanvas("fitresult", "");
    canvas->Divide(1,2);
    canvas->cd(1);

    TLegend* legend = new TLegend(0.73, 0.62, 0.98, 0.87);
    legend->SetTextFont(62);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);

    f_hdataSum->SetLineWidth(2);
    f_hdataSum->SetMarkerStyle(20);
    f_hdataSum->GetXaxis()->SetTitle("energy (keV)");
    f_hdataSum->GetYaxis()->SetTitle(Form("cts/(%d keV)",binning));
    f_hdataSum->Draw("P");
    f_hdataSum->Write();
    legend->AddEntry(f_hdataSum,"data","p");

    hMC->SetLineWidth(2);
    hMC->SetLineColor(1);
    hMC->Draw("samehisto");
    hMC->GetXaxis()->SetTitle("Energy (keV)");
    hMC->GetYaxis()->SetTitle(Form("Events/(%d keV)",binning));
    hMC->Write();
    legend->AddEntry(hMC,"model","l");

    //the single contributions
    string name;
    int index = 0;
    nParametersSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++)
    {
        // skip parameter if requested
        if( SkipParameter(p) ) { nParametersSkipped++; continue; }

        name = ParameterName(p);

        index = p - nParametersSkipped;

        p_MC.at(index)->SetLineWidth(2);

        int pp = index/8;
        p_MC.at(index)->SetLineColor( index%8 + 2 );
        p_MC.at(index)->SetLineStyle( pp + 1 );

        p_MC.at(index)->Draw("samehisto");
        p_MC.at(index)->GetXaxis()->SetTitle("energy (keV)");
        p_MC.at(index)->GetYaxis()->SetTitle( Form("cts/(%d keV)",binning) );
        legend->AddEntry(p_MC.at(index),Form("%s",name.c_str()),"l");

    // fix me
        p_MC.at(index)->Write(Form("%s",name.c_str()));
    }

    legend->Draw();

    canvas->cd(2);
    hresiduals->Add(f_hdataSum);
    hresiduals->Divide(hMC); //ratio
    hresiduals->SetLineWidth(2);
    hresiduals->SetMarkerStyle(20);
    hresiduals->GetXaxis()->SetTitle("energy (keV)");
    hresiduals->GetYaxis()->SetTitle(Form("residual counts/(%d keV)",binning));
    hresiduals->Draw();
    hresiduals->Write();

    TLine * l = new TLine( f_hemin, 1., f_hemax, 1. );
    l->SetLineWidth(1);
    l->SetLineStyle(2);
    l->Draw();

    canvas->Write();

    // events info
    cout << "---------------------------------------------" << endl;
    cout << "---------------------------------------------" << endl;
    nHistosRead = 0;
    nParametersSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++ )
    {
        // skip parameter if requested
        if( SkipParameter(p) ) { nParametersSkipped++; continue; }

        int ncorrelations = NumberOfCorrelatedHistos(p);

        cout << "Par: " << ParameterName(p) << endl;

        for( int c = 0; c < ncorrelations; c++ )
        {
            cout << "Hist: " << f_j_parconf["parameters"][p]["mc"][c]["histoname"].asString() << endl;
            cout << "\tevents (" << f_hemin << " - " << f_hemax << "): " << eventsMC[nHistosRead] << endl;
            cout << "\tevents (0 - 7500): " << eventsMC_all[nHistosRead++] << endl;
        }
    }
    cout << "Total events in MC: " << hMC->Integral() << endl;
    cout << "Events in data: " << f_hdataSum->Integral() << endl;
//    cout << "Estimate of BI: " << hMC_all->GetBinLowEdge( 1401 ) << "keV-" << hMC_all->GetBinLowEdge( 1501 ) << "keV" << endl;
//    cout << hMC_all->Integral( 1401, 1500 ) << "cts -> " << hMC_all->Integral( 1401, 1500 )/100./5.988 << "cts/(keV kg yr)" << endl;
//    cout << "Estimate of BI: " << hMC_all->GetBinLowEdge( 1415 ) << "keV-" << hMC_all->GetBinLowEdge( 1465 ) << "keV" << endl;
//    cout << hMC_all->Integral( 1415, 1464 ) << "cts -> " << hMC_all->Integral( 1415, 1464 )/50./5.988 << "cts/(keV kg yr)" << endl;
    cout << "---------------------------------------------" << endl;
    cout << "---------------------------------------------" << endl;

    rootOut->Close();

    if( f_j_masterconf["update-parameters"].asBool() )
        UpdateParameters( "updated-parameters.json" );

    return;
}

int GPIIBackgroundAlpha::UpdateParameters( string filename, double nerror )
{
    const vector<double> parameters = GetBestFitParameters();
    const vector<double> parameters_error = GetBestFitParameterErrors();

    Json::Value updatedpars = f_j_parconf;

    int nSkipped = 0;

    for( unsigned int p = 0; p < f_npars; p++ )
    {
        if(SkipParameter(p)) { nSkipped++; continue; }

        int index = p - nSkipped;

        double newmin = max( 0., parameters[index] - nerror * parameters_error[index] );
        double newmax = parameters[index] + nerror * parameters_error[index];

        updatedpars["parameters"][p]["min"] = newmin;
        updatedpars["parameters"][p]["max"] = newmax;
    }

    ofstream file( filename.c_str() );
    file << updatedpars;
    file.close();

    cout << "Updated Parameters written to file: " << filename << endl;

    return 0;
}


// ---------------------------------------------------------
Json::Value GPIIBackgroundAlpha::GetJsonValueFromFile( string filename )
{
    ifstream file( filename, ifstream::binary );

    if( !file.is_open() )
    {
        cout << "Cannot open file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    if( f_verbosity > 0 ) cout << "Reading JSON file:" << filename << endl;

    Json::Value val;

    file >> val;

    file.close();

    if( f_verbosity > 0 ) cout << val << "\n";

    return val;
}

// ---------------------------------------------------------
bool GPIIBackgroundAlpha::IsOn( GETRunConfiguration * RunConf, string det )
{
    int c = GetChannel( RunConf, det );

        if( RunConf -> IsOn( c ) )
        return true;

    return false;
}

// ---------------------------------------------------------
int GPIIBackgroundAlpha::GetChannel( GETRunConfiguration * RunConf, string det )
{
    GEChannel * channel = RunConf -> GetChannel( det.c_str() );
        int c = channel -> GetChannelNumber();

    return c;
}

// dir structure
// OUTDIR/parconf/detconf/runconf/precision/
// ---------------------------------------------------------
string GPIIBackgroundAlpha::GetOutputDirectory()
{
    string parconf = f_j_masterconf["parconf"].asString();
    int from = parconf.find_last_of('/')+1;
    int length = parconf.find_last_of('-') - from;
    string dir_parconf = parconf.substr( from, length );

    string detconf = f_j_masterconf["detconf"].asString();
    from = detconf.find_last_of('/')+1;
    length = detconf.find_last_of('-') - from;
    string dir_detconf = detconf.substr( from, length );

    string runconf = f_j_masterconf["runconf"].asString();
    from = runconf.find_last_of('/')+1;
    length = runconf.find_last_of('-') - from;
    string dir_runconf = runconf.substr( from, length );

    string precision = f_j_masterconf["precision"].asString();

    string OUTDIR = f_j_masterconf["output-directory"].asString();
    OUTDIR += "/";
    OUTDIR += dir_parconf;
    OUTDIR += "/";
    OUTDIR += dir_detconf;
    OUTDIR += "/";
    OUTDIR += dir_runconf;
    OUTDIR += "/";
    OUTDIR += precision;

    return OUTDIR;
}

// ---------------------------------------------------------
string GPIIBackgroundAlpha::GetOutputFilenameBase()
{
    string parconf = f_j_masterconf["parconf"].asString();
    int from = parconf.find_last_of('/')+1;
    int length = parconf.find_last_of('-') - from;
    string dir_parconf = parconf.substr( from, length );

    string detconf = f_j_masterconf["detconf"].asString();
    from = detconf.find_last_of('/')+1;
    length = detconf.find_last_of('-') - from;
    string dir_detconf = detconf.substr( from, length );

    string runconf = f_j_masterconf["runconf"].asString();
    from = runconf.find_last_of('/')+1;
    length = runconf.find_last_of('-') - from;
    string dir_runconf = runconf.substr( from, length );

    string precision = f_j_masterconf["precision"].asString();

    double binning = f_j_masterconf["histo"]["binning"].asDouble();

    bool fitoverflow = f_j_masterconf["fitoverflow"].asBool();
    string s_fitoverflow = "overflowbin";

    if( fitoverflow )   s_fitoverflow = "yof";
    else                s_fitoverflow = "nof";

    string filename_base = "BATOutput_";
    filename_base += dir_parconf;
    filename_base += "_";
    filename_base += dir_detconf;
    filename_base += "_";
    filename_base += dir_runconf;
    filename_base += "_";
    filename_base += precision;
    filename_base += "_";
    filename_base += to_string( (int)binning );
    filename_base += "_";
    filename_base += to_string( (int)f_hemin );
    filename_base += "_";
    filename_base += to_string( (int)f_hemax );
    filename_base += "_";
    filename_base += s_fitoverflow;

    return filename_base;
}

string GPIIBackgroundAlpha::GetRunConfLTName()
{
    string runconfname = f_j_masterconf["runconf"].asString();
    string detconfname = f_j_masterconf["detconf"].asString();

    string runLTconfname = "./data/";

    int from = runconfname.find_last_of('/') + 1;
    int to = runconfname.find_last_of('-');
    runLTconfname += runconfname.substr( from, to-from );
    runLTconfname += "-";
    from = detconfname.find_last_of('/') + 1;
    to = detconfname.find_last_of('-');
    runLTconfname += detconfname.substr( from, to-from );

    runLTconfname += "-RunLT.json";

    return runLTconfname;
}

string GPIIBackgroundAlpha::GetDetConfLTName()
{
    string runconfname = f_j_masterconf["runconf"].asString();
    string detconfname = f_j_masterconf["detconf"].asString();

    string detLTconfname = "./data/";

    int from = runconfname.find_last_of('/') + 1;
    int to = runconfname.find_last_of('-');
    detLTconfname += runconfname.substr( from, to-from );
    detLTconfname += "-";
    from = detconfname.find_last_of('/') + 1;
    to = detconfname.find_last_of('-');
    detLTconfname += detconfname.substr( from, to-from );

    detLTconfname += "-DetLT.json";

    return detLTconfname;
}

bool GPIIBackgroundAlpha::SkipParameter( int p )
{
    bool use = f_j_parconf["parameters"][p].get("use",true).asBool();
    return !use;
}

string GPIIBackgroundAlpha::ParameterName( int p )
{
    string name = f_j_parconf["parameters"][p]["name"].asString();
    return name;
}

int GPIIBackgroundAlpha::NumberOfCorrelatedHistos( int p )
{
    int n = f_j_parconf["parameters"][p]["mc"].size();

    return n;
}

double GPIIBackgroundAlpha::WeightOfHistogram( int p, int c )
{
    double weight = f_j_parconf["parameters"][p]["mc"][c].get("weight",1.0).asDouble();

    return weight;
}

double GPIIBackgroundAlpha::GeneratedPrimaries( int p, int c )
{
    double primaries = f_j_parconf["parameters"][p]["mc"][c]["primaries"].asDouble();

    return primaries;
}

bool GPIIBackgroundAlpha::FitOverflowBin()
{
    bool f = f_j_masterconf["fitoverflow"].asBool();

    if(f) BCLog::OutSummary( "Fitting WITH OVERFLOW bin" );

    return f;
}
