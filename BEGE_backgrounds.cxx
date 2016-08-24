// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

// C/C++ includes
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>

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

// BAT includes
#include <BAT/BCMath.h>
#include <BAT/BCParameter.h>

// gerda-ada includes
#include "FileMap.h"
#include "DataLoader.h"

// own includes
#include "RunPhaseII.h"
#include "BEGE_backgrounds.h"

using namespace std;
using namespace gada;

// ---------------------------------------------------------
BEGE_backgrounds::BEGE_backgrounds() : BCModel()
{
	f_ndets = 0;
	DefineParameters();
};

// ---------------------------------------------------------
BEGE_backgrounds::BEGE_backgrounds(const char * name) : BCModel(name)
{
	f_ndets = 0;
	DefineParameters();
};

// ---------------------------------------------------------
BEGE_backgrounds::~BEGE_backgrounds()
{
	for( auto i : f_hdata ) delete i;
	f_hdata.clear();

	for( auto i : f_MC ) delete i;
	f_MC.clear();

	f_MCname.clear();
	f_vdata.clear();
	f_lowerlimits.clear();
	f_upperlimits.clear();
	f_vMC.clear();
};

// ---------------------------------------------------------
void BEGE_backgrounds::DefineParameters()
{
        vector<double> RangeMin_enrBEGe = {100.,0.,0.,0.,0.};
        vector<double> RangeMax_enrBEGe = {250.,50.,50.,100.,150.};

	// Add all MC histograms as fit parameters
	for(int iMC = 0; iMC < (int)f_MC.size(); iMC++)
	{
                AddParameter( Form( "par_%d_%s",iMC,f_MCname.at(iMC).c_str() ), RangeMin_enrBEGe.at(iMC), RangeMax_enrBEGe.at(iMC) );
	}
}

// ---------------------------------------------------------
void BEGE_backgrounds::SetHistogramParameters( int hnumbins, double hemin, double hemax )
{
	f_hnumbins = hnumbins;
	f_hemin = hemin;
	f_hemax = hemax;

	double binsize = (f_hemax - f_hemin) / (double)f_hnumbins;

	// exclude 2039 +- 25 keV from fit
	if ( (2039. + 25. - f_hemin) > 0 )
	{
		int bin_start = ceil( (2039. - 25. - f_hemin) / binsize );
		int bin_stop = ceil( (2039. + 25. - f_hemin) / binsize );

		cout << "Skipping blinded region: bins " << bin_start << " - " << bin_stop << endl;

		for( int i = bin_start; i <= bin_stop; i++)
			f_binsToSkip.push_back(i);
	}
}

// FIX ME: Detectors cannot change the data channel!!!
// FIX ME: The sum histogram is still valid but the single histograms are not!!!
// Here goes a file list with keys
// Of each run we need to load the proper data based on the detector status meta data
// The meta data file contains also information about the detector type
// Here use only enrBEGe detectors
// ---------------------------------------------------------
int BEGE_backgrounds::ReadDataEnrBEGe( std::vector<int> runlist )
{
	// Prepare sum histograms
	f_hdataSum = new TH1D("henergySum",
			"Energy channels enrBEGe ON, mult==1, no veto (mu,LAr), no TP",
			f_hnumbins, f_hemin, f_hemax);

	int bins = int(f_hemax-f_hemin);

	f_hdataSum_fine = new TH1D("henergySum_fine",
			"Energy channels enrBEGe ON, mult==1, no veto (mu,LAr), no TP",
		    bins, f_hemin, f_hemax);

	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");

/*	double liveTimeAnalyzed = 0;
	double totExposure = 0;
	double totExposureBEGe = 0;
	double totExposureCoax = 0;
*/
	for( auto Run : runlist )
	{
		cout << "---------------------------------------------" << endl;
		cout << "Adding run " << Run << " to analysis." << endl;
		cout << "---------------------------------------------" << endl;

		RunPhaseII * runX = new RunPhaseII( Run, Form( "run%04d-phy-detStatus.txt", Run ),
				Form( "run%04d-phy-analysis.txt", Run ), Form( "run%04d-phy-allFiles.txt", Run ) );

		f_ndets = runX->GetDetectors().size();
		fDetectorDynamicRange.resize( f_ndets );

		cout << "Found " << f_ndets << " detectors in Run" << Run << endl;
/*		cout << "Livetime: " << runX->GetLiveTime() << endl;
		cout << "Exposure: " << runX->GetExposure() << endl;
		cout << "enrBEGe Exposure: " << runX->GetExposureBEGE() << endl;
		cout << "enrCoax Exposure: " << runX->GetExposureCOAX() << endl;

		liveTimeAnalyzed += runX->GetLiveTime();
		totExposure += runX->GetExposure();
		totExposureBEGe += runX->GetExposureBEGE();
		totExposureCoax += runX->GetExposureCOAX();
*/
		if( f_hdata.size() == 0 )
		{
			f_hdata = vector<TH1D*>( f_ndets, NULL );
			fDetectorLiveTime = vector<double>(f_ndets, 0);
		}

		string META_FILE = runX->GetDataKeysAnalysisFile();

		cout << "Adding key list " << META_FILE << endl;

		// Here the data loader could be included...
		gada::FileMap myMap;
		myMap.SetRootDir( GERDA_PHASEII_DATA );
		myMap.BuildFromListOfKeys( META_FILE );

		gada::DataLoader l;
		l.AddFileMap(&myMap);
		l.BuildTier3();
//		l.BuildTier4();

		TChain * chain = l.GetSharedMasterChain();
		int nentries = chain->GetEntries();

		cout << "There are " << nentries << " events in the chain!" <<endl;
  
		// fill the data in histograms
		int eventChannelNumber;
		unsigned long long timestamp;
		unsigned int decimalTimestamp;
		vector<int> * firedFlag = new vector<int>(f_ndets);
		int multiplicity;
		vector<double> * energy = new vector<double>(f_ndets);
		int isTP;
//		int isMuVetoed;
//		int isLArVetoed;
		int isVetoed;
		int isVetoedInTime;
		vector<int> * failedFlag = new vector<int>(f_ndets);

		chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
		chain -> SetBranchAddress("timestamp",&timestamp);
		chain -> SetBranchAddress("decimalTimestamp",&decimalTimestamp);
		chain -> SetBranchAddress("firedFlag", &firedFlag);
		chain -> SetBranchAddress("multiplicity",&multiplicity);
		chain -> SetBranchAddress("rawEnergyGauss",&energy);
//		chain -> SetBranchAddress("energy",&energy);
		chain -> SetBranchAddress("isTP",&isTP);
//		chain -> SetBranchAddress("isMuVetoed", &isMuVetoed);
//		chain -> SetBranchAddress("isLArVetoed",&isLArVetoed);
		chain -> SetBranchAddress("isVetoed", &isVetoed);
		chain -> SetBranchAddress("isVetoedInTime", &isVetoedInTime);
		chain -> SetBranchAddress("failedFlag",&failedFlag);

		// Read only BEGe data
		vector<int> channelsBEGeON;
		vector<int> channelsOFF;

		for( auto det : runX->GetDetectors() )
		{
			string type = det->GetDetectorType();
			string status = det->GetDetectorAnalysisStatus();
			string name = det->GetDetectorName();
			int channel = det->GetDataChannel();
			
			fDetectorDynamicRange.at(channel) = det->GetDynamicRange();
			fDetectorLiveTime.at(channel) = det->GetLiveTime();

			if( type == "enrBEGe" && status == "ON" )
			{
				cout << "+ " << type << " detector " << name << " in channel " << channel << endl;

				channelsBEGeON.push_back( channel );

				if( f_hdata.at(channel) == NULL )
				{
					TH1D* henergy = new TH1D(Form("henergy_%d",channel),
							Form("Energy channel %d, mult == 1, no veto",channel),
							f_hnumbins, f_hemin, f_hemax);

					f_hdata.at(channel) = henergy;
				}
			}
			else if( status == "OFF" )
			{
				channelsOFF.push_back( channel );
			}
		}

		// loop over all events
		for (int jentry = 0; jentry < nentries; jentry++)
		{
			if ( jentry%( (int)(nentries/10) ) == 0 )
				cout << " processing event " << jentry << " (" << (int)(jentry*100/nentries) << "%)" << endl;
			chain->GetEntry( jentry );

			// Cuts after Quality cuts
			// is test pulser
			// is vetoed by muon veto
			// has multiplicity != 1
			// is vetoed by LAr veto
			if ( isTP ) 				continue;
			if ( multiplicity != 1 ) 		continue;
//			if ( isMuVetoed ) 			continue;
//			if ( isLArVetoed ) 			continue;
			if ( isVetoed ) 			continue;
			if ( isVetoedInTime ) 			continue;

			// do not consider coincidences from ANG1 and RG3 at any time
			// make coincidence loop only over other detectors
			int coincidence = 0;

			for( int i = 0; i < f_ndets; i++ )
			{
				// continue if channel is switched off and not used for AC
				if( find( channelsOFF.begin(), channelsOFF.end(), i ) != channelsOFF.end() )
					continue;
				int hitDet = firedFlag->at(i);
				if( hitDet )
					coincidence += 1;
			}

			if( coincidence != 1 ) continue;

			for( auto BEGeChannel : channelsBEGeON )
			{
				if( failedFlag->at(BEGeChannel) != 0 )		continue;
				if( firedFlag->at(BEGeChannel) != 1 )		continue;

				// do not add events with energy greater than the dynamic range of the detector
				if( energy->at(BEGeChannel) >= fDetectorDynamicRange.at( BEGeChannel ) )
					continue;

				// skip events in blinded region
				if( ( energy->at(BEGeChannel) >= (2039. - 25.) ) && ( energy->at(BEGeChannel) <= (2039. + 25.) ) )
					continue;

				// fill the channel histogram
				f_hdata.at(BEGeChannel)->Fill( energy->at(BEGeChannel) );
				// fill also the sum histogram
				f_hdataSum->Fill( energy->at(BEGeChannel) );
				f_hdataSum_fine->Fill( energy->at(BEGeChannel) );
			}
		}

		// make sure the blinded region is really empty
		for( auto BEGeChannel : channelsBEGeON )
		{
			for( auto b : f_binsToSkip )
				f_hdata.at(BEGeChannel)->SetBinContent(b,0);
		}

		for( auto b : f_binsToSkip )
			f_hdataSum->SetBinContent(b,0);

		delete runX;
	}

	// ----------------------------------------------
	// --- fill the data vector for faster access ---
	// ----------------------------------------------

	FillDataArray();

	// ----------------------------------------------
	// --- control output ---
	// ----------------------------------------------
	int channel = 0;

	cout << "---------------------------------------------" << endl;
	cout << "All selected runs combined:" << endl;
	for( auto LT : fDetectorLiveTime )
	{	
		cout << " " << channel << ": " << LT << endl;
		channel++;
	}

	cout << "---------------------------------------------" << endl;
	cout << "---------------------------------------------" << endl;

	return 0;
}

// ---------------------------------------------------------

int BEGE_backgrounds::FillDataArray()
{
  
  for(int ibin = 1; ibin <= f_hnumbins; ibin++)
  {
	  double value = f_hdataSum->GetBinContent( ibin );
	  f_vdata.push_back(value);
  }
  
  // fill also the arrays with upper and lower limits 
  // of the bins in the data histograms
  // dont't know what they should be good for...
  for(int ibin = 1; ibin <= f_hnumbins; ibin++)
  {
      double lowerlimit = f_hdataSum->GetBinLowEdge( ibin );
      double upperlimit = lowerlimit + f_hdataSum->GetBinWidth( ibin );

      f_lowerlimits.push_back(lowerlimit);
      f_upperlimits.push_back(upperlimit);
  }

  return 0;
}


// This part is a bit tricky
// Either we put the full MC files here (but they are kind of large)
// Or we treat them beforehand by smearing and histograming
// That could be difficult however based on the binning chosen for the data
// it is possible to make errors with the MC spectra
// Best is a treatment with only smearing but the full energy of the events is
// kept and histogramed only here
// ---------------------------------------------------------
int BEGE_backgrounds::ReadMCAlpha()
{
	string MC_SMEARED_DIR = getenv("MC_SMEARED_DIR");

	if( MC_SMEARED_DIR.empty() )
	{
		cout << "ERROR: environment variable MC_SMEARED_DIR not set!" << endl;
		return -1;
	}

	// Reading MC histograms with a binning of 1keV
  
	//-------------
	// Po210_pPlus
	//-------------
	f_MC_FileName = MC_SMEARED_DIR;
	f_MC_FileName += "/histograms_Po210_onPplusSurface.root";
//	AddMCSingle( "Po210_pPlus_dl0nm", "hist_dl0nm" );
//	AddMCSingle( "Po210_pPlus_dl100nm", "hist_dl100nm" );
//	AddMCSingle( "Po210_pPlus_dl200nm", "hist_dl200nm" );
//	AddMCSingle( "Po210_pPlus_dl300nm", "hist_dl300nm" );
//	AddMCSingle( "Po210_pPlus_dl400nm", "hist_dl400nm" );
//	AddMCSingle( "Po210_pPlus_dl500nm", "hist_dl500nm" );
	AddMCSingle( "Po210_pPlus_dl600nm", "hist_dl600nm" );
//	AddMCSingle( "Po210_pPlus_dl700nm", "hist_dl700nm" );
//	AddMCSingle( "Po210_pPlus_dl800nm", "hist_dl800nm" );
//	AddMCSingle( "Po210_pPlus_dl900nm", "hist_dl900nm" );
//	AddMCSingle( "Po210_pPlus_dl1000nm", "hist_dl1000nm" );

//      //------------------
//	// Ra226chain_pPlus
//	//------------------
	f_MC_FileName = MC_SMEARED_DIR;
//	f_MC_FileName += "/Ra226chain_onPplusSurface_PhaseI.root";
//	AddMCSingle("Ra226_pPlus","h_Ra226_onPplusSurface");
//	AddMCSingle("Rn222_pPlus","h_Rn222_onPplusSurface");
//	AddMCSingle("Po214_pPlus","h_Po214_onPplusSurface");
//	AddMCSingle("Po218_pPlus","h_Po218_onPplusSurface");

	f_MC_FileName += "/histograms_Ra226_onPplusSurface.root";
//	AddMCSingle("Ra226_chain_pPlus","h_Ra226chain_onPplusSurface");

//	AddMCSingle( "Ra226_pPlus_dl0nm", "hist_dl0nm" );
//	AddMCSingle( "Ra226_pPlus_dl100nm", "hist_dl100nm" );
//	AddMCSingle( "Ra226_pPlus_dl200nm", "hist_dl200nm" );
//	AddMCSingle( "Ra226_pPlus_dl300nm", "hist_dl300nm" );
//	AddMCSingle( "Ra226_pPlus_dl400nm", "hist_dl400nm" );
//	AddMCSingle( "Ra226_pPlus_dl500nm", "hist_dl500nm" );
	AddMCSingle( "Ra226_pPlus_dl600nm", "hist_dl600nm" );
//	AddMCSingle( "Ra226_pPlus_dl700nm", "hist_dl700nm" );
//	AddMCSingle( "Ra226_pPlus_dl800nm", "hist_dl800nm" );
//	AddMCSingle( "Ra226_pPlus_dl900nm", "hist_dl900nm" );
//	AddMCSingle( "Ra226_pPlus_dl1000nm", "hist_dl1000nm" );

	f_MC_FileName = MC_SMEARED_DIR;
	f_MC_FileName += "/histograms_Rn222_onPplusSurface.root";

//	AddMCSingle( "Rn222_pPlus_dl0nm", "hist_dl0nm" );
//	AddMCSingle( "Rn222_pPlus_dl100nm", "hist_dl100nm" );
//	AddMCSingle( "Rn222_pPlus_dl200nm", "hist_dl200nm" );
//	AddMCSingle( "Rn222_pPlus_dl300nm", "hist_dl300nm" );
//	AddMCSingle( "Rn222_pPlus_dl400nm", "hist_dl400nm" );
//	AddMCSingle( "Rn222_pPlus_dl500nm", "hist_dl500nm" );
	AddMCSingle( "Rn222_pPlus_dl600nm", "hist_dl600nm" );
//	AddMCSingle( "Rn222_pPlus_dl700nm", "hist_dl700nm" );
//	AddMCSingle( "Rn222_pPlus_dl800nm", "hist_dl800nm" );
//	AddMCSingle( "Rn222_pPlus_dl900nm", "hist_dl900nm" );
//	AddMCSingle( "Rn222_pPlus_dl1000nm", "hist_dl1000nm" );


//	//--------------------
//	// Ra226chain_inLArBH
//	//--------------------
	f_MC_FileName = MC_SMEARED_DIR;
	f_MC_FileName += "/Ra226chain_inLArBH_PhaseI.root";
//	AddMCSingle("Ra226_chain_inLArBH","h_Ra226chain_inLArBH");

	AddMCSingle("Ra226_inLArBH","h_Ra226_inLArBH");
//	AddMCSingle("Rn222_inLArBH","h_Rn222_inLArBH");
//	AddMCSingle("Po214_inLArBH","h_Po214_inLArBH");
//	AddMCSingle("Po218_inLArBH","h_Po218_inLArBH");


//	//------------------
//	// Ra226chain_pPlus
//	//------------------
//	// --- Ra226 ---
//	f_MC_FileName = "";
//	AddMC("Ra226_pPlus");
//
//	// --- Rn222 ---
//	f_MC_FileName = "";
//	AddMC("Rn222_pPlus");
//
//	// --- Po218 ---
//	f_MC_FileName = "";
//	AddMC("Po218_pPlus");
//
//	// --- Po214 ---
//	f_MC_FileName = "";
//	AddMC("Po214_pPlus");
//
//
//	//--------------------
//	// Ra226chain_inLArBH
//	//--------------------
//	// --- Ra226 ---
//	f_MC_FileName = "";
//	AddMC("Ra226_inLArBH");
//
//	// --- Rn222 ---
//	f_MC_FileName = "";
//	AddMC("Rn222_inLArBH");
//
//	// --- Po218 ---
//	f_MC_FileName = "";
//	AddMC("Po218_inLArBH");
//
//	// --- Po214 ---
//	f_MC_FileName = "";
//	AddMC("Po214_inLArBH");

	// copy MC information in arrays to make program faster
	FillMCArrays();

	return 0;
}

// ---------------------------------------------------------

int BEGE_backgrounds::AddMC( string name )
{
	cout << "---------------------------------------------" << endl;
	cout << "Adding MC component " << name << endl;
	cout << "---------------------------------------------" << endl;

	TH1D* henergy = new TH1D(Form("h_%s",name.c_str()),
			Form("%s",name.c_str()),
			f_hnumbins, f_hemin, f_hemax);

	int bins = int( f_hemax - f_hemin );

	string namefine = name;
	namefine.append("_fine");

	TH1D* henergy_fine = new TH1D(Form("h_%s",namefine.c_str()),
			Form("%s",namefine.c_str()),
			bins, f_hemin, f_hemax);

	int allbins = 7500;

	string nameall = name;
	nameall.append("_all");

	TH1D* henergy_all = new TH1D(Form("h_%s",nameall.c_str()),
			Form("%s",nameall.c_str()),
			allbins, 0., 7500.);

	// initialize the histogram arrays
	for( int i = 1; i <= f_hnumbins; i++ )
		henergy->SetBinContent( i, 0 );
	for( int i = 1; i <= bins; i++ )
		henergy_fine->SetBinContent( i, 0 );
	for( int i = 1; i <= allbins; i++ )
		henergy_all->SetBinContent( i, 0 );

	henergy->GetXaxis()->SetCanExtend( false );
	henergy_fine->GetXaxis()->SetCanExtend( false );
	henergy_all->GetXaxis()->SetCanExtend( false );

	// loop over all channels
	TFile * MCfile = new TFile( f_MC_FileName.c_str() );

	for( int idet = 0; idet < f_ndets; idet++ )
    {
		// Fill only MC histograms that have a corresponding data histogram
		if( f_hdata.at(idet) != NULL )
		{
			double livetime = fDetectorLiveTime.at(idet);

			cout << "+ MC spectrum h" << idet << " with lifetime weight " << livetime << endl;

			TH1D* MC_raw = (TH1D*)MCfile->Get( Form( "h%i", idet ) );
			MC_raw->Scale(livetime);

			int start_deleting = MC_raw->FindBin( (double)fDetectorDynamicRange.at( idet ) + 0.5  );
			int last_bin = MC_raw->GetNbinsX();

			for( int b = start_deleting; b <= last_bin; b++ )
			{
				MC_raw->SetBinContent( b, 0 );
			}

			for( int i = 1; i < last_bin; i++ )
			{
				double bincontent = MC_raw->GetBinContent( i );
				double bincenter = MC_raw->GetBinCenter( i );

				if( henergy->FindBin( bincenter ) > 0 && henergy->FindBin( bincenter ) <= f_hnumbins )
					henergy->SetBinContent( henergy->FindBin( bincenter ),
							henergy->GetBinContent( henergy->FindBin( bincenter ) ) + bincontent );

				if( henergy_fine->FindBin( bincenter ) > 0 && henergy_fine->FindBin( bincenter ) <= bins )
					henergy_fine->SetBinContent( henergy_fine->FindBin( bincenter ),
							henergy_fine->GetBinContent( henergy_fine->FindBin( bincenter ) ) + bincontent );

				if( henergy_all->FindBin( bincenter ) > 0 && henergy_all->FindBin( bincenter ) <= allbins )
					henergy_all->SetBinContent( henergy_all->FindBin( bincenter ),
							henergy_all->GetBinContent( henergy_all->FindBin( bincenter ) ) + bincontent );
			}
		}
	}

	MCfile->Close();

	double intInBlindedRegion = 0; 
	if( f_binsToSkip.size() > 0 )
		henergy->Integral( f_binsToSkip.front(), f_binsToSkip.back() );
	double scaling = henergy->Integral() - intInBlindedRegion;

	henergy->Scale( 1./scaling );
	henergy_fine->Scale( 1./scaling );
	henergy_all->Scale( 1./scaling );

	f_MC.push_back( henergy );
	f_MCfine.push_back( henergy_fine );
	f_MCall.push_back( henergy_all );
	f_MCname.push_back( name );

	return 0;
}


// ---------------------------------------------------------

int BEGE_backgrounds::AddMCSingle( string name, string histoname )
{
	cout << "---------------------------------------------" << endl;
	cout << "Adding MC component " << name << endl;
	cout << "---------------------------------------------" << endl;

	// set the first BEGE channel
	TH1D* henergy = new TH1D(Form("h_%s",name.c_str()),
			Form("%s",name.c_str()),
			f_hnumbins, f_hemin, f_hemax);

	int bins = int( f_hemax - f_hemin );

	string namefine = name;
	namefine.append("_fine");

	TH1D* henergy_fine = new TH1D(Form("h_%s",namefine.c_str()),
			Form("%s",namefine.c_str()),
			bins, f_hemin, f_hemax);

	int allbins = 7500;

	string nameall = name;
	nameall.append("_all");

	TH1D* henergy_all = new TH1D(Form("h_%s",nameall.c_str()),
			Form("%s",nameall.c_str()),
			allbins, 0., 7500.);

	// initialize the histogram arrays
	for( int i = 1; i <= f_hnumbins; i++ )
		henergy->SetBinContent( i, 0 );
	for( int i = 1; i <= bins; i++ )
		henergy_fine->SetBinContent( i, 0 );
	for( int i = 1; i <= allbins; i++ )
		henergy_all->SetBinContent( i, 0 );

	henergy->GetXaxis()->SetCanExtend( false );
	henergy_fine->GetXaxis()->SetCanExtend( false );
	henergy_all->GetXaxis()->SetCanExtend( false );

	// loop over all channels
	TFile * MCfile = new TFile( f_MC_FileName.c_str() );

	// get histo from file
	TH1D* MC_raw = (TH1D*)MCfile->Get( histoname.c_str() );
	int last_bin = MC_raw->GetNbinsX();

	cout << "+ MC spectrum " << histoname << endl;

	for( int i = 1; i <= last_bin; i++ )
	{
		double bincontent = MC_raw->GetBinContent( i );
		double bincenter = MC_raw->GetBinCenter( i );

		if( henergy->FindBin( bincenter ) > 0 && henergy->FindBin( bincenter ) <= f_hnumbins )
			henergy->SetBinContent( henergy->FindBin( bincenter ),
					henergy->GetBinContent( henergy->FindBin( bincenter ) ) + bincontent );

		if( henergy_fine->FindBin( bincenter ) > 0 && henergy_fine->FindBin( bincenter ) <= bins )
			henergy_fine->SetBinContent( henergy_fine->FindBin( bincenter ),
					henergy_fine->GetBinContent( henergy_fine->FindBin( bincenter ) ) + bincontent );

		if( henergy_all->FindBin( bincenter ) > 0 && henergy_all->FindBin( bincenter ) <= allbins )
			henergy_all->SetBinContent( henergy_all->FindBin( bincenter ),
					henergy_all->GetBinContent( henergy_all->FindBin( bincenter ) ) + bincontent );
	}

	cout << "MC histo loaded" << endl;

	double intInBlindedRegion = 0;
	if( f_binsToSkip.size() > 0 )
		intInBlindedRegion = henergy->Integral( f_binsToSkip.front(), f_binsToSkip.back() );
	double scaling = henergy->Integral() - intInBlindedRegion;

	cout << "Scaling MC " << 1./scaling << endl;

	henergy->Scale( 1./scaling );
	henergy_fine->Scale( 1./scaling );
	henergy_all->Scale( 1./scaling );

	f_MC.push_back( henergy );
	f_MCfine.push_back( henergy_fine );
	f_MCall.push_back( henergy_all );
	f_MCname.push_back( name );

	MCfile->Close();

	return 0;
}


// ---------------------------------------------------------

int BEGE_backgrounds::FillMCArrays()
{

  for( int iMC = 0; iMC < (int)f_MC.size(); iMC++ )
  {
      for( int ibin = 1; ibin <= f_hnumbins; ibin++ )
      {
    	  double value = f_MC.at(iMC)->GetBinContent( ibin );
    	  f_vMC.push_back(value);
      }
  }

  return 0;
}

// ---------------------------------------------------------
    
double BEGE_backgrounds::LogLikelihood(const std::vector <double> & parameters)
{

  // This methods returns the logarithm of the conditional probability
  // p(data|parameters). This is where you have to define your model.
  //
  // For this method to not be alpha specific we have to exclude the blinding region from the likelihood

  double logprob = 0.;

  for( int ibin = 0; ibin < f_hnumbins; ibin++)
  {
	  // skip blinded region in the likelihood
	  if( find( f_binsToSkip.begin(), f_binsToSkip.end(), ibin ) != f_binsToSkip.end() )
		  continue;

	  double lambda = 0.;

	  for( int iMC = 0; iMC < (int)f_MC.size(); iMC++)
	  {
		  lambda += parameters.at(iMC) * f_vMC[iMC*f_hnumbins + ibin];
	  }
	    
      double bincontent = f_vdata[ibin];
      
      double sum = bincontent*log(lambda) - lambda - BCMath::LogFact((int)bincontent);
      
      logprob += sum;
  }

  return logprob;
}

// ---------------------------------------------------------
double BEGE_backgrounds::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

  double logprob = 0.;

  // skip blinded region in the likelihood
  int bins = f_binsToSkip.size();
  double binsize = (f_hemax - f_hemin) / (double)f_hnumbins;
  double skippedRange = binsize * (double)bins;

  for( int iMC = 0; iMC < (int)f_MC.size(); iMC++ )
  {
      double range = GetParameter(iMC)->GetRangeWidth();
      // flat prior for all contributions
      logprob += log( 1. / ( range - skippedRange ) );
  }
  
  return logprob;
}
// ---------------------------------------------------------

double BEGE_backgrounds::EstimatePValue()
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

  vector<double> mean( f_hnumbins, 0 );
  vector<int> nom( f_hnumbins, 0 );

  vector<double> parameters = GetBestFitParameters();

  for( int ibin = 0; ibin < f_hnumbins; ibin++ )
  {
	  //TODO
	  if( find( f_binsToSkip.begin(), f_binsToSkip.end(), ibin ) != f_binsToSkip.end() )
		  continue;

      double lambda = 0.;

      for( int iMC = 0; iMC < (int)f_MC.size(); iMC++ )
      {
    	  lambda += parameters.at(iMC)*f_vMC[iMC*f_hnumbins+ibin];
      }

      int counter = ibin;

      mean.at(counter) = std::max( lambda, 1e-8 );
      nom.at(counter) = int( mean.at(counter) );
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
      for( int ibin = 0; ibin < f_hnumbins; ibin++)
      {
    	  if( find( f_binsToSkip.begin(), f_binsToSkip.end(), ibin ) != f_binsToSkip.end() )
    		  continue;

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

  cout << "p-value is " << pvalue << endl;

  return pvalue;
}

// ---------------------------------------------------------

void BEGE_backgrounds::DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename)
{
  TFile* rootOut = new TFile( rootfilename, "recreate" );

  if(!rootOut->IsOpen()) cout<<"No rootfile opened!"<<endl;

  rootOut->cd();

  TH1D* hMC = new TH1D("hMC", "model", f_hnumbins, f_hemin, f_hemax);

  int bins = int( f_hemax - f_hemin );
  TH1D* hMC_fine = new TH1D("hMC_fine", "model", bins, f_hemin, f_hemax);

  int binsall = 7500;
  TH1D* hMC_all = new TH1D("hMC_all", "model", binsall, 0, 7500.);

  TH1D* hresiduals = new TH1D("hresiduals", "residuals", f_hnumbins, f_hemin, f_hemax);

  vector<double> eventsMC;
  vector<double> eventsMC_fine;
  vector<double> eventsMC_all;

  for( int iMC = 0; iMC < (int)f_MC.size(); iMC++ )
  {
	  // prepare the histograms
      f_MC.at(iMC)->Scale( parameters.at(iMC) );
      f_MCfine.at(iMC)->Scale( parameters.at(iMC) );
      f_MCall.at(iMC)->Scale( parameters.at(iMC) );

      f_MC.at(iMC)->Write();
      f_MCfine.at(iMC)->Write();
      f_MCall.at(iMC)->Write();

      hMC->Add( f_MC.at(iMC) );
      hMC_fine->Add( f_MCfine.at(iMC) );
      hMC_all->Add( f_MCall.at(iMC) );

      // count the events
      eventsMC.push_back( f_MC.at(iMC)->Integral() );
      eventsMC_all.push_back( f_MCall.at(iMC)->Integral() );
  }

  // write the single detector data spectra
  for(int idet=0; idet<f_ndets; idet++)
  {
	  if( f_hdata.at(idet) != NULL )
	  {
		  double binwidth = f_hdata.at(idet)->GetBinWidth(1);

		  f_hdata.at(idet)->SetLineWidth(2);
		  f_hdata.at(idet)->GetXaxis()->SetTitle("Energy (keV)");
		  f_hdata.at(idet)->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth));
		  f_hdata.at(idet)->Write();
	  }
  }
  // write the finer binning data histogram
  f_hdataSum_fine->SetLineWidth(2);
  f_hdataSum_fine->GetXaxis()->SetTitle("Energy (keV)");
  f_hdataSum_fine->GetYaxis()->SetTitle("Events/keV");
  f_hdataSum_fine->Write();

  // write the finer binning MC histogram
  hMC_fine->SetLineWidth(2);
  hMC_fine->GetXaxis()->SetTitle("Energy (keV)");
  hMC_fine->GetYaxis()->SetTitle("Events/keV");
  hMC_fine->Write();

  // write the complete spectrum MC histogram
  hMC_all->SetLineWidth(2);
  hMC_all->GetXaxis()->SetTitle("Energy (keV)");
  hMC_all->GetYaxis()->SetTitle("Events/keV");
  hMC_all->Write();

  // write the fit results

  TCanvas* canvas=new TCanvas("canvas", "");
  canvas->cd();

  TLegend* legend = new TLegend(0.73, 0.62, 0.98, 0.87);
  legend->SetTextFont(62);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  double binwidth_hsum = f_hdataSum->GetBinWidth(1);

  f_hdataSum->SetLineWidth(2);
  f_hdataSum->SetMarkerStyle(20);
  f_hdataSum->GetXaxis()->SetTitle("Energy (keV)");
  f_hdataSum->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth_hsum));
  f_hdataSum->Draw("P");
  f_hdataSum->Write();
  legend->AddEntry(f_hdataSum,"data","p");

  double binwidth_hMC = hMC->GetBinWidth(1);

  hMC->SetLineWidth(2);
  hMC->SetLineColor(1);
  hMC->Draw("same");
  hMC->GetXaxis()->SetTitle("Energy (keV)");
  hMC->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth_hMC));
  hMC->Write();
  legend->AddEntry(hMC,"model","l");

  //the single contributions

  string name;
  
  for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
  {
      f_MC.at(iMC)->SetLineWidth(2);

      if(iMC < 10) f_MC.at(iMC)->SetLineColor(iMC+2);

      if(iMC >= 8)
      {
    	  f_MC.at(iMC)->SetLineColor(iMC-8+2);
    	  f_MC.at(iMC)->SetLineStyle(2);
      }

      f_MC.at(iMC)->Draw("same");
      f_MC.at(iMC)->GetXaxis()->SetTitle("Energy (keV)");
      f_MC.at(iMC)->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth_hMC));
      name=f_MCname.at(iMC);
      legend->AddEntry(f_MC.at(iMC),Form("%s",name.c_str()),"l");

      f_MC.at(iMC)->Write();
  }

  legend->Draw();
  canvas->Write();
      
  hresiduals->Add(f_hdataSum);
  hresiduals->Divide(hMC); //ratio
  hresiduals->SetLineWidth(2);
  hresiduals->SetMarkerStyle(20);
  hresiduals->GetXaxis()->SetTitle("Energy (keV)");
  hresiduals->GetYaxis()->SetTitle(Form("residual counts/(%d keV)",(int)binwidth_hMC));
  hresiduals->Write();

  // events info
  cout << "---------------------------------------------" << endl;
  cout << "---------------------------------------------" << endl << endl;
  for( int iMC = 0; iMC < (int)eventsMC.size(); iMC++ )
  {
      cout << "Events: " << endl;
      cout << f_MCname.at(iMC) << ": " << eventsMC.at(iMC) << endl;
      cout << "In total spectrum from (0 - 7500) keV" << endl;
      cout << f_MCname.at(iMC) << ": " << eventsMC_all.at(iMC) << endl;
  }
  cout << "Total events in MC: " << hMC->Integral() << endl;
  cout << endl;
  cout << "Events in data: " << f_hdataSum->Integral() << endl;
  cout << endl;
  cout << "Estimate of BI: " << hMC_all->GetBinLowEdge( 1401 ) << "keV-" << hMC_all->GetBinLowEdge( 1501 ) << "keV" << endl;
  cout << hMC_all->Integral( 1401, 1500 ) << "cts -> " << hMC_all->Integral( 1401, 1500 )/100./5.988 << "cts/(keV kg yr)" << endl;
  cout << "Estimate of BI: " << hMC_all->GetBinLowEdge( 1415 ) << "keV-" << hMC_all->GetBinLowEdge( 1465 ) << "keV" << endl;
  cout << hMC_all->Integral( 1415, 1464 ) << "cts -> " << hMC_all->Integral( 1415, 1464 )/50./5.988 << "cts/(keV kg yr)" << endl;
  cout << "---------------------------------------------" << endl;
  cout << "---------------------------------------------" << endl << endl;

  rootOut->Close();
}
// ---------------------------------------------------------

