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
#include "BEGE_backgrounds.h"

using namespace std;

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
	// Add all MC histograms as fit parameters
	for(int iMC = 0; iMC < (int)f_MC.size(); iMC++)
	{
		AddParameter( Form("par_%d",iMC), 0., 100. );
	}
}

// ---------------------------------------------------------
void BEGE_backgrounds::SetHistogramParameters( int hnumbins, double hemin, double hemax )
{
	f_hnumbins = hnumbins;
	f_hemin = hemin;
	f_hemax = hemax;
}

// FIX ME: Detectors cannot change the data channel!!!
// FIX ME: The sum histogram is still valid but the single histograms are not!!!
// Here goes a file list with keys
// Of each run we need to load the proper data based on the detector status meta data
// The meta data file contains also information about the detector type
// Here use only enrBEGe detectors
// ---------------------------------------------------------
int BEGE_backgrounds::ReadData( string runlist )
{
	// Prepare sum histograms
	f_hdataSum = new TH1D("henergySum",
			"Energy channels enrBEGe ON, mult==1, no veto",
			f_hnumbins, f_hemin, f_hemax);

	int bins = int(f_hemax-f_hemin);

	f_hdataSum_fine = new TH1D("henergySum_fine",
			"Energy channels enrBEGe ON, mult==1, no veto",
		    bins, f_hemin, f_hemax);

	string GERDA_PHASEII_DATA = getenv("GERDA_PHASEII_DATA");

	ifstream RunList( runlist );

	if( !RunList.is_open() )
	{
		cout << "No run list found." << endl;
		return 1;
	}

	int Run; RunList >> Run;

	while( !RunList.eof() )
	{
		RunPhaseII * runX = new RunPhaseII( Run, Form( "run%04d-phy-detStatus.txt", Run ),
				Form( "run%04d-phy-analysis.txt", Run ), Form( "run%04d-phy-allFiles.txt", Run ) );

		f_ndets = runX->GetDetectors().size();

		if( f_hdata.size() == 0 )
		{
			f_hdata = vector<TH1D*>( f_ndets, NULL );
		}

		string META_FILE = runX->GetDataKeysAnalysisFile();

		// Here the data loader could be included...
		gada::FileMap myMap;
		myMap.SetRootDir( GERDA_PHASEII_DATA );
		myMap.BuildFromListOfKeys( META_FILE );

		gada::DataLoader l;
		l.AddFileMap(&myMap);
		l.BuildTier3();

		TChain * chain = l.GetMasterChain();
		int nentries = chain->GetEntries();

		cout << "There are " << nentries << " events in the chain!" <<endl;
  
		// fill the data in histograms
		int eventChannelNumber;
		int isvetoed;
		int isvetoedInTime;
		int isTP;
		int firedChannels;
		double ch[f_ndets];
		int firedFlag[f_ndets];

		chain -> SetBranchAddress("eventChannelNumber", &eventChannelNumber);
		chain -> SetBranchAddress("isVetoed",&isvetoed);
		chain -> SetBranchAddress("isVetoedInTime",&isvetoedInTime);
		chain -> SetBranchAddress("firedFlag", firedFlag);
		chain -> SetBranchAddress("multiplicity",&firedChannels);
		chain -> SetBranchAddress("energy",ch);
		chain -> SetBranchAddress("isTP",&isTP);

		// Modify all this part in order to read only BEGe data
		// Maybe using the ChannelMapping file we also use in MaGe
		vector<int> channelsBEGeON;
		vector<int> channelsOFF;

		for( auto det : runX->GetDetectors() )
		{
			string type = det->GetDetectorType();
			string status = det->GetDetectorAnalysisStatus();
			int channel = det->GetDataChannel();

			if( type == "enrBEGe" && status == "ON" )
			{
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

		delete runX;

		// loop over all events
		for (int jentry = 0; jentry < nentries; jentry++)
		{
			//if (jentry%100000==0) cout<<" now event "<<jentry<<endl;
			chain->GetEntry(jentry);

			if (isTP) continue;
			if (firedChannels != 1) continue;
			if (isvetoed) continue;
			if (isvetoedInTime) continue;

			// do not consider coincidences from ANG1 and RG3 at any time
			// make coincidence loop only over other detectors
			// --> channel 0==ANG1, channel 7==RG3
			int coincidence = 0;
			for( int i = 0; i < f_ndets; i++)
			{
				// continue if channel is switched off and not used for AC
				if( find( channelsOFF.begin(), channelsOFF.end(), i ) != channelsOFF.end() )
					continue;
				int hitDet = firedFlag[i];
				if( hitDet == 1 )
					coincidence+=1;
			}

			if( coincidence != 1 ) continue;

			for( auto BEGeChannel : channelsBEGeON )
			{
				if( firedFlag[BEGeChannel] == 1 )
				{
					// fill the channel histogram
					f_hdata.at(BEGeChannel)->Fill( ch[BEGeChannel] );
					// fill also the sum histogram
					f_hdataSum->Fill( ch[BEGeChannel] );
					f_hdataSum_fine->Fill( ch[BEGeChannel] );
				}
			}
		}
	}

	// ----------------------------------------------
	// --- fill the data vector for faster access ---
	// ----------------------------------------------

	FillDataArray();

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
// Or we treat them beforehand by smearing and histogramming
// That could be difficult however based on the binning chosen for the data
// it is possible to make errors with the MC spectra
// Best is a treatment with only smearing but the full energy of the events is
// kept and histogrammed only here
// ---------------------------------------------------------
int BEGE_backgrounds::ReadMC()
{
	// Reading MC histograms with a binning of 1keV
  
	//-------------
	// Po210_pPlus
	//-------------
	f_chain = new TChain("smearedTree");
	for(int i=0; i<100; i++)
    {
		f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/"
				"Po210_pPlus/Po210_onPplusContactSurface_GERDAPhaseI_1000000evts_%04d_smeared.root",i));
    }
	AddMC("Po210_pPlus");
	f_chain->Clear();

  //------------------
  // Ra226chain_pPlus
  //------------------
  // --- Ra226 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_pPlus/"
    		  "Ra226_onPplusContactSurface_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Ra226_pPlus");
  f_chain->Clear();

  // --- Rn222 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_pPlus/"
    		  "Rn222_onPplusContactSurface_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Rn222_pPlus");
  f_chain->Clear();

  // --- Po218 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_pPlus/"
    		  "Po218_onPplusContactSurface_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Po218_pPlus");
  f_chain->Clear();

  // --- Po214 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_pPlus/"
    		  "Po214_onPplusContactSurface_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Po214_pPlus");
  f_chain->Clear();


  //--------------------
  // Ra226chain_inLArBH
  //--------------------
  // --- Ra226 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_inLArBH/Ra226_inANG3BoreHole_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Ra226_inLArBH");
  f_chain->Clear();

  // --- Rn222 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_inLArBH/Rn222_inANG3BoreHole_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Rn222_inLArBH");
  f_chain->Clear();

  // --- Po218 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_inLArBH/Po218_inANG3BoreHole_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Po218_inLArBH");
  f_chain->Clear();

  // --- Po214 ---
  f_chain = new TChain("smearedTree");
  for(int i=0; i<10; i++)
    {
      f_chain->AddFile(Form("/raid4/gerda/hemmer/BEGE_backgrounds/MC_files/Ra226chain_inLArBH/Po214_inANG3BoreHole_GERDAPhaseI_1000000evts_000%d_smeared.root",i));
    }
  AddMC("Po214_inLArBH");
  f_chain->Clear();


  // copy MC information in arrays to make program faster
  FillMCArrays();

  return 0;
}

// ---------------------------------------------------------
int BEGE_backgrounds::AddMC(string name)
{

  int nentries = f_chain->GetEntries();

  cout << "There are " << nentries << " events in the chain!" << endl;
  
  // fill the data in histograms

  const int nchannels_tot = 40;

  int firedChannels;
  float det_energy[nchannels_tot];
  int firedFlag[nchannels_tot];

  f_chain->SetBranchAddress("firedFlag", firedFlag);
  f_chain->SetBranchAddress("multiplicity",&firedChannels);
  f_chain->SetBranchAddress("det_energy",det_energy);


  // set the first BEGE channel
  
  TH1D* henergy = new TH1D(Form("h_%s",name.c_str()),
			   Form("%s",name.c_str()),
			   f_hnumbins, f_hemin, f_hemax);

  int bins=int(f_hemax-f_hemin);

  string namefine=name;
  namefine.append("_fine");

  TH1D* henergy_fine = new TH1D(Form("h_%s",namefine.c_str()),
			   Form("%s",namefine.c_str()),
			   bins, f_hemin, f_hemax);

  int allbins=int(7500.-570.);

  string nameall=name;
  nameall.append("_all");

  TH1D* henergy_all = new TH1D(Form("h_%s",nameall.c_str()),
			   Form("%s",nameall.c_str()),
			   allbins, 570., 7500.);

  const int firstBEGE=9;

  // loop over all events
  for (int jentry=0; jentry<nentries; jentry++) 
    {
      //if (jentry%100000==0) cout<<" now event "<<jentry<<endl;
      f_chain->GetEntry(jentry);
      if (firedChannels != 1) continue; 

      if(firedFlag[firstBEGE]==1)
	{
	  henergy_all->Fill(det_energy[firstBEGE]);

	  if(det_energy[firstBEGE]>=2019. && det_energy[firstBEGE]<=2059.)
	    continue;
	  else
	    {
	      henergy->Fill(det_energy[firstBEGE]);
	      henergy_fine->Fill(det_energy[firstBEGE]);
	    }
	}
    }
  

  double scaling = henergy->Integral();
  henergy->Scale(1./scaling);
  henergy_fine->Scale(1./scaling);
  henergy_all->Scale(1./scaling);
  f_MC.push_back(henergy);
  f_MCfine.push_back(henergy_fine);
  f_MCall.push_back(henergy_all);
  f_MCname.push_back(name);

  return 0;
}

// ---------------------------------------------------------

int BEGE_backgrounds::FillMCArrays()
{

  for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
    {
      for(int ibin=0; ibin<f_hnumbins; ibin++)
	{
	  double value=f_MC.at(iMC)->GetBinContent(ibin+1);
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
  
  double logprob = 0.;

//       index = 
// 	GetParameter(Form("ExpPar1_%d",idet))->GetIndex();
//       double ExpPar1 = parameters.at(index);

//       index = 
// 	GetParameter(Form("ExpPar2_%d",idet))->GetIndex();
//       double ExpPar2 = parameters.at(index);

//       index = 
// 	GetParameter(Form("Constant_%d",idet))->GetIndex();
//       double Constant = parameters.at(index);

  for(int ibin=0; ibin<f_hnumbins; ibin++)
    {
      double lambda=0.;

      for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
	{
	  lambda+=parameters.at(iMC)*f_vMC[iMC*f_hnumbins+ibin];
	}	  
// 	  double lowerlimit = f_lowerlimits[ibin];
// 	  double upperlimit = f_upperlimits[ibin];

// 	  double alpha = ((ExpPar1/ExpPar2) 
// 			 * (TMath::Exp(ExpPar2*upperlimit)
// 			    -TMath::Exp(ExpPar2*lowerlimit)))
// 	    + Constant * (upperlimit-lowerlimit);
// 	  double alpha = Constant * (upperlimit-lowerlimit);
	    
      double bincontent = f_vdata[ibin];
      
      double sum = -lambda 
	+ bincontent*log(lambda) 
	- BCMath::LogFact((int)bincontent);
      
      logprob += sum;
    }

  return logprob;
}

// ---------------------------------------------------------

double BEGE_backgrounds::LogAPrioriProbability(std::vector <double> parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

  double logprob = 0.;

  for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
    {
      double range = GetParameter(iMC)->GetRangeWidth();
      // flat prior for all contributions
      logprob += log(1./range);
    }
//       // exponential function + 
//       // constant function to account 
//       // for alpha contribution
//       // flat priors:

//       index = 
// 	GetParameter(Form("ExpPar1_%d",idet))->GetIndex();
//       double d_ExpPar1 = GetParameter(index)->GetRangeWidth();
//       //flat prior for ExpPar1
//       logprob += log(1./d_ExpPar1);

//       index = 
// 	GetParameter(Form("ExpPar2_%d",idet))->GetIndex();
//       double d_ExpPar2 = GetParameter(index)->GetRangeWidth();
//       //flat prior for ExpPar2
//       logprob += log(1./d_ExpPar2);

//       index = 
// 	GetParameter(Form("Constant_%d",idet))->GetIndex();
//       double d_Constant = GetParameter(index)->GetRangeWidth();
//       //flat prior for Constant
//       logprob += log(1./d_Constant);

  
   return logprob;
}
// ---------------------------------------------------------

double BEGE_backgrounds::EstimatePValue()
{
  //Allen's routine for the evaluation of p-value
  //This is derived from PRD 83 (2011) 012004, appendix
  //taken from Luciano
  double logp0=LogLikelihood(GetBestFitParameters()); //likelihood at the mode

  /*
    Now we initialize the Markov Chain by setting the number of entries in a bin to
    the integer 
    part of the mean in that bin.  This will give the maximum possible likelihood.
  */
  double sumlog=0;


  /* mean is the array where you have the expected mean in each bin (calculated
     from the parameter values at the mode. Nom is the nominal value of entries in
     each bin (integer part)
  */

  vector<double> mean(f_hnumbins,0);
  vector<int> nom(f_hnumbins,0);

  vector<double> parameters = GetBestFitParameters();


  for(int ibin=0; ibin<f_hnumbins; ibin++)
    {
      double lambda=0.;

      for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
	{
	  lambda+=parameters.at(iMC)*f_vMC[iMC*f_hnumbins+ibin];
	}	  
// 	  double lowerlimit = f_lowerlimits[ibin];
// 	  double upperlimit = f_upperlimits[ibin];

// 	  double alpha = ((ExpPar1/ExpPar2) 
// 			 * (TMath::Exp(ExpPar2*upperlimit)
// 			    -TMath::Exp(ExpPar2*lowerlimit)))
// 	    + Constant * (upperlimit-lowerlimit);
// 	  double alpha = Constant * (upperlimit-lowerlimit);
	    

      int counter = ibin;

      mean.at(counter) = std::max(lambda,1e-8);
      nom.at(counter) = int(mean.at(counter));
      sumlog += BCMath::LogPoisson(nom[counter],mean[counter]);

    }

  cout << "Logprob for best: " << sumlog << endl;

  /*
  Now we run the Markov chain to generate new configurations.  One iteration means
  a loop over all the bins, with an attempt to vary each bin up or down by one unit.  We
  accept/reject at each step  and compare the data logprob to the simulated at the end of each iteration.
  */
  const int nloops = 100000;
  int Pgood = 0;
  for (int iloop=0;iloop<nloops;iloop++)
    {
      for (int ibin=0;ibin<f_hnumbins;ibin++)
	{
	  int counter = ibin;
	  if (rand()>RAND_MAX/2) // Try to increase the bin content by 1 
	    {
	      double r=mean[counter]/(nom[counter]+1);
	      double rtest=double(rand())/RAND_MAX;
	      if (rtest<r) //Accept
		{
		  nom[counter]=nom[counter]+1;
		  sumlog+=log(r);
		}
	    }
	  else // Try to decrease the bin content by 1 
	    {
	      double r=nom[counter]/mean[counter];
	      double rtest=double(rand())/RAND_MAX;
	      if (rtest<r) //Accept
		{
		  nom[counter]=nom[counter]-1;
		  sumlog+=log(r);
		}
	    }
	}
      if (sumlog<logp0) Pgood++;
    }
  double pvalue = double(Pgood)/double(nloops);
  cout << "p-value is " << pvalue << endl;
  return pvalue;
}

// ---------------------------------------------------------

void BEGE_backgrounds::DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename)
{
  TFile* rootOut=new TFile(rootfilename,"recreate");
  if(!rootOut->IsOpen())
    cout<<"No rootfile opened!"<<endl;
  rootOut->cd();

  TH1D* hMC = new TH1D("hMC",
		       "model",
		       f_hnumbins, f_hemin, f_hemax);

  int bins=int(f_hemax-f_hemin);
  TH1D* hMC_fine = new TH1D("hMC_fine",
			    "model",
			    bins, f_hemin, f_hemax);

  int binsall=int(7500.-570.);
  TH1D* hMC_all = new TH1D("hMC_all",
			   "model",
			   binsall, 570., 7500.);

  TH1D* hresiduals = new TH1D("hresiduals",
			      "residuals",
			      f_hnumbins, f_hemin, f_hemax);


  vector<double> eventsMC;
  vector<double> eventsMC_all;
  //   vector<double> eventsMC_fine;

  for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
    {
      // prepare the histograms
      f_MC.at(iMC)->Scale(parameters.at(iMC));
      f_MCfine.at(iMC)->Scale(parameters.at(iMC));
      f_MCall.at(iMC)->Scale(parameters.at(iMC));
      //f_MC.at(iMC)->Write();
      f_MCfine.at(iMC)->Write();
      f_MCall.at(iMC)->Write();
      hMC->Add(f_MC.at(iMC));
      hMC_fine->Add(f_MCfine.at(iMC));
      hMC_all->Add(f_MCall.at(iMC));

      // count the events
      eventsMC.push_back(f_MC.at(iMC)->Integral());
      eventsMC_all.push_back(f_MCall.at(iMC)->Integral());
      
      //       int lowerBin=f_MCfine.at(iMC)->FindBin(1000.);
      //       int higherBin=f_MCfine.at(iMC)->FindBin(1450.);
      
      //       eventsMC_fine.push_back(f_MCfine.at(iMC)->Integral(lowerBin, higherBin-1));
    }


  double binwidth=f_hdata.at(0)->GetBinWidth(1);

  // write the single detector data spectra
  for(int idet=0; idet<f_ndets; idet++)
    {
      f_hdata.at(idet)->SetLineWidth(2);
      f_hdata.at(idet)->GetXaxis()->SetTitle("Energy (keV)");
      f_hdata.at(idet)->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth));
      f_hdata.at(idet)->Write();
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

  f_hdataSum->SetLineWidth(2);
  f_hdataSum->SetMarkerStyle(20);
  f_hdataSum->GetXaxis()->SetTitle("Energy (keV)");
  f_hdataSum->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth));
  f_hdataSum->Draw("P");
  f_hdataSum->Write();
  legend->AddEntry(f_hdataSum,"data","p");

  hMC->SetLineWidth(2);
  hMC->SetLineColor(1);
  hMC->Draw("same");
  hMC->GetXaxis()->SetTitle("Energy (keV)");
  hMC->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth));
  hMC->Write();
  legend->AddEntry(hMC,"model","l");

  //the single contributions

  string name;
  
  for(int iMC=0; iMC<(int)f_MC.size(); iMC++)
    {
      f_MC.at(iMC)->SetLineWidth(2);
      if(iMC<10)
	{
	  f_MC.at(iMC)->SetLineColor(iMC+2);
	}
      if(iMC>=8)
	{
	  f_MC.at(iMC)->SetLineColor(iMC-8+2);
	  f_MC.at(iMC)->SetLineStyle(2);
	}
      f_MC.at(iMC)->Draw("same");
      f_MC.at(iMC)->GetXaxis()->SetTitle("Energy (keV)");
      f_MC.at(iMC)->GetYaxis()->SetTitle(Form("Events/(%d keV)",(int)binwidth));
      name=f_MCname.at(iMC);
      legend->AddEntry(f_MC.at(iMC),Form("%s",name.c_str()),"l");

      f_MC.at(iMC)->Write();
    }

  legend->Draw();
  canvas->Write();
      
  hresiduals->Add(f_hdataSum, hMC, 1., -1.);
  hresiduals->SetLineWidth(2);
  hresiduals->SetMarkerStyle(20);
  hresiduals->GetXaxis()->SetTitle("Energy (keV)");
  hresiduals->GetYaxis()->SetTitle(Form("residual counts/(%d keV)",(int)binwidth));
  hresiduals->Write();

  
  // events info
  cout<<"---------------------------------------------"<<endl;
  cout<<"---------------------------------------------"<<endl<<endl;
  for(int iMC=0; iMC<(int)eventsMC.size(); iMC++)
    {
      cout<<"Events: "<<endl;
      cout<<f_MCname.at(iMC)<<": "<<eventsMC.at(iMC)<<endl;
      cout<<"In total spectrum from (570 - 7500) keV"<<endl;
      cout<<f_MCname.at(iMC)<<": "<<eventsMC_all.at(iMC)<<endl;
    }
  cout<<"Total events in MC: "<<hMC->Integral()<<endl;
  cout<<endl;
  cout<<"Events in data: "<<f_hdataSum->Integral()<<endl;
  cout<<endl;
  cout<<"---------------------------------------------"<<endl;
  cout<<"---------------------------------------------"<<endl<<endl;

//   cout<<"---------------------------------------------"<<endl;
//   cout<<"---------------------------------------------"<<endl<<endl;
//   double sum=0.;
//   for(int iMC=0; iMC<(int)eventsMC_fine.size(); iMC++)
//     {
//       cout<<"Events between 1 and 1.45 MeV: "<<endl;
//       cout<<f_MCname.at(iMC)<<": "<<eventsMC_fine.at(iMC)<<endl;
//       sum+=eventsMC_fine.at(iMC);
//     }
//   cout<<"In total "<<sum<<" events in this range"<<endl;
//   cout<<endl;
//   cout<<"---------------------------------------------"<<endl;
//   cout<<"---------------------------------------------"<<endl<<endl;

  rootOut->Close();
}
// ---------------------------------------------------------

