// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__GPIIBACKGROUNDALPHA__H
#define __BAT__GPIIBACKGROUNDALPHA__H

// C/C++ includes
#include <vector>
#include <map>

// jsoncpp
#include "json/value.h"

// ROOT includes
#include "TH1D.h"
#include "TChain.h"

// BAT includes
#include <BAT/BCModel.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>

// gerda ada
#include "GETRunConfiguration.hh"

// This is a GPIIBackgroundAlpha header file.
// Model source code is located in file GPIIBackgroundAlpha/GPIIBackgroundAlpha.cxx

// ---------------------------------------------------------
class GPIIBackgroundAlpha : public BCModel
{
   public:

      // Constructors and destructor
      GPIIBackgroundAlpha();
      GPIIBackgroundAlpha(const char * name);
      ~GPIIBackgroundAlpha();

      // Methods to overload, see file GPIIBackgroundAlpha.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector <double> & parameters);

      // Histograms
      void SetHistogramParameters(int hNbins, double hMin, double hMax)
      {
          f_hnumbins = hNbins;
          f_hemin = hMin;
          f_hemax = hMax;
      }
      int InitializeDataHistograms( std::vector<std::string> detectorlist );

      // Read Data
      int ReadData( std::string runlist, std::string data_set,
          std::string detectorlistname, bool useDetectorList );
      int ReadRunData( std::string keylist, std::vector<std::string> detectorlist );
      int FillDataArray();

      // Read MC pdfs
      void SetParConfigFile( std::string name );
      std::string GetParConfigFile(){ return f_parConfigFile; };

      int InitializeMCHistograms();
      int ReadMC();
      int AddMC( std::string name );
      int AddMCSingle( std::string name, std::string histoname );
      int FillMCArrays();

      std::vector<std::string> GetMCParNames();

      double getndets() {return f_ndets;};

      double EstimatePValue();
      void DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename);

      Json::Value GetJsonValueFromFile( std::string filename );

      bool IsOn( GETRunConfiguration * RunConf, std::string det );
      int GetChannel( GETRunConfiguration * RunConf, std::string det );

      void SetVerbosity( int v ) { f_verbosity = v; return; };
      int GetVerbosity() { return f_verbosity; };

 private:

      int f_ndets;
      int f_hnumbins;
      double f_hemin;
      double f_hemax;

      std::string f_parConfigFile;
      Json::Value f_j_parameters;

      std::vector<double> f_RunLiveTime;
      std::map<std::string,double> f_DetectorLiveTime;
      //std::vector<int> fDetectorDynamicRange;

      std::map<std::string,TH1D*> f_hdata;
      TH1D* f_hdataSum;
      TH1D* f_hdataSum_fine;
      TH1D* f_hdataSum_all;

      std::string f_MC_FileName;

      std::vector<TH1D*> f_MC;
      std::vector<TH1D*> f_MC_fine;
      std::vector<TH1D*> f_MC_all;

      std::vector<double> f_vdata;
      std::vector<double> f_vMC;

      std::vector<double> f_lowerlimits;
      std::vector<double> f_upperlimits;

      int f_verbosity;

};
// ---------------------------------------------------------

#endif
