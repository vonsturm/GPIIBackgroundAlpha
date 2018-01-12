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
#include <BAT/BCLog.h>

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
      GPIIBackgroundAlpha(string masterconfname);
      ~GPIIBackgroundAlpha();

      // Methods to overload
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector <double> & parameters);
      double EstimatePValue();

      // Histograms
      void SetHistogramParameters(double hBinning, double hMin, double hMax);
      void SetHistogramParameters(int hNbins, double hMin, double hMax)
      {
          f_hnumbins = hNbins;
          f_hemin = hMin;
          f_hemax = hMax;
      };

      // Read Data
      int InitializeDataHistograms();
      int ReadData(); //WRAPPER
      int ReadDataFromHistogram( std::string infilename );
      int ReadDataFromEvents( std::string outfilename );
      int WriteDataToFileForFastAccess( std::string outfilename );
      int ReadRunData( std::string keylist, std::vector<std::string> detectorlist );
      int FillDataArray();

      // Read MC pdfs
      int InitializeMCHistograms();
      int ReadMC();
      int AddMC( std::string name );
      int AddMCSingle( std::string name, std::string histoname );
      int FillMCArrays();

      //---- GETTERS AND SETTERS ----
      void SetVerbosity( int v ) { f_verbosity = v; return; };
      int GetVerbosity() { return f_verbosity; };
      //
      void SetPrecision( std::string precisionString );
      //
      void SetNdetectors( int n ){ return; };
      int GetNdetectors(){ return f_ndets; };
      //
      std::vector<std::string> GetMCParNames();

      //---- NICE TO HAVE ----
      Json::Value GetJsonValueFromFile( std::string filename );

      bool IsOn( GETRunConfiguration * RunConf, std::string det );
      int GetChannel( GETRunConfiguration * RunConf, std::string det );

      // FIX ME
      void DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename);

 private:

     //----- METHODS -----

     // Set json config file
     void SetMasterConf( std::string masterconfname ); // wrapps all of them
     void SetDetConf( std::string detconfname )
     {
         f_j_detconf = GetJsonValueFromFile( detconfname );
         f_ndets = f_j_detconf["detectors"].size();
         return;
     };
     void SetRunConf( std::string runconfname )
     {
         f_j_runconf = GetJsonValueFromFile( runconfname );
         f_nruns = f_j_runconf["runs"].size();
         return;
     };
     void SetParConf( std::string parconfname )
     {
         f_j_parconf = GetJsonValueFromFile( parconfname );
         f_npars = f_j_parconf["parameters"].size();
         return;
     };
     void SetEnvConf( std::string envconfname ){ f_j_envconf = GetJsonValueFromFile( envconfname ); return; };

     void UnwrapMasterConf();

     //----- PARAMETERS -----

     std::string f_masterconfname;

     Json::Value f_j_masterconf;
     Json::Value f_j_detconf;
     Json::Value f_j_runconf;
     Json::Value f_j_parconf;
     Json::Value f_j_envconf;

     int f_ndets;
     int f_nruns;
     int f_npars;

     int f_hnumbins;
     double f_hemin;
     double f_hemax;

     std::vector<double> f_RunLiveTime;
     std::map<std::string,double> f_DetectorLiveTime;
     //std::vector<int> fDetectorDynamicRange;

     std::map<std::string,TH1D*> f_hdata;
     TH1D* f_hdataSum;
     TH1D* f_hdataSum_fine;
     TH1D* f_hdataSum_all;

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
