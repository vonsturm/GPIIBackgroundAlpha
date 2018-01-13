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
      GPIIBackgroundAlpha(std::string masterconfname);
      ~GPIIBackgroundAlpha();

      // Methods to overload
      void DefineParameters(); // Priors are defined here
//      double LogAPrioriProbability(const std::vector<double> &parameters); // ONLY NEEDED IF PRIORS ARE MORE COMPLICATED
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
      int ReadRunData( std::string keylist );
      int FillDataArray();

      // Read MC pdfs
      int InitializeMCHistograms();
      int ReadMC();
      int ReadSingleMC( int par_index, int local_index, int global_index, std::string histoname, std::string filename );
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

      // Output directory and filename
      std::string GetOutputDirectory();
      std::string GetOutputFilenameBase();
      std::string GetRunConfLTName();
      std::string GetDetConfLTName();

      bool SkipParameter( int p );
      std::string ParameterName( int p );

      bool GetWriteMCMCChain(){ return f_j_masterconf["MCMC-write-chain"].asBool(); };

      void DumpHistosAndInfo( std::string rootfilename );

      void UpdateParameters();

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

     unsigned int f_ndets;
     unsigned int f_nruns;
     unsigned int f_npars;

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
