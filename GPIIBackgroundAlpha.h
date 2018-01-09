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

// ROOT includes
#include "TH1D.h"
#include "TChain.h"

// BAT includes
#include <BAT/BCModel.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>


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
      int InitializeHistograms( std::vector<std::string> detectorlist );

      // Read Data
      int ReadData( std::string runlist, std::string data_set,
          std::string detectorlistname, bool useDetectorList, int verbosity = 0 );
      int ReadRunData( std::string keylist, std::vector<std::string> detectorlist );
      int FillDataArray();

      // Read MC pdfs
      void SetParConfigFile( std::string name ){ f_parConfigFile = name; return; };
      std::string GetParConfigFile(){ return f_parConfigFile; };

      int ReadMC();
      int AddMC( std::string name );
      int AddMCSingle( std::string name, std::string histoname );
      int FillMCArrays();

      std::vector<std::string> GetMCParNames(){ return f_MCname; };

      double getndets() {return f_ndets;};

      double EstimatePValue();
      void DumpHistosAndInfo(std::vector<double> parameters, char* rootfilename);

 private:

      int f_ndets;
      int f_hnumbins;
      double f_hemin;
      double f_hemax;

      std::string f_parConfigFile;

      std::vector<double> f_RunLiveTime;
      std::map<std::string,double> f_DetectorLiveTime;
      //std::vector<int> fDetectorDynamicRange;

      std::map<std::string,TH1D*> f_hdata;
      TH1D* f_hdataSum;
      TH1D* f_hdataSum_fine;

      std::string f_MC_FileName;

      std::vector<TH1D*> f_MC;
      std::vector<TH1D*> f_MCfine;
      std::vector<TH1D*> f_MCall;
      std::vector<std::string> f_MCname;

      std::vector<double> f_vdata;
      std::vector<double> f_vMC;

      std::vector<double> f_lowerlimits;
      std::vector<double> f_upperlimits;

};
// ---------------------------------------------------------

#endif
